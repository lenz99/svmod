#!/usr/bin/env Rscript
# mkuhn, 20140624, 20150304
#
# Example script how to use one-class SVM to learn what is a novelty from simulated mapping data.
# Using my own tuning and evaluation code (and not the mlr code)
#
# INPut: patient data and feature data
# OUTput: different evaluation metrics (saved as RDS)

library(logging); basicConfig()
library(dplyr)
library(e1071)
library(scatterplot3d)

library(devtools)
dev_mode(on=TRUE)
svmod_pck <- require(svmod)

stopifnot( svmod_pck )


baseDir <- svmod::BASEDIR
myNGSProp <- ngs.prop()
mySampleProp <- sample.prop()


patDataFile <- paste0(baseDir, "patient/patData.rds")
patData <- readRDS(file=patDataFile)
loginfo("Loaded patient data ::%s:: from %d patients with in total %d regions.", patDataFile, length(unique(patData$patId)), NROW(patData))
# add information to patient data
#patData$patIdNbr <- getPatNbr(patData$patId)

patData <- patData[order(patData$patId, patData$targetStartPos, patData$targetEndPos), ]


## FEATURE DATA
mappingDir <- getVirtualPatientPath(baseDir, what="mapping", .sample.prop = mySampleProp, .ngs.prop = myNGSProp)
loginfo("Using feature file from %s.", mappingDir)
featMat0 <- readRDS(file=paste0(mappingDir, "features.rds"))  #"featuresPlusY.rds"
featDat0 <- data.frame(featMat0) %>% arrange_(~patId, ~targetStartPos, ~targetEndPos)
loginfo("Found %d feature entries in ::%s::. Every simulated patient position should correspond to a feature entry.", NROW(featDat0), mappingDir)

# add simulation info
featDat <- patData %>%
  dplyr::select(patId, chrom, targetStartPos, targetEndPos, status,SVtype) %>% 
  dplyr::inner_join(y = featDat0,
                    by = c("patId", "targetStartPos", "targetEndPos"))


# ## FEATURE FIXING
# # dichotimize feature "region.length"
# featureDat[, "feat.region.length"] <- ifelse(featureDat[, "feat.region.length"]<801, yes=0, no=1) #short: 0, normal: 1
# # region to target: on border:0, within 400: 1, further away: 2
# featureDat[, "feat.region.targetDist.5p"] <- ifelse(featureDat[, "feat.region.targetDist.5p"] == 0, yes=0,
#                                                     ifelse(featureDat[, "feat.region.targetDist.5p"] <= 400, yes=1, no=2))
# featureDat[, "feat.region.targetDist.3p"] <- ifelse(featureDat[, "feat.region.targetDist.3p"] == 0, yes=0,
#                                                     ifelse(featureDat[, "feat.region.targetDist.3p"] <= 400, yes=1, no=2))


# feature selection
# ad-hoc feature selection based on correlation!
featDat %<>%
  dplyr::select(patId:endPos, ends_with(".d"), ends_with(".t"), -feat.ins.q1.prop.d, -feat.del.q1.prop.d, -feat.clipped.q1.prop.d, -feat.cov.max.d, -feat.cov.sd.max.d, -feat.cov.ratio.d) %>% 
  dplyr::filter_(~complete.cases(.))

# drop constant features from feature set
featCols.const <- apply(featDat[,-c(1:6)], 2, sd, na.rm=TRUE) == 0
stopifnot( ! any(featCols.const) )
# if ( any(featCols.const) ){
#   loginfo("Drop features %s. They have no variance.", names(which(featCols.const)))
#   featDat <- featDat[, ! 6+featCols.const] # shift 6 first columns
# }

loginfo("Feature selection and complete data filter yield %d features and %d entries.", NCOL(featDat), NROW(featDat))

# result variable
featDat.status <- ifelse(featDat$status == 'germline', 0L, 1L) # 0='germline', 1='differential SV'
featDat.SVtype <- featDat$SVtype

# remove patient info
featDat <- dplyr::select(featDat, -c(patId:endPos))

ssi <- subsampleInd(dat.status = (featDat.status == 0L), prop.outl = 0.1)
                                                    

# train / validation sets
set.seed(123456)
TRAIN.PROP <- .7
loginfo("Split data into training (%.1f%%) and validation data.", TRAIN.PROP*100)

trainInd <- sort(sample(NROW(featDat), ceiling(NROW(featDat)*TRAIN.PROP)))
featDat.train <- featDat[trainInd,]
featDat.val   <- featDat[-trainInd,]


SVtype.train <- featDat.SVtype[trainInd]
SVtype.val <- featDat.SVtype[-trainInd]
status.train <- featDat.status[trainInd]
status.val <- featDat.status[-trainInd]


# as matrix
featMat.train <- as.matrix(featDat.train)
featMat.val <- as.matrix(featDat.val)



# learning algorithm: one-class SVM
loginfo("Learn algorithm on training data: novelty detection one-class SVM")


loginfo("Example: some one-class SVM with linear kernel (=no kernel)")
svmMod.Fdiff.lin1 <- applyOneClassSVM(X=featMat.train, nu=0.1, kernel="linear")
evaluateOneClassSVM(svmMod.Fdiff.lin1, X=featMat.train, y=status.train, do.plot="2d")
evaluateOneClassSVM(svmMod.Fdiff.lin1, newdata=featMat.val, y=status.val, do.plot="2d")

loginfo("Example: some one-class SVM with polynomial kernel")
svmMod.Fdiff.poly1 <- applyOneClassSVM(X=featMat.train, nu=0.1, kernel="polynomial", degree=3, gamma=0.1)
evaluateOneClassSVM(svmMod.Fdiff.poly1, X=featMat.train, y=status.train, do.plot="no")

loginfo("Example: some one-class SVM with radial kernel")
svmMod.Fdiff.rad1 <- applyOneClassSVM(X=featMat.train, nu = 0.05, kernel = "radial", gamma = 0.1)
evaluateOneClassSVM(svmMod.Fdiff.rad1, X = featMat.train, y=status.train, do.plot = "3d")
eval.asvm.rad1 <- evaluateOneClassSVM(svmMod.Fdiff.rad1, newdata = featMat.val, y=status.val, do.plot = "3d")



loginfo("Tune hyper-parameters of one-class SVM with radial kernel towards accuracy (ACC)")
tuneObj.acc <- ocSVMTune(train.x = featMat.train, train.status = !as.logical(status.train),
                     ranges=list(nu=c(0.05, 0.1, 0.5), gamma = 2^(-3:-2)))
tuneObj.acc
summary(tuneObj.acc)
plot(tuneObj.acc)

saveRDS(tuneObj.acc, "~/seq1/SVdata/pres/data/tuneACC.rds")
# evaluation on separate validation data to be meaningful in error measures
eval.svm.acc <- evaluateOneClassSVM(tuneObj.acc$best.model, tunedFor = "Acc", newdata = featMat.val, y=status.val, do.plot="3d")


loginfo("Tune hyper-parameter of one-class SVM with radial kernel towards Recall (SEN)")
tuneObj.recall <- ocSVMTune(train.x = featMat.train, train.status = !as.logical(status.train),
                         ranges=list(nu=c(0.05, 0.15), gamma = 2^(-5:-2)), tune.measure="recall")
tuneObj.recall
plot(tuneObj.recall)
summary(tuneObj.recall)

saveRDS(tuneObj.recall, "~/seq1/SVdata/pres/data/tuneSEN.rds")

# evaluation on separate validation data to be generalizable
eval.svm.recall <- evaluateOneClassSVM(tuneObj.recall$best.model, tunedFor = "Recall", newdata = featMat.val, y=status.val, do.plot="3d")


loginfo("Tune hyper-parameter of one-class SVM with radial kernel towards F1-score (F1)")
tuneObj.f1 <- ocSVMTune(train.x = featureDat.diff.train, train.status = !as.logical(isDSV.train),
                            ranges=list(nu=c(0.05, 0.15), gamma = 2^(-5:-2)), tune.measure="f1")
saveRDS(tuneObj.f1, "~/seq1/SVdata/pres/data/tuneF1.rds")
tuneObj.f1
plot(tuneObj.f1)
summary(tuneObj.f1)
# 0 error can occur here (if CV-validation is on too small folds!?)

# evaluation on separate validation data to be generalizable
eval.svm.f1 <- evaluateOneClassSVM(tuneObj.f1$best.model, newdata=featureDat.diff.val, y=isDSV.val, do.plot="3d")


loginfo("Use minimal SVM, focus on clipped reads")
svmMod.clip.rad1 <- applyOneClassSVM(X=featureDat.diff.train[, c("feat.clipped.q1.bases.d", "feat.cov.sd.d")],
                                     nu = 0.05, kernel = "radial", gamma = 0.1)


selectedFeatures <- c("feat.clipped.q1.bases.d", "feat.cov.sd.d", "feat.pe.isize.sd.d", "feat.cov.roll.maxAvg.d")
tuneObj_sel.f1 <- ocSVMTune(train.x = featureDat.diff.train[, selectedFeatures],
                             train.status = !as.logical(isDSV.train),
                        ranges=list(nu=c(0.05, 0.15), gamma = 2^(-5:-2)), tune.measure="f1")
saveRDS(tuneObj_sel.f1, "~/seq1/SVdata/pres/data/tune_selF1.rds")
tuneObj_sel.f1
plot(tuneObj_sel.f1)
summary(tuneObj_sel.f1)
eval.svm_sel.f1 <- evaluateOneClassSVM(tuneObj_sel.f1$best.model, newdata=featureDat.diff.val[, selectedFeatures], y=isDSV.val,
                                       do.plot='3d', tunedFor="f1")
