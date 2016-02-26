#!/usr/bin/env Rscript
# mkuhn, 20140624
#
# Example script how to use one-class SVM to learn what is a novelty from simulated mapping data.


library(logging); basicConfig()
library(e1071)
library(scatterplot3d)

library(devtools)
dev_mode(on=TRUE)

library(svmod)


baseDir <- svmod::BASEDIR
myNGSProp <- ngs.prop()
mySampleProp <- sample.prop()


patDataFile <- paste0(baseDir, "patient/patData.rds")
patData <- readRDS(file=patDataFile)
loginfo("Loaded patient data ::%s:: from %d patients with in total %d regions.", patDataFile, length(unique(patData$patId)), NROW(patData))
# add information to patient data
#patData$patIdNbr <- getPatNbr(patData$patId)
patData$isDSV <- ifelse(patData$SVtype != 'no', 1, 0)

patData <- patData[order(patData$patId, patData$startPos, patData$endPos), ]
isDSV <- patData$isDSV

## FEATURE DATA
mappingDir <- getVirtualPatientPath(baseDir, what="mapping", sample.prop = mySampleProp, ngs.prop = myNGSProp) #paste0(baseDir, "mapping/", toString(mySampleProp), "/", toString(myNGSProp), "/")
loginfo("Using feature file from %s.", mappingDir)
featureDat <- readRDS(file=paste0(mappingDir, "features.rds"))  #"featuresPlusY.rds"
featureDat <- featureDat[order(featureDat[, "patId"], featureDat[, "startPos"], featureDat[, "endPos"]),]
loginfo("Found %d feature entries in ::%s::. Every simulated patient position should correspond to a feature entry.", NROW(featureDat), mappingDir)
# stopifnot( colnames(featureDat)[NCOL(featureDat)] == 'isDSV' )
# stopifnot( all( featureDat[, 'isDSV'] %in% c(0,1)) )
# 
# isDSV <- featureDat[, 'isDSV']
# featureDat <- featureDat[, -NCOL(featureDat)]

# # sort patData in same order as featureDat
# patData2 <- merge(patData, y=featureDat[,c("patId", "startPos", "endPos")], by.x=c("patIdNbr", "startPos", "endPos"), by.y=c("patId", "startPos", "endPos"),
#                   sort = TRUE)



# drop patID columns.
# Ordering ensures that patData and featureData files are aligned.
# ZZZ think of targetDist columns, & col5: feat.region.length. Dichotimize?
featureDat <- featureDat[, -c(1:5)]

## FEATURE FIXING
# dichotimize feature "region.length"
featureDat[, "feat.region.length"] <- ifelse(featureDat[, "feat.region.length"]<801, yes=0, no=1) #short: 0, normal: 1
# region to target: on border:0, within 400: 1, further away: 2
featureDat[, "feat.region.targetDist.5p"] <- ifelse(featureDat[, "feat.region.targetDist.5p"] == 0, yes=0,
                                                    ifelse(featureDat[, "feat.region.targetDist.5p"] <= 400, yes=1, no=2))
featureDat[, "feat.region.targetDist.3p"] <- ifelse(featureDat[, "feat.region.targetDist.3p"] == 0, yes=0,
                                                    ifelse(featureDat[, "feat.region.targetDist.3p"] <= 400, yes=1, no=2))
                                                    
# drop constant features from feature matrix
featureCols.const <- apply(featureDat, 2, sd, na.rm=TRUE) == 0
if ( any(featureCols.const) ){
  loginfo("Drop features %s. They have no variance.", names(which(featureCols.const)))
  featureDat <- featureDat[, ! featureCols.const]
}

# drop the baseline .1-features
featureDat.diff <- featureDat[, ! grepl("[.]1$", colnames(featureDat))] #, ends_with(".d"))


# train / validation sets
set.seed(123456)
TRAIN.PROP <- .7
loginfo("Split data into training (%.1f%%) and validation data.", TRAIN.PROP*100)

trainInd <- sort(sample(NROW(featureDat), ceiling(NROW(featureDat)*TRAIN.PROP)))
featureDat.train <- featureDat[trainInd,]
featureDat.val   <- featureDat[-trainInd,]
patDat.train <- patData[trainInd,]
patDat.val <- patData[-trainInd,]

featureDat.diff.train <- featureDat.diff[trainInd, ]
featureDat.diff.val   <- featureDat.diff[-trainInd, ]

isDSV.train <- isDSV[trainInd]
isDSV.val   <- isDSV[-trainInd]



# learning algorithm: one-class SVM
loginfo("Learn algorithm on training data: novelty detection one-class SVM")


loginfo("Example: some one-class SVM with linear kernel (=no kernel)")
svmMod.Fdiff.lin1 <- applyOneClassSVM(X=featureDat.diff.train, nu=0.1, kernel="linear")
evaluateOneClassSVM(svmMod.Fdiff.lin1, X=featureDat.diff.train, y=isDSV.train, do.plot="2d")
evaluateOneClassSVM(svmMod.Fdiff.lin1, newdata=featureDat.diff.val, y=isDSV.val, do.plot="2d")

loginfo("Example: some one-class SVM with polynomial kernel")
svmMod.Fdiff.poly1 <- applyOneClassSVM(X=featureDat.diff.train, nu=0.1, kernel="polynomial", degree=3, gamma=0.1)
evaluateOneClassSVM(svmMod.Fdiff.poly1, X=featureDat.diff.train, y=isDSV.train, do.plot="no")

loginfo("Example: some one-class SVM with radial kernel")
svmMod.Fdiff.rad1 <- applyOneClassSVM(X=featureDat.diff.train, nu = 0.05, kernel = "radial", gamma = 0.1)
evaluateOneClassSVM(svmMod.Fdiff.rad1, X = featureDat.diff.train, y=isDSV.train, do.plot = "3d")
eval.asvm.rad1 <- evaluateOneClassSVM(svmMod.Fdiff.rad1, newdata = featureDat.diff.val, y=isDSV.val, do.plot = "3d")



loginfo("Tune hyper-parameters of one-class SVM with radial kernel towards accuracy (ACC)")
tuneObj.acc <- ocSVMTune(train.x = featureDat.diff.train, train.status = !as.logical(isDSV.train),
                     ranges=list(nu=c(0.05, 0.15), gamma = 2^(-5:-2)))
tuneObj.acc
summary(tuneObj.acc)
plot(tuneObj.acc)

saveRDS(tuneObj.acc, "~/seq1/SVdata/pres/data/tuneACC.rds")
# evaluation on separate validation data to be meaningful in error measures
eval.svm.acc <- evaluateOneClassSVM(tuneObj.acc$best.model, tunedFor = "Acc", newdata = featureDat.diff.val, y=isDSV.val, do.plot="3d")


loginfo("Tune hyper-parameter of one-class SVM with radial kernel towards Recall (SEN)")
tuneObj.recall <- ocSVMTune(train.x = featureDat.diff.train, train.status = !as.logical(isDSV.train),
                         ranges=list(nu=c(0.05, 0.15), gamma = 2^(-5:-2)), tune.measure="recall")
tuneObj.recall
plot(tuneObj.recall)
summary(tuneObj.recall)

saveRDS(tuneObj.recall, "~/seq1/SVdata/pres/data/tuneSEN.rds")

# evaluation on separate validation data to be generalizable
eval.svm.recall <- evaluateOneClassSVM(tuneObj.recall$best.model, tunedFor = "Recall", newdata = featureDat.diff.val, y=isDSV.val, do.plot="3d")


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
