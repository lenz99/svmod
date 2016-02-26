
#' Own MCC-like measure, implemented as a Measure object of the mlr framework. There was a bug in MCC in mlr-2.4
my.mcc <- mlr::makeMeasure(id = "my.mcc", minimize = FALSE, properties = c("classif"), worst = -1L, best = 1L,
                           note = "modified MCC (cleared denominator if a row or column of the 2x2 confusion table is zero)",
                           fun = function(task, model, pred, feats, extra.args) {
                             #   tb = table(pred$data$response, pred$data$truth)
                             #   1 - sum(diag(tb)) / sum(tb)
                             truth <- pred$data$truth; response <- pred$data$response
                             negative <- pred$task.desc$negative
                             positive <- pred$task.desc$positive
                             
                             tn = mlr::measureTN(truth, response, negative) 
                             tp = mlr::measureTP(truth, response, positive) 
                             fn = mlr::measureFN(truth, response, negative) 
                             fp = mlr::measureFP(truth, response, positive) 
                             
                             #   #denom <- prod(table(truth, response))
                             #   denom <- (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
                             #   
                             #   if ( ! exists("denom") || is.null(denom) || denom <= 0) {
                             #     #cat("my.MCC: Denumerator is <= 0!\n")
                             #     return(0L)
                             #   }
                             #   (tp * tn - fp * fn)/sqrt(denom)
                             
                             # original Matthews formula (wiki)
                             n <- tn + tp + fn + fp
                             s <- (tp + fn) / n
                             p <- (tp + fp) / n
                             
                             if ( s == 0L || s == 1L || p == 0L || p == 1L) return(0L)
                             (tp / n - s * p) / sqrt(p * s * (1-s) * (1-p))
                           })



#' Evaluate some learners on a feature set, based on simulation results.
#' 
#' The feature set from simulation run is specified by the sample- and NGS-properties. It writes out the benchmark object in the corresponding performance/ -directory (if requested).
#' One needs also to give the margin left and right of the simulated SV that was used prior to extract features from the mapping.
#' \describe{
#'   \item{with \code{useSimSVInfo=TRUE}}{the simulated SV-information is directly used, i.e. no SV-screening step has taken place but features have been extracted directly at the simulated SV-locus.}
#'   \item{with \code{useSimSVInfo=FALSE}}{the SV-screening step has been performed that resulted in screening negative regions where SV-injected gave rise to positive and negative training cases.
#'   The screening positive regions form the test set where the learners naturally need to be evaluated.}
#' }
#' 
#' @param sample.prop sample properties underlying the features
#' @param ngs.prop NGS properties underlying the features
#' @param patIds numeric vector of patient IDs to restrict features to only given patients. Defaults to \code{NULL}, i.e. all available features
#' @param k number of sub-sampling iterations for tuning (inner loop). Defaults to 10.
#' @param prop.sSV proportion of somatic SVs in the genomic regions in training data
#' @param margin.bp margin size (in bp) that was used for feature extraction (either flanking the exact simulated SV or flanking the screened regions). The margin is part of the feature filename.
#' @param useSimSVInfo flag. Do we want to use the knowledge of SV-presence and SV-coordinates for the evaluation? Defaults to \code{FALSE}, i.e. use read coverage and clipped base peaks to learn regions where features are extracted.
#' @param saveEvaluation flag if the benchmark is to be written out to disk in the performance/-directory.
#' @param verbose logical flag: show info during tuning/evaluation
#' @return dataframe with benchmark performance results.
#' @export
evalLearners_OLD <- function(baseDir=BASEDIR, sample.prop, ngs.prop, patIds=NULL,
                         #parallelMode=c("socket", "multicore"), nbrCores=2L, parallelLevel=c("mlr.tuneParams", "mlr.resample", "mlr.benchmark", "mlr.selectFeatures"),
                         tuneMode=c("grid", "random", "gensa"),
                         k=10L, prop.sSV=0.05, margin.bp=25L, useSimSVInfo=FALSE, saveEvaluation=TRUE, verbose=FALSE){
  
  # parallelMode <- match.arg(parallelMode)
  # parallelLevel <- match.arg(parallelLevel)
  tuneMode <- match.arg(tuneMode)
  
  # example call
#   for (j in c(60, 40, 20)){
#     cat("COV",j,"\n")
#     for (i in c(10, 25, 33, 50, 75, 90, 95, 99)){
#       cat("COV",j," Tumorload ",i, "\n")
#       sp <- sample.prop(cov=j, tumorLoad=i); print(sp)
#       el <- try(evalLearners(sample.prop=sp, ngs.prop=ngs.prop(), tuneMode="grid", k=5L, verbose=FALSE))
#       #evalLearners(baseDir=BASEDIR, sample.prop=sp, ngs.prop=ngs.prop(), prop.sSV=0.05, margin.bp=25L, saveEvaluation=TRUE)
#     }#rof i
#   }#rof j
  
  stopifnot( is.numeric(margin.bp) )
  
  # for quick start only to set arguments
  if (FALSE){
    baseDir <- BASEDIR
    sample.prop <- sample.prop(cov=60L, tumorLoad=90L)
    ngs.prop <- ngs.prop()
    patIds <- NULL
    tuneMode <- "random"
    #nbrCores=3L; parallelMode <- "socket";  parallelLevel <- "mlr.tuneParams"
    prop.sSV=0.05; margin.bp=25L; useSimSVInfo <- FALSE;
    k <- 5L; saveEvaluation=FALSE; verbose=FALSE
  }#fi
  
  
  
  # Read in data -----
  
  featFileBase <- paste0(if (isTRUE(useSimSVInfo)) "features_sim_margin" else "features_screen_margin", margin.bp, ".rds")  #ZZZ fix naming here!
  featFile <- file.path(getVirtualPatientPath(basedir = baseDir, what = "feature", .sample.prop = sample.prop, .ngs.prop = ngs.prop), featFileBase )
  
  if (! isTRUE(file.exists(featFile)) ){
    logging::logerror("Feature data file %s not found!", featFile)
    return(invisible(NULL))
  }
  
  data <- readRDS(featFile)
  

  
  
  
  # mkuhn, 20150504: optional patient filter on feature dataframe
  if ( ! is.null(patIds) && is.numeric(patIds) && length(patIds) >= 1L ){
    data <- dplyr::filter_(data, ~patId %in% patIds)
    logging::loginfo("Restricted to %d feature rows that belong to given patients.", NROW(data))
  }
  
  # check if patient info is available within feature data
  if ( ! (is.data.frame(data) && NROW(data) >= 1L && "status" %in% names(data) && is.factor(data$status) && nlevels(data$status) == 2L) ){
    logging::logwarn("Something wrong with feature data in %s. Status missing?", featFile)
    return(invisible(NULL))
  }
  
  
  
  
  # Task ----
  
  # prepare data for learner
  dat.task <- data %>% 
    dplyr::select_(~one_of("status", "type"), ~ends_with(".d"), ~ends_with(".t")) %>% 
    dplyr::filter_(~complete.cases(.))
  
  # # mkuhn, 2015-12-28: when we had NA in status in test data
  # if (isTRUE(useSimSVInfo)){
  #   dat.task %<>% dplyr::filter_(~complete.cases(.))
  # } else {
  #   # for SV-injection we have NA-status for all test-data. Hence, complete cases except for status here!
  #   dat.task <- dat.task[complete.cases(dat.task[, names(dat.task) != 'status']),]
  # }
  
  # mkuhn, 2015-04-23: re-name factor levels as the learners expect IN and OUT
  if ( identical(levels(dat.task$status), c("germline", "DSV")) ){
    dat.task$status <- factor(dat.task$status, levels = c("germline", "DSV"), labels = STATUS_LEVELS)
    logging::logwarn("Recode status vector in: ", paste(STATUS_LEVELS, collapse=" & "))
  }
  
  stopifnot( identical(levels(dat.task$status), STATUS_LEVELS) )
  
  
  
  doSubSampling <- FALSE #ZZZ # mkuhn, 2015-12-30 decide here what is best?!
  if ( isTRUE(doSubSampling) ){
    # subsampling: training data with fixed proportion of somatic SVs
    logging::loginfo("Subsampling from available data to get somatic SV proportion prop.sSV=%.3f", prop.sSV)
    #ssInd <- subsampleInd(dat.status = dat.task$status == STATUS_LEVELS[1], prop.outl = prop.sSV)
    ssInd <- subsampleInd(dat.status = dat.task$status, prop.outl = prop.sSV)
    stopifnot( length(ssInd$na.ind) == 0L ) # we should not have NAs any more (mkuhn, 2015-12-29)
    dat.task <- dat.task[c(ssInd$train.ind), ]
  }
  
  # mkuhn, 2015-12-28:
  # when using screening and SV-injection for training data
  # we use fixed holdout for evaluation based on type-column below
  if (! isTRUE(useSimSVInfo) ){
    stopifnot( "type" %in% names(dat.task), setequal(c("train", "test"), levels(dat.task$type)) )
    trainInd <- which(dat.task$type == 'train')
    testInd <- which(dat.task$type == 'test')
    
    stopifnot( length(trainInd) + length(testInd) == NROW(dat.task) )
  }#fi
  
  
  # drop column type for task data (if available)
  dat.task %<>% dplyr::select_(~-one_of("type")) #dat.task2[, names(dat.task2) != 'type']
  
# mkuhn, 2015-11-30: I do not want to handle parallel stuff per function but once per package! cf setParallel()
#   # mkuhn,2015-04-23: turn on parallel computing (if requested)
#   if (is.numeric(nbrCores) && nbrCores >= 2L){
#     #library(parallelMap)
#     on.exit(parallelMap::parallelStop())
#     
# #     # register the mlr-levels
# #     parallelMap::parallelRegisterLevels(package="mlr", levels=c("benchmark", "resample", "selectFeatures", "tuneParams"))
#     #parallelMap::parallelGetRegisteredLevels()
#     if (parallelMode == 'socket') parallelMap::parallelStartSocket(cpus = round(nbrCores), level=parallelLevel) else ##, level = "mlr.benchmark")
#       parallelMap::parallelStartMulticore(cpus = round(nbrCores), level=parallelLevel) ##, level = "mlr.benchmark")
#     
#     #parallelMap::parallelLibrary("ParamHelpers", "FSelector", "e1071") #, "mlr")
#   }#fi nbrCores
#   
  
  # define the Task
  somSV_Task <- mlr::makeClassifTask(id = "somaticSV", data = dat.task, target = "status") %>%
    # feature selection: criterion-based on Task level
    mlr::removeConstantFeatures(perc = 0.001, tol=1e-5, show.info=verbose)
  
  
  
  # Learners ----
  
  # define learners
  ocsvm2_Learner <- mlr::makeLearner("classif.ocsvm2", kernel="radial", config=list(on.learner.error="warn"))
  #svm2_Learner <- mlr::makeLearner("classif.svm", type="nu-classification", kernel="radial", config=list(on.learner.error="warn"))
  svm2_Learner <- mlr::makeLearner("classif.svm", type="C-classification", kernel="radial", config=list(on.learner.error="warn"), predict.type = "prob")
  #ZZZ mkuhn, 2015-12-30: more learners? kernel="poly"
  
  
  # feature selection via FILTER (criterion based) fused to learner
  
  # fuse learner with a feature filter method, based on information gain (entropy)
  ocsvm2_Learner.f <- mlr::makeFilterWrapper(ocsvm2_Learner, fw.method = "rf.importance") 
  svm2_Learner.f <- mlr::makeFilterWrapper(svm2_Learner, fw.method = "rf.importance") %>% 
    mlr::setPredictType(predict.type = "prob")
  #"kruskal.test") #  "information.gain") #, fw.abs=5L)
  
  
#   ## TEST
#   tuning.ps <-  ParamHelpers::makeParamSet(
#     ParamHelpers::makeNumericParam("fw.perc", lower=0.1, upper = 0.9),
#     ParamHelpers::makeDiscreteParam("nu", values = c(0.05, 0.25)), #c(0.01, 0.05, 0.1, 0.25, 0.5)),
#     ParamHelpers::makeDiscreteParam("gamma", values = c(0.0625, 0.125)) #c(0.0625, 0.125, 0.5, 1, 2)) #-3:1, trafo = function(x) 2^x)
#   )
#   
#   resampDesc.inner <-  mlr::makeResampleDesc("CV", iters = 5L, stratify = TRUE)
#   #resampDesc.inner <-  mlr::makeResampleDesc("Subsample", split=4/5, iters = 10L, stratify = TRUE)
#   #r <- resample(learner = svm2_Learner.f, task = somSV_Task, resampling = resampDesc.inner, show.info = TRUE, models = TRUE)
#   
#   # tuning grid search
#   tuneCtrl <-  mlr::makeTuneControlGrid(resolution = 5L)
#   
#   
# #   iris.task <- mlr::makeClassifTask(id="iris", data=iris, target="Species")
# #   lrn = mlr::makeFilterWrapper(learner = "classif.fnn", fw.method = "information.gain", fw.abs = 2)
# #   rdesc = mlr::makeResampleDesc("CV", iters = 10)
# #   r = mlr::resample(learner = lrn, task = iris.task, resampling = rdesc, show.info = FALSE, models = TRUE)
# #   r$aggr
#   
#   ## TEST END
  
  
  
  # Tuning ----

  tuning.ps.ocsvm2 <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.33, upper = 0.9)
    , ParamHelpers::makeNumericParam("nu", lower=0.05, upper = 0.5)
    #, ParamHelpers::makeNumericParam("cost", lower=0.5, upper = 15)
    #, ParamHelpers::makeDiscreteParam("kernel", values = c("radial", "polynomial"))
    , ParamHelpers::makeNumericParam("gamma", lower=2**-7, upper = 2**4) 
    #, ParamHelpers::makeNumericParam("gamma", lower=2**-8, upper = 2**8, requires = quote(kernel != 'linear')) #trafo = function(x) 2^x, 
    #, ParamHelpers::makeIntegerParam("degree", lower = 2L, upper = 5L, requires = expression(kernel=="polynomial"))
    #    , ParamHelpers::makeDiscreteParam("nu", values = c(0.05, 0.25)) #c(0.01, 0.05, 0.1, 0.25, 0.5)),
    #    , ParamHelpers::makeDiscreteParam("gamma", values = c(0.0625, 0.125)) #c(0.0625, 0.125, 0.5, 1, 2)) #-3:1, trafo = function(x) 2^x)
  )
  
  tuning.ps.svm2 <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.33, upper = 0.9)
    #, ParamHelpers::makeNumericParam("nu", lower=0.05, upper = 0.5)
    , ParamHelpers::makeNumericParam("cost", lower=0.5, upper = 15)
    #, ParamHelpers::makeDiscreteParam("kernel", values = c("radial", "polynomial"))
    , ParamHelpers::makeNumericParam("gamma", lower=2**-7, upper = 2**4) 
    #, ParamHelpers::makeNumericParam("gamma", lower=2**-8, upper = 2**8, requires = quote(kernel != 'linear')) #trafo = function(x) 2^x, 
    #, ParamHelpers::makeIntegerParam("degree", lower = 2L, upper = 5L, requires = expression(kernel=="polynomial"))
#    , ParamHelpers::makeDiscreteParam("nu", values = c(0.05, 0.25)) #c(0.01, 0.05, 0.1, 0.25, 0.5)),
#    , ParamHelpers::makeDiscreteParam("gamma", values = c(0.0625, 0.125)) #c(0.0625, 0.125, 0.5, 1, 2)) #-3:1, trafo = function(x) 2^x)
  )
  
  
  # tune control
  
  tuneCtrl.ocsvm2 <- switch(tuneMode,
                            random = mlr::makeTuneControlRandom(maxit = 50L),  ##random search
                            gensa = mlr::makeTuneControlGenSA(maxit = 50L),
                            grid = , ##grid search as default
                            mlr::makeTuneControlGrid(resolution = 4L) #c(fw.perc=2L, nu=4L, gamma=3L))
  )
   
  tuneCtrl.svm2 <- switch(tuneMode,
                            random = mlr::makeTuneControlRandom(maxit = 50L),
                            gensa = mlr::makeTuneControlGenSA(maxit = 50L),
                            grid = ,
                            mlr::makeTuneControlGrid(resolution = 2L, tune.threshold = TRUE) #c(fw.perc=2L, cost=4L, gamma=3L))
  )
  
  ## irace
  #tuneCtrl = mlr::makeTuneControlIrace(maxExperiments = 200L)
  
  # tuning grid search
  #tuneCtrl <-  mlr::makeTuneControlGrid(resolution = 5L)
  
  
  # re-sampling for tuning
  #resampDesc.inner <-  mlr::makeResampleDesc("CV", iters = 5L, stratify = TRUE)
  resampDesc.inner <-  mlr::makeResampleDesc("Subsample", split=.9, iters = k, stratify = TRUE)
  
  # # tune the hyperparameters (by hand)
  # res <- tuneParams(ocsvm2_Learner, task = somSV_Task, resampling = resampDesc.inner, par.set = tuning.ps,
  #                   control = tuneCtrl, measures = list(mlr::acc, setAggregation(mlr::acc, test.sd)), show.info = FALSE)
  # # tune by hand
  # svm2_tuneObj <- tuneParams(learner = svm2_Learner.f, task = somSV_Task, resampling = resampDesc.inner, par.set = tuning.ps, measures = list(my.mcc, mlr::acc, mlr::tpr), control = tuneCtrl, show.info=TRUE)
  
  # tuned learners
  ocsvm2_Learner.f.t <- mlr::makeTuneWrapper(ocsvm2_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                             par.set = tuning.ps.ocsvm2, control = tuneCtrl.ocsvm2, show.info=verbose)
  svm2_Learner.f.t <- mlr::makeTuneWrapper(svm2_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                           par.set = tuning.ps.svm2, control = tuneCtrl.svm2, show.info=verbose)
  
  
  # mkuhn, 2015-12-29: quick evaluation on training data
  if (FALSE){
    trained <- mlr::train(svm2_Learner.f.t, task=somSV_Task, subset = trainInd)
    predObj <- predict(trained, task = somSV_Task, subset = trainInd)
    predObj <- predict(trained, task = somSV_Task, subset = testInd)
    #plotLearnerPrediction(svm2_Learner.f.t, task=somSV_Task) #takes time again!
    table( mlr::getPredictionTruth(predObj), mlr::getPredictionResponse(predObj))
    mlr::getConfMatrix(predObj)
    mlr::performance(predObj, measures = list(my.mcc, mlr::mcc, mlr::acc, mlr::tpr, mlr::tnr))
  }
  
  
  
  # Evaluation (benchmarking) -----
  
  # choose evaluation-resampling (outer loop)
  resampDesc.outer <- if (isTRUE(useSimSVInfo)){
    logging::loginfo("Using subsampling for evaluation (outer loop)..")
    #resampDesc.outer <-  mlr::makeResampleDesc("CV", iters=2L, stratify = TRUE)
    mlr::makeResampleDesc("Subsample", split=.9, iters=k, stratify = TRUE)
  } else {
    logging::loginfo("Using screen-positive regions as fixed holdout test data for evaluation..")
    stopifnot( length(trainInd) + length(testInd) == NROW(dat.task) )
    mlr::makeFixedHoldoutInstance(train.inds=trainInd, test.inds=testInd, size=NROW(dat.task))
  }
  
  
#   # mkuhn, 20150428: testing for parallel
#   r <- mlr::resample(svm2_Learner.f.t, task = somSV_Task, resampling = resampDesc.outer, measures = list(my.mcc, mlr::mcc, mlr::acc), show.info=TRUE, models=TRUE, extract=mlr::getTuneResult)
#   return(r)
  
  #list(ocsvm2_Learner.f.t, svm2_Learner.f.t)
  bmark <- mlr::benchmark(learners = list(svm2_Learner.f.t), show.info = verbose,
                          tasks = somSV_Task, resamplings = resampDesc.outer, measures = list(my.mcc, mlr::mcc, mlr::acc, mlr::tpr))                   
  
  
  # write out benchmark
  if (isTRUE(saveEvaluation)){
    benchmarkFile <- file.path(getVirtualPatientPath(basedir = baseDir, what = "performance", .sample.prop = sample.prop, .ngs.prop = ngs.prop, create=TRUE),
                               paste0(PERF_PREFIX, "_learner", if (! is.null(margin.bp)) paste0("_margin:", margin.bp), ".rds")) #"_ssv:", round(prop.sSV*100,0L),".rds"))
    saveRDS(bmark, file=benchmarkFile)
    loginfo("Saved benchmark object to %s.", benchmarkFile)
  }
  
  
  mlr::getBMRPerformances(bmark, as.df=TRUE)
}


#' Prepares the feature data frame for task.
#' @param strategy how to divide into training and test data
prepareFeatureData <- function(data, strategy=c("all", "pat1"), doSubSampling=FALSE, subSampling.prop.sSV=NULL, verbose=FALSE){
  strategy <- match.arg(strategy)
  
  stopifnot( is.data.frame(data), NROW(data) > 0L )
  dat.task <- data %>% 
    dplyr::filter_(~complete.cases(.))
  
  
  # mkuhn, 2015-04-23: re-name factor levels as the learners expect IN and OUT
  if ( identical(levels(dat.task$status), c("germline", "DSV")) ){
    dat.task$status <- factor(dat.task$status, levels = c("germline", "DSV"), labels = STATUS_LEVELS)
    logging::logwarn("Recode status vector in: %s", paste(STATUS_LEVELS, collapse=" & "))
  }
  
  stopifnot( identical(levels(dat.task$status), STATUS_LEVELS) )
  
  # down-/up- weighting cases?
  if ( isTRUE(doSubSampling) ){
    stopifnot( ! is.null(subSampling.prop.sSV) )
    
    # subsampling: training data with fixed proportion of somatic SVs
    logging::loginfo("Subsampling from available data to get somatic SV proportion prop.sSV=%.3f", subSampling.prop.sSV)
    #ssInd <- subsampleInd(dat.status = dat.task$status == STATUS_LEVELS[1], prop.outl = prop.sSV)
    ssInd <- subsampleInd(dat.status = dat.task$status, prop.outl = subSampling.prop.sSV)
    stopifnot( length(ssInd$na.ind) == 0L ) # we should not have NAs any more (mkuhn, 2015-12-29)
    dat.task <- dat.task[c(ssInd$train.ind), ]
  }
  
  
  stopifnot( "type" %in% names(dat.task), setequal(c("train", "test"), levels(dat.task$type)) )
  
  
  
  if (strategy == 'all'){
    trainInd <- which(dat.task$type == 'train')
    testInd <- which(dat.task$type == 'test')
    stopifnot( length(trainInd) + length(testInd) == NROW(dat.task), length(trainInd) > 0L )
  } else {
    # pat1
    # use training cases of pat1 for training
    # use rest for testing (test of pat1 and stuff of pat2..)
    
    # fix outcome for test cases
    dat.train <- dplyr::filter_(dat.task, ~type == 'train')
    dat.test0 <- dplyr::filter_(dat.task, ~type == 'test')
    stopifnot( NROW(dat.train) + NROW(dat.test0) == NROW(dat.task))
    
    patData <- readRDS(file.path(getVirtualPatientPath(what = "patient"), "patData.rds"))
    
    dat.test <- addStatusToTestFromSim(featDat.test = dat.test0, patData = patData, myCHROM = CHROM, 
                           evalMargin = 300L, featNames = names(dat.train))
    
    dat.task <- rbind(dat.train, dat.test)
    
    stopifnot( strategy == 'pat1', sum(dat.task$patId == 1L) > 0L )
    logging::loginfo("Taking %d training entries of pat1 for training. Rest for testing..", sum(dat.task$patId == 1L & dat.task$type == 'train'))
    
    trainInd <- which(dat.task$type == 'train' & dat.task$patId == 1L)
    testInd <- setdiff( 1:NROW(dat.task), trainInd )
    
  }
  
  
  # drop column "type" for task data (if available)
  dat.task %<>% 
    dplyr::select_(~one_of("status"), ~ends_with(".r"), ~ends_with(".d"), ~ends_with(".t")) 
  
  
  # define the Task
  taskObj <- mlr::makeClassifTask(id = "somaticSV", data = dat.task, target = "status", positive = "OUT") %>%
    # feature selection: criterion-based on Task level
    # drop near constant features (with at max 1/1000 observations different from constant)
    mlr::removeConstantFeatures(perc = 1e-3, tol=1e-4, show.info=verbose)
  
  
  list(task=taskObj, trainInd=trainInd, testInd=testInd)
}



makeLearnerList <- function(short=FALSE){
  # define learners
  ocsvm2_Learner <- mlr::makeLearner("classif.ocsvm2", config=list(on.learner.error="warn"))
  #svm2_Learner <- mlr::makeLearner("classif.svm", type="nu-classification", kernel="radial", config=list(on.learner.error="warn"))
  svm2_Learner <- mlr::makeLearner("classif.svm", type="nu-classification", config=list(on.learner.error="warn"), predict.type = "prob")
  kknn_Learner <- mlr::makeLearner("classif.kknn", scale=TRUE, config=list(on.learner.error="warn"))
  pamr_Learner <- mlr::makeLearner("classif.pamr", config=list(on.learner.error="warn"))
  glmnet_Learner <- mlr::makeLearner("classif.glmnet", standardize=TRUE, config=list(on.learner.error="warn"))
  nnet_Learner <- mlr::makeLearner("classif.nnet", config=list(on.learner.error="warn"))
  logreg_Learner <- mlr::makeLearner("classif.logreg", config=list(on.learner.error="warn"))
  rf_Learner <- mlr::makeLearner("classif.randomForest", config=list(on.learner.error="warn"))
  
  
  # feature selection via FILTER (criterion based) fused to learner
  
  # fuse learner with a feature filter method, based on information gain (entropy)
  myFilterCriterion <- "cforest.importance"
  ocsvm2_Learner.f <- mlr::makeFilterWrapper(ocsvm2_Learner, fw.method = myFilterCriterion) 
  svm2_Learner.f <- mlr::makeFilterWrapper(svm2_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  kknn_Learner.f <- mlr::makeFilterWrapper(kknn_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  pamr_Learner.f <- mlr::makeFilterWrapper(pamr_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  glmnet_Learner.f <- mlr::makeFilterWrapper(glmnet_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  nnet_Learner.f <- mlr::makeFilterWrapper(nnet_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  logreg_Learner.f <- mlr::makeFilterWrapper(logreg_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  # rf filter prior to rf?
  rf_Learner.f <- mlr::makeFilterWrapper(rf_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  #"kruskal.test") #  "information.gain") #, fw.abs=5L)
  
  # list of filtered learners
  list(ocsvm2=ocsvm2_Learner.f, svm2=svm2_Learner.f, knn=kknn_Learner.f, #pamr_Learner.f,
       elasticNet=glmnet_Learner.f, #nnet_Learner.f,
       logReg=logreg_Learner.f, randomForest=rf_Learner.f) 
}



#' Evaluate some learners on a feature set, based on simulation results.
#' 
#' This is a 2nd attempt: started off as a copy from evalLearners(). but has more learners.
#' 
#' The feature set from simulation run is specified by the sample- and NGS-properties. It writes out the benchmark object in the corresponding performance/ -directory (if requested).
#' One needs also to give the margin that was used left and right of the screened or simulated SV-regions for features extraction from the mapping.
#' \describe{
#'   \item{with \code{useSimSVInfo=TRUE}}{the simulated SV-information is directly used, i.e. no SV-screening step is done but use features extracted directly at the simulated SV-locus.}
#'   \item{with \code{useSimSVInfo=FALSE}}{the SV-screening step has been performed that resulted in screening negative regions that give rise to positive (via SV-injection) and negative training cases.
#'   The screening positive regions form the test set where the learners naturally need to be evaluated.}
#' }
#' Patient IDs do not play a role for simulation runs, except you can restrict number of features to come from certain simulation patients.
#' On real data the algorithm is run per patient ("personalized").
#' 
#' The \code{sources}-list defines the feature-data source:
#' sample.prop sample properties underlying the features
#' ngs.prop NGS properties underlying the features
#' patIds numeric vector of patient IDs to restrict features to only given patients. Defaults to \code{NULL}, i.e. all available features from any patient.
#' margin.bp margin size (in bp) that was used for feature extraction (either flanking the exact simulated SV or flanking the screened regions). The margin is part of the feature filename and is used to identify the features to use.
#' useSimSVInfo flag. For simulation case: Do we want to use the knowledge of SV-presence and SV-coordinates from simulation for the evaluation? Defaults to \code{FALSE}, i.e. use read coverage and clipped base peaks screening to learn regions where features are extracted.
#' 
#' @param sources list that contains information for which data (training and test) to do the evaluation. E.g. \code{isSimulation=} governs if it is simulated data or biotec data
#' @param k number of sub-sampling iterations for tuning (inner loop). Defaults to 10.
#' @param tuneMode character. name of supported tuning method
#' @param tune.niter number of iterations for tuning (no effect on grid, as this is set independent of tune.niter)
#' @param prop.sSV proportion of somatic SVs in the genomic regions in training data
#' @param saveEvaluation flag if the benchmark is to be written out to disk in the performance/-directory.
#' @param shortEval flag if you want only quick evaluation
#' @param suffix character filename suffix, optional
#' @param verbose logical flag: show info during tuning/evaluation
#' @return dataframe with benchmark performance results.
#' @export
evalLearners2 <- function(baseDir=BASEDIR, sources=list(isSimulation=TRUE, sampleProp=sample.prop(), ngsProp=ngs.prop(), patIds=NULL, useSimSVInfo=FALSE, margin.bp=25L), 
                         #parallelMode=c("socket", "multicore"), nbrCores=2L, parallelLevel=c("mlr.tuneParams", "mlr.resample", "mlr.benchmark", "mlr.selectFeatures"),
                         tuneMode=c("grid", "random", "gensa", "irace"), tune.niter=100L, shortEval=FALSE,
                         k=10L, prop.sSV=0.05, saveEvaluation=TRUE, suffix=NULL, verbose=FALSE){
  
  # parallelMode <- match.arg(parallelMode)
  # parallelLevel <- match.arg(parallelLevel)
  tuneMode <- match.arg(tuneMode)
  
  # mkuhn, 2015-12-30 decide here what is best?! Do I want to sub-sample the majority class? or use mlr-Wrapper? or nothing?
  # default: we do not do case-upweighting or down-weighting by hand.
  doSubSampling <- FALSE
  
  
  # for quick start only to set arguments
  if (FALSE){
    baseDir <- BASEDIR
    sources <- list(isSimulation=TRUE, sampleProp=sample.prop(cov=60L, tumorLoad=90L), ngsProp=ngs.prop(), patIds=NULL, useSimSVInfo=FALSE, margin.bp=25L)
    tuneMode <- "random" #"irace"
    tune.niter <- 100L
    #nbrCores=3L; parallelMode <- "socket";  parallelLevel <- "mlr.tuneParams"
    prop.sSV=0.05; 
    shortEval <- TRUE
    k <- 5L; saveEvaluation=FALSE; suffix <- NULL; verbose=TRUE
  }#fi
  
  stopifnot( is.list(sources), all(c("isSimulation", "useSimSVInfo", "margin.bp") %in% names(sources)), is.numeric(sources[["margin.bp"]]) )
  
  
  
  
  # Read in data -----
  
  patStr <- ""
  
  # feature RDS-file
  featFile <- if ( isTRUE(sources[["isSimulation"]]) ){
    
    featFileBase <- sprintf("features_sim_%s_margin%d.rds", if (isTRUE(sources[["useSimSVInfo"]])) "noScreen" else "screen", sources[["margin.bp"]])
    #paste0(if (isTRUE(sources[["useSimSVInfo"]])) "features_sim_margin" else "features_screen_margin", sources[["margin.bp"]], ".rds")
    file.path(getVirtualPatientPath(basedir = baseDir, what = "feature", .sample.prop = sources[["sampleProp"]], .ngs.prop = sources[["ngsProp"]]), featFileBase )
  } else {
    
    # biotec data: need a patient ID
    stopifnot( is.numeric(sources[["patIds"]]), length(sources[["patIds"]]) == 1L )
    patStr <- paste0("_P", sources[["patIds"]][1])
    file.path(baseDir, "biotec", "feature", paste0("features_real_screen", patStr, "_margin", sources[["margin.bp"]], ".rds"))
  }
    
  
  if (! isTRUE(file.exists(featFile)) ){
    logging::logerror("Feature data file %s not found!", featFile)
    return(invisible(NULL))
  }
  
  data <- readRDS(featFile)
  stopifnot( exists("data"), is.data.frame(data) )  
  
  
    
    
  
  # mkuhn, 20150504: optional patient filter on feature dataframe (useful only for sim-cases as they have feature file with multiple simulated patients.
  # biotec data is anyway per individual patient
  if ( isTRUE(sources[["isSimulation"]]) && ! is.null(sources[["patIds"]]) && is.numeric(sources[["patIds"]]) && length(sources[["patIds"]]) >= 1L ){
    data <- dplyr::filter_(data, ~patId %in% sources[["patIds"]])
    logging::loginfo("Restricted to %d features to given patients.", NROW(data))
  }
  
  # check if patient info is available within feature data
  if ( ! (is.data.frame(data) && NROW(data) >= 1L && "status" %in% names(data) && is.factor(data$status) && nlevels(data$status) == 2L) ){
    logging::logwarn("Something wrong with feature data in %s. Status missing?", featFile)
    return(invisible(NULL))
  }
    
  
  
  
  
  
  # Task ----
  
  # prepare data for learner
  taskList <- prepareFeatureData(data, doSubSampling = doSubSampling, subSampling.prop.sSV = prop.sSV)
  somSV_Task <- taskList[["task"]]
  trainInd <- taskList[["trainInd"]]
  testInd <- taskList[["testInd"]]  
  
  
  # mkuhn, 2015-12-28:
  # with real data OR with simulation with screening step and SV-injection for training data
  # we use screen-positive data as fixed holdout for evaluation based on type-column below
  useScreenPosTestdataForEvaluation <- ! ( isTRUE(sources[["isSimulation"]]) && isTRUE(sources[["useSimSVInfo"]]) )
  
  
  if ( isTRUE(useScreenPosTestdataForEvaluation) ){
    logging::loginfo("Using screen-positive test data for evaluation")

    if ( length(testInd) == 0L ){
      logging::logwarn("Test data is missing. But I needed it as fixed holdout data..")
      return(invisible(NULL))
    }
    
  }#fi
  
  
  
  
  # Learners ----
  
  # define learners
  ocsvm2_Learner <- mlr::makeLearner("classif.ocsvm2", config=list(on.learner.error="warn"))
  #svm2_Learner <- mlr::makeLearner("classif.svm", type="nu-classification", kernel="radial", config=list(on.learner.error="warn"))
  svm2_Learner <- mlr::makeLearner("classif.svm", type="nu-classification", config=list(on.learner.error="warn"), predict.type = "prob")
  kknn_Learner <- mlr::makeLearner("classif.kknn", scale=TRUE, config=list(on.learner.error="warn"))
  pamr_Learner <- mlr::makeLearner("classif.pamr", config=list(on.learner.error="warn"))
  glmnet_Learner <- mlr::makeLearner("classif.glmnet", standardize=TRUE, config=list(on.learner.error="warn"))
  nnet_Learner <- mlr::makeLearner("classif.nnet", config=list(on.learner.error="warn"))
  logreg_Learner <- mlr::makeLearner("classif.logreg", config=list(on.learner.error="warn"))
  rf_Learner <- mlr::makeLearner("classif.randomForest", config=list(on.learner.error="warn"))
  
  
  # feature selection via FILTER (criterion based) fused to learner
  
  # fuse learner with a feature filter method, based on information gain (entropy)
  myFilterCriterion <- "cforest.importance"
  ocsvm2_Learner.f <- mlr::makeFilterWrapper(ocsvm2_Learner, fw.method = myFilterCriterion) 
  svm2_Learner.f <- mlr::makeFilterWrapper(svm2_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  kknn_Learner.f <- mlr::makeFilterWrapper(kknn_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  pamr_Learner.f <- mlr::makeFilterWrapper(pamr_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  glmnet_Learner.f <- mlr::makeFilterWrapper(glmnet_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  nnet_Learner.f <- mlr::makeFilterWrapper(nnet_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  logreg_Learner.f <- mlr::makeFilterWrapper(logreg_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  # rf filter prior to rf?
  rf_Learner.f <- mlr::makeFilterWrapper(rf_Learner, fw.method = myFilterCriterion) %>% 
    mlr::setPredictType(predict.type = "prob")
  #"kruskal.test") #  "information.gain") #, fw.abs=5L)
  
  
  
  
  
  # Tuning ----
  
  tuning.ps.ocsvm2 <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
    , ParamHelpers::makeNumericParam("nu", lower=0.05, upper = 0.5)
    #, ParamHelpers::makeNumericParam("cost", lower=0.5, upper = 15)
    , ParamHelpers::makeDiscreteParam("kernel", values = c("linear", "radial", "polynomial"))
    #, ParamHelpers::makeNumericParam("gamma", lower=2**-7, upper = 2**4) 
    , ParamHelpers::makeNumericParam("gamma", lower=-8, upper = 8, requires = quote(kernel != 'linear'), trafo = function(x) 2^x)
    , ParamHelpers::makeIntegerParam("degree", lower = 2L, upper = 3L, requires = quote(kernel=="polynomial"))
    #    , ParamHelpers::makeDiscreteParam("nu", values = c(0.05, 0.25)) #c(0.01, 0.05, 0.1, 0.25, 0.5)),
    #    , ParamHelpers::makeDiscreteParam("gamma", values = c(0.0625, 0.125)) #c(0.0625, 0.125, 0.5, 1, 2)) #-3:1, trafo = function(x) 2^x)
  )
  
  tuning.ps.svm2 <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
    , ParamHelpers::makeNumericParam("nu", lower=0.05, upper = 0.75)
    #, ParamHelpers::makeNumericParam("cost", lower=0.5, upper = 15)
    , ParamHelpers::makeDiscreteParam("kernel", values = c("linear", "radial", "polynomial"))
    , ParamHelpers::makeNumericParam("gamma", lower=-8, upper = 8, requires = quote(kernel != 'linear'), trafo = function(x) 2^x)
    , ParamHelpers::makeIntegerParam("degree", lower = 2L, upper = 3L, requires = quote(kernel=="polynomial"))
  )
  
  
  tuning.ps.kknn <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
    , ParamHelpers::makeIntegerParam("k", lower=2, upper = 15)
    , ParamHelpers::makeDiscreteParam("kernel", values = c("rectangular", "epanechnikov")) #"optimal", # optimal is not feasible for parameter 'kernel'
  )
  
  tuning.ps.pamr <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
  )
  
  tuning.ps.glmnet <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
    , ParamHelpers::makeNumericParam("alpha", lower=0.1, upper = 1)
    , ParamHelpers::makeNumericParam("lambda.min.ratio", lower=-10, upper = -3, trafo = function(x) 2^x)
  )
  
  tuning.ps.nnet <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
    , ParamHelpers::makeIntegerParam("size", lower=0, upper = 9)
    , ParamHelpers::makeNumericParam("decay", lower=-2, upper = 2)
  )
  tuning.ps.logreg <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
  )
  tuning.ps.rf <-  ParamHelpers::makeParamSet(
    ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
    , ParamHelpers::makeIntegerParam("ntree", lower = 5, upper = 10, trafo = function(x) 2^x)
    , ParamHelpers::makeIntegerParam("mtry", lower = 5, upper = 15)
  )
  
  
  myResolution <- 3L
  if (isTRUE(shortEval)){
    k <- min(k,3L)
    myResolution <- min(myResolution, 2L)
  }
  
  # tune control: tune.threshold available if learner supports probability prediction
  
  tuneCtrl.notuneTh <- switch(tuneMode,
                            random = mlr::makeTuneControlRandom(maxit = tune.niter, tune.threshold = FALSE),  ##random search
                            gensa = mlr::makeTuneControlGenSA(maxit = tune.niter, tune.threshold = FALSE),
                            irace = mlr::makeTuneControlIrace(maxExperiments = tune.niter, tune.threshold = FALSE),
                            grid = , ##grid search as default
                            mlr::makeTuneControlGrid(resolution = myResolution, tune.threshold = FALSE) #c(fw.perc=2L, nu=4L, gamma=3L))
  )
  
  tuneCtrl.tuneTh <- switch(tuneMode,
                          random = mlr::makeTuneControlRandom(maxit = tune.niter, tune.threshold = TRUE),
                          gensa = mlr::makeTuneControlGenSA(maxit = tune.niter, tune.threshold = TRUE),
                          irace = mlr::makeTuneControlIrace(maxExperiments = tune.niter, tune.threshold = TRUE),
                          grid = ,
                          mlr::makeTuneControlGrid(resolution = myResolution, tune.threshold = TRUE) #c(fw.perc=2L, cost=4L, gamma=3L))
  )
  
  
  
  # re-sampling for tuning
  #resampDesc.inner <-  mlr::makeResampleDesc("CV", iters = 5L, stratify = TRUE)
  resampDesc.inner <-  mlr::makeResampleDesc("Subsample", split=.75, iters = k, stratify = TRUE)
  
  # # tune the hyperparameters (by hand)
  # res <- tuneParams(ocsvm2_Learner, task = somSV_Task, resampling = resampDesc.inner, par.set = tuning.ps,
  #                   control = tuneCtrl, measures = list(mlr::acc, setAggregation(mlr::acc, test.sd)), show.info = FALSE)
  # # tune by hand
  # svm2_tuneObj <- tuneParams(learner = svm2_Learner.f, task = somSV_Task, resampling = resampDesc.inner, par.set = tuning.ps, measures = list(my.mcc, mlr::acc, mlr::tpr), control = tuneCtrl, show.info=TRUE)
  
  # tuned learners
  ocsvm2_Learner.f.t <- mlr::makeTuneWrapper(ocsvm2_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                             par.set = tuning.ps.ocsvm2, control = tuneCtrl.notuneTh, show.info=verbose)
  svm2_Learner.f.t <- mlr::makeTuneWrapper(svm2_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                           par.set = tuning.ps.svm2, control = tuneCtrl.tuneTh, show.info=verbose)
  kknn_Learner.f.t <- mlr::makeTuneWrapper(kknn_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                           par.set = tuning.ps.kknn, control = tuneCtrl.tuneTh, show.info=verbose)
  pamr_Learner.f.t <- mlr::makeTuneWrapper(pamr_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                           par.set = tuning.ps.pamr, control = tuneCtrl.tuneTh, show.info=verbose)
  glmnet_Learner.f.t <- mlr::makeTuneWrapper(glmnet_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                           par.set = tuning.ps.glmnet, control = tuneCtrl.tuneTh, show.info=verbose)
  nnet_Learner.f.t <- mlr::makeTuneWrapper(nnet_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                             par.set = tuning.ps.nnet, control = tuneCtrl.tuneTh, show.info=verbose)
  logreg_Learner.f.t <- mlr::makeTuneWrapper(logreg_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                           par.set = tuning.ps.logreg, control = tuneCtrl.tuneTh, show.info=verbose)
  rf_Learner.f.t <- mlr::makeTuneWrapper(rf_Learner.f, resampDesc.inner, measures = list(my.mcc),
                                             par.set = tuning.ps.rf, control = tuneCtrl.tuneTh, show.info=verbose)
  
  
  learnerList <- list(ocsvm2_Learner.f.t, svm2_Learner.f.t, kknn_Learner.f.t,
             #pamr_Learner.f.t, nnet_Learner.f.t, # these failed!
             glmnet_Learner.f.t,
             logreg_Learner.f.t, rf_Learner.f.t)
  
  
  
  if (isTRUE(shortEval)){
    
    learnerList <- list(#ocsvm2_Learner.f.t, svm2_Learner.f.t,
      kknn_Learner.f.t,
      #pamr_Learner.f.t, nnet_Learner.f.t, # these failed!
      #glmnet_Learner.f.t,
      logreg_Learner.f.t, rf_Learner.f.t)
  }
  
  
  # mkuhn, 2015-12-29: quick evaluation on training or test data
  if ( FALSE && exists("trainInd") && exists("testInd") ){
    # quick check if filter-wrapper works properly..
    kknnL <- mlr::makeLearner("classif.kknn")
    kknnL.f <- mlr::makeFilterWrapper(kknnL, fw.method = myFilterCriterion) %>% 
      mlr::setPredictType(predict.type = "prob") %>% mlr::setHyperPars(fw.perc = .8)
    
    trained <- mlr::train(kknnL, task=somSV_Task, subset = trainInd)
    trained.f <- mlr::train(kknnL.f, task=somSV_Task, subset = trainInd)
    
    trained <- mlr::train(learnerList[[3L]], task=somSV_Task, subset = trainInd)
    #predObj <- predict(trained, task = somSV_Task, subset = trainInd)
    predObj <- predict(trained, task = somSV_Task, subset = testInd)
    #plotLearnerPrediction(svm2_Learner.f.t, task=somSV_Task) #takes time again!
    table( mlr::getPredictionTruth(predObj), mlr::getPredictionResponse(predObj))
    mlr::getConfMatrix(predObj)
    mlr::performance(predObj, measures = list(my.mcc, mlr::mcc, mlr::acc, mlr::tpr, mlr::tnr))
  }
  
  
  
  # Evaluation (benchmarking) -----
  
  # choose evaluation-resampling (outer loop)
  resampDesc.outer <- if ( isTRUE(useScreenPosTestdataForEvaluation) ){
    logging::loginfo("Using screen-positive regions as fixed holdout test data for evaluation..")
    stopifnot( length(trainInd) + length(testInd) == NROW(mlr::getTaskData(somSV_Task)) )
    mlr::makeFixedHoldoutInstance(train.inds=trainInd, test.inds=testInd, size=NROW(mlr::getTaskData(somSV_Task)))
  } else {
    logging::loginfo("Using subsampling for evaluation (outer loop)..")
    #resampDesc.outer <-  mlr::makeResampleDesc("CV", iters=2L, stratify = TRUE)
    mlr::makeResampleDesc("Subsample", split=.8, iters=k, stratify = TRUE)
  }
  
  
  #   # mkuhn, 20150428: testing for parallel
  #   r <- mlr::resample(svm2_Learner.f.t, task = somSV_Task, resampling = resampDesc.outer, measures = list(my.mcc, mlr::mcc, mlr::acc), show.info=TRUE, models=TRUE, extract=mlr::getTuneResult)
  #   return(r)
  
  # # TESTING: bmark as list of individual benchmark objects
  # bmark <- foreach::foreach(l=learnerList, .combine = list, .multicombine = TRUE, .errorhandling = "pass") %do% {
  #   cat("\nLearner ", l$id)
  # 
  #   mlr::benchmark(learners = l, show.info = verbose, 
  #                         tasks = somSV_Task, resamplings = resampDesc.outer, measures = list(my.mcc, mlr::mcc, mlr::acc, mlr::tpr))                   
  # }
  
  # mkuhn, 2016-01-11: create a single benchmark object
  bmark <- mlr::benchmark(learners = learnerList, show.info = verbose, 
                   tasks = somSV_Task, resamplings = resampDesc.outer, measures = list(my.mcc, mlr::mcc, mlr::acc, mlr::tpr, mlr::tnr))
  
  
  # Write out benchmark ------
  if (isTRUE(saveEvaluation)){
    
    evalFilename <- paste0(PERF_PREFIX, "_learner", patStr, if (! is.null(sources[["margin.bp"]])) paste0("_margin:", sources[["margin.bp"]]),
           "_tune:", tuneMode,
           if (isTRUE(doSubSampling)) paste0("_ssv:", round(prop.sSV*100,0L)),
           if (! is.null(suffix) && is.character(suffix) && nzchar(suffix)) paste0("_", suffix),".rds")
    
    
    benchmarkDir <- if ( isTRUE(sources[["isSimulation"]]) ) {
      file.path(getVirtualPatientPath(basedir = baseDir, what = "performance", .sample.prop = sources[["sampleProp"]], .ngs.prop = sources[["ngsProp"]], create=TRUE))
    } else {
      file.path(baseDir, "biotec", "performance")
    }
    
    
    try(dir.create(benchmarkDir, showWarnings = FALSE, recursive = TRUE))
    saveRDS(bmark, file=file.path(benchmarkDir, evalFilename))
    loginfo("Saved benchmark object to file %s in %s.", evalFilename, benchmarkDir)
  }
  
  
  ## TESTING
  #mlr::getBMRPerformances(bmark, as.df=TRUE)
  bmark #list
}




#' Evaluate tuned learners from benchmark-experiment form simulation on new test data
#' 
#' It uses new positive and negative cases from unseen patients for the evaluation.
#' It needs both the feature dataframe and the performance results (for the tuned parameters)
#' @return dataframe with performance results of tuned learners on new test data
evalTunedLearners <- function(sampleProp, ngsProp=ngs.prop(), featureMargin){
  
  featPath <- getVirtualPatientPath(what = "feature", .sample.prop=sampleProp, .ngs.prop = ngsProp)
  perfPath <- getVirtualPatientPath(what = "performance", .sample.prop=sampleProp, .ngs.prop = ngsProp)
  
  featFile <- file.path(featPath, paste0("features_sim_screen_margin", featureMargin, ".rds"))
  perfFile <- file.path(perfPath, paste0("perf_learner_margin:", featureMargin, "_tune:grid.rds"))
  
  if ( ! file.exists(featFile) ){
    logging::logwarn("Feature file %s for %s with fMargin %d is missing!", featFile, toString(sampleProp), featureMargin)
    return(NULL)
  }
  
  if (!  file.exists(perfFile)){
    logging::logwarn("Performance file %s for %s with fMargin %d is missing!", perfFile, toString(sampleProp), featureMargin)
    return(NULL)
  }
  
  
  featObj <- readRDS(featFile)
  bmrObj <- readRDS(perfFile)
  
  learnerList <- makeLearnerList()
  
  bmLearners <- mlr::getBMRLearners(bmrObj)
  stopifnot( length(bmLearners) == length(learnerList) )
  
  taskList <- prepareFeatureData(data = featObj, strategy = "pat1")
  taskObj <- taskList[["task"]]
  trainInd <- taskList[["trainInd"]]
  testInd <- taskList[["testInd"]]
  
  # learner index lInd
  foreach::foreach(lInd=seq_along(learnerList), .combine = 'rbind') %do% { 
    tuneResult <- getBMRTuneResults(bmrObj)[[1L]][[lInd]][[1L]]
    # set tuned parameters to learner
    bmLearner <- learnerList[[lInd]] %>% setHyperPars(par.vals = tuneResult$x)
    
    #perfDF <- NULL
    bmLearner.trained <- mlr::train(bmLearner, task = taskObj, subset = trainInd)
    bmLearner.pred <- predict(bmLearner.trained, newdata = mlr::getTaskData(taskObj, subset = testInd))
    if ( ! is.null(tuneResult$threshold) && is.numeric(tuneResult$threshold)){
      bmLearner.pred <-  setThreshold(bmLearner.pred, threshold = tuneResult$threshold)
    }
    #mlr::getConfMatrix(bmLearner.pred)
    
    fw.perc <- tuneResult$x[["fw.perc"]]
    if ( is.null(fw.perc)) fw.perc <- -1L
    
    perf.mcc <- performance(bmLearner.pred, measures = my.mcc)
    perf.tn <- performance(bmLearner.pred, measures = mlr::tn)
    perf.tp <- performance(bmLearner.pred, measures = mlr::tp)
    perf.fn <- performance(bmLearner.pred, measures = mlr::fn)
    perf.fp <- performance(bmLearner.pred, measures = mlr::fp)
    perf.tpr <- performance(bmLearner.pred, measures = mlr::tpr)
    perf.tnr <- performance(bmLearner.pred, measures = mlr::tnr)
    perf.acc <- performance(bmLearner.pred, measures = mlr::acc)
    
    
    #perfDF <- rbind(
    dplyr::data_frame(Coverage = sampleProp$cov, TumorLoad=sampleProp$tumorLoad, FeatureMargin=featureMargin,
                      Learner = names(learnerList)[lInd], nbrTrain=length(trainInd), nbrTest=length(testInd),
                      fw.perc=fw.perc, 
                      TN=perf.tn, TP=perf.tp, FN=perf.fn, FP=perf.fp,
                      MCC=perf.mcc, TPR=perf.tpr, TNR=perf.tnr, ACC=perf.acc)
    #  perfDF)
  }
}



#' Helper function to evaluate different tuned learners.
#' @param nIter just a counter that is reported in the dataframe
#' @return dataframe with 
evalTunedLearnersLoop <- function(covFilter=NULL, tlFilter=NULL, fMarginFilter=NULL, nIter=1L){
  
  inputDF <- expand.grid(Coverage=c(20, 40, 60, 80, 100), Tumorload=c(90, 75, 50, 33, 25, 10), featureMargin=c(25, 50))
  
  if (! is.null(covFilter) && is.numeric(covFilter) ){
    inputDF <- dplyr::filter_(inputDF, ~Coverage %in% covFilter)
  }
  
  if (! is.null(tlFilter) && is.numeric(tlFilter) ){
    inputDF <- dplyr::filter_(inputDF, ~Tumorload %in% tlFilter)
  }
  
  if (! is.null(fMarginFilter) && is.numeric(fMarginFilter) ){
    inputDF <- dplyr::filter_(inputDF, ~featureMargin %in% fMarginFilter)
  }
  
  logging::loginfo("Running over %d configurations..", NROW(inputDF))
  sampleProps <- apply(inputDF, 1, function(row) { sample.prop(cov = row[1], tumorLoad = row[2]) })
  fMargins <- inputDF$featureMargin
  
  stopifnot( length(sampleProps) == length(fMargins) )
  
  resDF <- foreach::foreach(sp=sampleProps, fm=fMargins, .combine = 'rbind') %dopar% {
    evalTunedLearners(sampleProp = sp, featureMargin = fm)
  }
  
  resDF$nIter <- nIter
  
  resDF
  
}


#' Helper: Quickly see what files are there..
#' @return dataframe with status information which tuned learners on simulation data are available
checkTunedLearners <- function(){
  covs=c(20, 40, 60, 80, 100); TLs=c(90, 75, 50, 33, 25, 10); featureMargin=c(25, 50)
  inputDF <- expand.grid(Coverage=covs, Tumorload=TLs, featureMargin=featureMargin)
  sampleProps <- apply(inputDF, 1, function(row) { sample.prop(cov = row[1], tumorLoad = row[2])})
  fMargins <- inputDF$featureMargin
  
  stopifnot( length(sampleProps) == length(fMargins) )
  
  ngsProp <- ngs.prop()
  
  checkDF <- foreach::foreach(sp=sampleProps, fm=fMargins, .combine = 'rbind') %do% {
    featPath <- getVirtualPatientPath(what = "feature", .sample.prop=sp, .ngs.prop = ngsProp)
    perfPath <- getVirtualPatientPath(what = "performance", .sample.prop=sp, .ngs.prop = ngsProp)
    
    featFile <- file.path(featPath, paste0("features_sim_screen_margin", featureMargin, ".rds"))
    perfFile <- file.path(perfPath, paste0("perf_learner_margin:", featureMargin, "_tune:grid.rds"))
    
    dplyr::data_frame(
      Coverage = sp$cov, TumorLoad=sp$tumorLoad, FeatureMargin=featureMargin,
      featFile=file.exists(featFile), perfFile=file.exists(perfFile))
  }
  
  checkDF
}



#' Evaluate screening, based on simulation data, for a given sample and NGS-setting.
#' 
#' The feature set from simulation run is specified by the sample- and NGS-properties.
#' It writes out the benchmark object in the corresponding performance/ -directory.
#' In case the feature data stems from a simulation run one needs also to give the margin left and right of the simulated SV
#' that was used prior to extract features from the mapping.
#' 
#' @param sample.prop sample properties underlying the features
#' @param ngs.prop NGS properties underlying the features
#' @param patIds numeric vector of patient IDs to restrict features to only given patients. Defaults to \code{NULL}, i.e. no filter of features
#' @param myTargetBedfile file name of enrichment bed file 
#' @param margin.bp evaluation margin (in bp) used for evaluation of screen: how many bp do I extend the peak position when I decide if it is a hit or not. whole REGION match is also done. margin becomes a part of the file name.
#' @param saveEvaluation flag if the evaluation is to be written out to disk in the performance/-directory.
#' @return dataframe with performance of screen. Or \code{NULL} if not patients selected.
#' @export
evalScreen <- function(baseDir=BASEDIR, patIds=NULL, myTargetBedFile=file.path(baseDir, "refGen", "TruSeq_exome_targeted_regions.hg19.bed"),
                       myChrom=CHROM, sample.prop, ngs.prop, margin.bp=25L, saveEvaluation=FALSE){
  
  stopifnot( file.exists(baseDir), file.exists(myTargetBedFile) )
  stopifnot( is.numeric(margin.bp) && margin.bp >= 0L ) # || margin.bp == 'REGION')
  
  # patient information: which regions were simulated and where are the differential SVs?!
  patDataFile <- file.path(getVirtualPatientPath(baseDir, what = "patient"), "patData.rds")
  
  patData <- readRDS(patDataFile)
  stopifnot( NROW(patData) >= 1L )
  
  patIds <- if (!is.null(patIds) && is.numeric(patIds)) intersect(patIds, unique(patData$patId)) else unique(patData$patId)
  
  if ( length(patIds) == 0L ){
    logging::logwarn("No patients for screen evaluation after patID-filtering. Quitting.")
    return(NULL)
  }#fi
  
  
  logging::loginfo("In %s, found data for %d patients.", patDataFile, length(patIds))
  
  
  
  # mapping directory for given setting of coverage and TL
  mappingDir <- getVirtualPatientPath(baseDir, what = "mapping", .sample.prop=sample.prop, .ngs.prop=ngs.prop)
  
  
  
  # iterate over patients:
  # get for each patient information about his covered regions. Then we assess in those regions clipped base peaks
  evalScreenDat <- foreach(pId=patIds, .combine = 'rbind', .errorhandling = "stop") %dopar% { ##} %do% { 
                           #.combine = list, .errorhandling = "pass", .multicombine = TRUE) %dopar% { #
    
    
    pData <- patData %>%
      dplyr::filter_(~patId == pId, ~chrom == myChrom)
    if (NROW(pData) == 0L) return(NULL)
    
    # these are the selected target regions for the patient where we can expect coverage
    pDataRanges <- GenomicRanges::GRanges(seqnames = myChrom,
                                          ranges = IRanges::IRanges(start=pData$targetStartPos, end=pData$targetEndPos))
    
    
    # get simulated mapping for the specified setting
    pat.n.bamFile <- list.files(path=mappingDir, pattern=paste0("^", getPatStr(pId),"[_]N.s.bam$"), full.names = TRUE)
    pat.t.bamFile <- list.files(path=mappingDir, pattern=paste0("^", getPatStr(pId),"[_]T.s.bam$"), full.names = TRUE)
    
    if ( checkFile(pat.n.bamFile) && checkFile(pat.t.bamFile) ){
      #stopifnot( checkFile(pat.n.bamFile), checkFile(pat.t.bamFile) )
      
      
      # get covered regions in WT and MUT BAM file (GRanges object)
      coveredRegions  <- findCoveredRegionsInTarget(bamFile_WT = pat.n.bamFile, bamFile_MUT = pat.t.bamFile, 
                                                    chrom = myChrom, targetBedFile = myTargetBedFile, extra.margin = 1000L, 
                                                    k = 23L, minWidth = 50L, minCov_WT = 25L, minAvgCov_MUT = 5L) #minGapWidth = , 
      
      logging::loginfo("Found %d covered regions for patient %s.", NROW(coveredRegions), getPatStr(pId))
      
      if (NROW(coveredRegions) > 0L){
        # assess each covered region of the tumor sample with respect to clipped base peaks
        clippedRes_MUT <- foreach::foreach(i=seq_along(coveredRegions), .combine='rbind') %do% {
          #cat("\n ** i=",i, "**\n")
          
          # assess clipped-base peaks in covered region: either 0 or 1 or 2 entries (entry = 3p-5p peak pair) that correspond to highest potential SV
          regionStart <- BiocGenerics::start(coveredRegions)[i]
          regionEnd <- BiocGenerics::end(coveredRegions)[i]
          
          peaks_MUT <- assessClippedBasesInRegion(bamFile = pat.t.bamFile, chrom = myChrom, minThresholdAbs = 3L, minThreshold = 0.05, do.plot=FALSE,
                                                  .targetStartPos = regionStart, .targetEndPos = regionEnd)
          
          cInfo <- peaks_MUT[["clippedInfo"]]
          #pInfo <- peaks_MUT[["peakInfo"]]
          #negCandPos <- peaks_MUT[["nonPeakPos"]]
          
          #cat(cInfo, "\n")
          cInfo
        }#foreach clippedRes
        
        
        ## DEBUG: names
        cat(sprintf("clippedRes_MUT is %s with %d entries and with names %s\n", class(clippedRes_MUT), NROW(clippedRes_MUT), paste(names(clippedRes_MUT), collapse="*")))
        ## DEBUG
        
        stopifnot( NROW(coveredRegions) == NROW(clippedRes_MUT), "nbrPeaks" %in% colnames(clippedRes_MUT) )
        posCandidateRegion.ind <- which(clippedRes_MUT[, "nbrPeaks"]  > 0L)
        negCandidateRegion.ind <- which(clippedRes_MUT[, "nbrPeaks"] == 0L)
        logging::logdebug("Assessed each covered region of the patient regarding clipped peaks: found %d negative and %d positive candidate regions",
                          length(negCandidateRegion.ind), length(posCandidateRegion.ind))
        
        
        # covered regions that overlap the simulated patient regions: for each covered region the patient region index (findOverlaps returns indices of subject)
        pDataHit.ind <- IRanges::findOverlaps(query=coveredRegions, subject = pDataRanges, select="arbitrary")
        logging::logdebug("Covered regions and corresponding patient region index vector with %d NA values (%s).", sum(is.na(pDataHit.ind)), paste(pDataHit.ind, collapse="-"))
        
        # mkuhn, 2016-01-26: I expect that each covered regions has its origin in a simulated patient region
        stopifnot( length(pDataHit.ind) == NROW(clippedRes_MUT) )
        
        # evaluate only covered patient regions: use given margin and always the whole region match
        cbind(patId=pId, chrom=myChrom, covRegStart=GenomicRanges::start(coveredRegions), covRegEnd=GenomicRanges::end(coveredRegions),
              clippedRes_MUT[,c("nbrPeaks", "minScore", "posDiff", "pos.5p", "pos.3p", "pos.left", "pos.right")],
              pDataHit.ind=pDataHit.ind, pData[pDataHit.ind, c("targetStartPos", "targetEndPos", "SVtype", "SVlength", "SVstart")], stringsAsFactors=FALSE) %>% 
          dplyr::mutate_(SVtruth=~factor( SVtype != 'no', levels = c(FALSE, TRUE), labels = c("noSV", "SV")),
                         SVleft =~ifelse(is.na(SVstart) | SVstart == 0L, yes = -1L, no = targetStartPos + SVstart),
                         SVright=~ifelse(SVleft < 0L, yes = -1L, no = targetStartPos + SVstart + SVlength),
                         # screenRegOverlap is T/F if clipped-peak overlaps with simulated SV, otherwise (i.e. if no peak) NA 
                         screenRegOverlap=~pos.left - margin.bp < SVright & pos.right + margin.bp > SVleft,
                         #screen match? "no peak" in screening is match if noSV or this region was not in the simualted patient data  (hence, no SV simulated)
                         screenCall=~factor(ifelse(nbrPeaks > 0L, yes=screenRegOverlap, no=SVtruth == 'noSV' | is.na(SVtruth)), levels = c(FALSE, TRUE), labels = c("noMatch", "Match")),
                         # whole region match
                         screenCallReg=~factor(nbrPeaks > 0 && ! is.na(pos.left) && SVtruth == 'SV' || nbrPeaks == 0L && is.na(pos.left) && SVtruth == 'noSV', levels = c(FALSE, TRUE), labels = c("noMatch", "Match")) )
      } else {
        #else: NROW(coveredRegions) == 0
        NULL
      }
    } else {
      #else: BAM files missing
      logging::logwarn("sorted BAM files for patient %s not available!", pId)
      NULL  
    }
    
  }# hcaerof
  
  logging::loginfo("There are %d evaluations of covered regions.", NROW(evalScreenDat))
  
  if (isTRUE(saveEvaluation)){
    # write out screening performance
    screenPerfDir <- getVirtualPatientPath(what = "performance", .sample.prop=mySampleProp, .ngs.prop=myNGSProp)
    dir.create(screenPerfDir, showWarnings = FALSE, recursive = TRUE)
    screenPerfFilename <- paste0(PERF_PREFIX, "_screen_chr:", myChrom, "_margin:", margin.bp, ".rds")
    
    logging::loginfo("Write out performance of screening to file %s in directory %s.", screenPerfFilename, screenPerfDir)
    saveRDS(evalScreenDat, file.path(screenPerfDir, screenPerfFilename))
  }
  
  return(evalScreenDat)
}





#' Collects all benchmark performance objects based on simulation runs and gathers the performances in a common dataframe. OLD: this function is specific for tcSVM-learner. There is a general one now.
#' 
#' Sample and NGS-information is parsed from path of the benchmark files. In the filename of the benchmark file on finds also information on the proportion of somatic SVs in the training data.
#' @param cov optional filter on only some coverage values from sampling settings in simulations
#' @param margin.bp optional filter on only specific margin values from feature extraction in simulations
#' @return dataframe with benchmark performances from the different available scenarios
#' @note There is a general performance evaluation harvester function. So, this function is obsolete.
#' @export
harvestLearnerPerformanceEvaluations_OLD <- function(baseDir=BASEDIR, cov=NULL, margin.bp=NULL){
  perfBasedir <- file.path(baseDir, "performance_tcsvm")
  
  perfFiles <- list.files(path=perfBasedir, pattern = "benchmark_margin[[:digit:]]+_ssv[[:digit:]]+[.]rds", recursive = TRUE, full.names = TRUE)
  
  if (! is.null(cov)){
    perfFiles <- grep(pattern=paste0("cov(", paste0(cov, collapse="|"), ")"), perfFiles, value = TRUE)
  }
  
  
  # extract only those files that match the given margin values
  if (! is.null(margin.bp)){
    perfFiles <- grep(pattern=paste0("_margin(", paste0(margin.bp, collapse="|"), ")"), perfFiles, value = TRUE)
  }
  
  if (length(perfFiles) == 0L) return(invisible(NULL))
  

  # mkuhn, 2015-12-04 handle parallelism on a package basis, not per function  
  # # register number of cores
  # nbrCores <- max(1L, round(nbrCores))
  # doParallel::registerDoParallel(cores=nbrCores)
  # on.exit(doParallel::stopImplicitCluster())
  
  # load all benchmark objects from the performance/-directory in turn
  allPerformance <- foreach::foreach(pf=perfFiles, .combine = 'rbind') %dopar% {
    # mkuhn, 20150424: parse information from path and file name
    
    # ssv <- as.numeric(sub("benchmark_ssv", "", sub(".rds", "", basename(pf), fixed=TRUE), fixed=TRUE))
    pfbase_split <- strsplit(sub("[.]rds$", "", basename(pf)), split="_", fixed=TRUE)[[1L]]
    stopifnot( length(pfbase_split) == 3L && pfbase_split[1L] == 'benchmark' )
    
    margin <- as.numeric(sub("^margin", "", pfbase_split[2L])) 
    ssv <- as.numeric(sub("^ssv", "", pfbase_split[3L]))
    
    sample.prop <- parseSampleProp(dirname(pf))
    ngs.prop <- parseNGSProp(dirname(pf))
    
    perf <- readRDS(pf)
    cbind(cov=sample.prop$cov, tumorLoad=sample.prop$tumorLoad, 
          pe.ins.mean=ngs.prop$pe.ins.mean, pe.ins.sd=ngs.prop$pe.ins.sd, read.length=ngs.prop$read.length,
          margin.bp=margin, ssv.prop=ssv,
          mlr::getBMRPerformances(perf, as.df = TRUE))
  }
  
  return(allPerformance)
}


#' Collects all performance evaluation information based on simulation runs and gathers the performances in a common dataframe.
#' 
#' Different things can be evaluated, e.g. the ability of the screen to find the SVs, or the different learners.
#' 
#' Sample and NGS-information is parsed from path of the benchmark files.
#' In the filename of the performance evaluation files one typically finds additional information on the properties (e.g. proportion of somatic SVs in the training data or margins used).
#' @param what character that describes what performance evaluation is requested. E.g. "screen" or "learner_tcsvm".
#' @param cov.values numeric vector. optional filter on only some coverage values from sampling settings in simulation.
#' @param tl.values numeric vector. optional filter on tumor load values
#' @param margin.values numeric vector. optional filter on only specific margin values
#' @param keepExtensions flag. Do we tolerate performance files with an extension at all? Defaults to \code{FALSE}.
#' @param extensions character. File name extensions just before .rds in the end of filename. If extension filtering its name:value components are also parsed into the performance dataframe.
#' @return dataframe with performances from the different available scenarios
#' @export
harvestPerformanceEvaluations <- function(baseDir=BASEDIR, what=c("screen", "learner"), cov.values=NULL, tl.values=NULL, margin.values=NULL, tune.values=NULL,
                                          keepExtensions=FALSE, extensions=NULL){ #}, nbrCores=2L){
  what <- match.arg(what)
  
  stopifnot( isTRUE(keepExtensions) || is.null(extensions) )
  
  perfBasedir <- file.path(baseDir, "performance")
  if (! isTRUE(file.exists(perfBasedir)) ) return(invisible(NULL))
  
  stopifnot( is.character(what), length(what) == 1L )
  
  perfFiles <- list.files(path=perfBasedir, pattern = paste0("^", PERF_PREFIX, "_", what, ".*rds$"), recursive = TRUE, full.names = TRUE)
  
  # special filters: cov, TL and margin
  if (! is.null(cov.values) && is.numeric(cov.values) && length(cov.values) >= 1L ){
    perfFiles <- grep(pattern=paste0("cov(", paste0(cov.values, collapse="|"), ")"), perfFiles, value = TRUE)
  }
  
  # Tumor load (TL)
  if (! is.null(tl.values) && is.numeric(tl.values) && length(tl.values) >= 1L ){
    perfFiles <- grep(pattern=paste0("TL(", paste0(tl.values, collapse="|"), ")"), perfFiles, value = TRUE)
  }
  
  # margin: extract only those files that match the given margin values
  if (! is.null(margin.values) && is.numeric(margin.values) && length(margin.values) >= 1L ){
    perfFiles <- grep(pattern=paste0("_margin:(", paste0(margin.values, collapse="|"), ")"), perfFiles, value = TRUE)
  }
  
  # tune.values: extract only those files that match the given margin values
  if (! is.null(tune.values) && is.character(tune.values) && length(tune.values) >= 1L ){
    perfFiles <- grep(pattern=paste0("_tune:(", paste0(tune.values, collapse="|"), ")"), perfFiles, value = TRUE)
  }
  
  # extensions: extract only those files that match the given margin values
  if (! is.null(extensions) && is.character(extensions) && length(extensions) >= 1L ){
    perfFiles <- grep(pattern=paste0("(", paste0(extensions, collapse="|"), ").rds$"), perfFiles, value = TRUE)
  }
  
  
  if (length(perfFiles) == 0L) {
    loginfo("No performance RDS files remained to be harvested..")
    return(invisible(NULL))
  }
  
  
  # what to extract from the screen performance evaluation data.
  # Special functions for each type of performance data (if necessary)
  # len50_120: describes a range of SV-length
  summScreenPerf <- function(perfObj, SVlengthCondition=c("SVall", "SVlen50_120", "SVlen150_250")){
    SVlengthCondition <- match.arg(SVlengthCondition)
    stopifnot( is.data.frame(perfObj), NROW(perfObj) >= 1L, all(c("SVtruth", "screenCall", "screenCallReg", "SVlength") %in% names(perfObj)) )
    
    perfObj <- switch(SVlengthCondition,
                      SVlen50_120=dplyr::filter_(perfObj, ~SVlength == 0L | SVlength >= 50L & SVlength <= 120L),
                      SVlen150_250=dplyr::filter_(perfObj, ~SVlength == 0L | SVlength >= 150L & SVlength <= 250L),
                      SVall=,
                      perfObj)
    
    confTabMargin <- xtabs(~SVtruth + screenCall, data=perfObj)
    stopifnot( identical(rownames(confTabMargin), c("noSV", "SV")) )
    stopifnot( identical(colnames(confTabMargin), c("noMatch", "Match")) )
    confTabRegion <- xtabs(~SVtruth + screenCallReg, data=perfObj)
    stopifnot( identical(colnames(confTabRegion), c("noMatch", "Match")) )
    stopifnot( identical(rownames(confTabRegion), c("noSV", "SV")) )
    
    retDF <- data.frame(nbrPat=length(table(perfObj$patId)), nbrRegEval=sum(confTabMargin), nbrSVReg=sum(perfObj$SVtruth == 'SV'), meanSVlen=mean(perfObj$SVlength[perfObj$SVlength > 0L]),
                        senMargin= prop.table(confTabMargin, 1L)["SV", "Match"], specMargin= prop.table(confTabMargin, 1L)["noSV", "Match"],
                        senRegion= prop.table(confTabRegion, 1L)["SV", "Match"], specRegion= prop.table(confTabRegion, 1L)["noSV", "Match"])
    
    names(retDF) <- paste(names(retDF), SVlengthCondition, sep=".")
    return(retDF)
  }
  
  
  # load all benchmark objects from the performance/-directory in turn
  allPerformance <- foreach::foreach(pf=perfFiles, .combine = 'rbind') %do% {
    # mkuhn, 20150424: parse information from path and file name
    #browser()
    # ssv <- as.numeric(sub("benchmark_ssv", "", sub(".rds", "", basename(pf), fixed=TRUE), fixed=TRUE))
    pfbase_split <- strsplit(sub("[.]rds$", "", basename(pf)), split="_", fixed=TRUE)[[1L]]
    stopifnot( length(pfbase_split) >= 2L, pfbase_split[1L] == PERF_PREFIX, pfbase_split[2L] == what )
    
    # remaining options
    pfbase_split <- pfbase_split[-c(1L,2L)]
    
    hasExtension <- length(pfbase_split) > 2L
    
    if ( ! hasExtension || isTRUE(keepExtensions) ){
      pfopt_split <- strsplit(pfbase_split, split = ":")
      
      # keep only things that have two components as in name:value
      pfopt_split <- pfopt_split[which(lengths(pfopt_split) == 2L)]
      
      stopifnot( all(lengths(pfopt_split) == 2L) )
      
      names(pfopt_split) <- sapply(pfopt_split, head, n=1)
      # names optList (potentially with different data types)
      optDF <- data.frame(lapply(pfopt_split, function(x){ z <- x[2]; if (z == gsub("[^0-9]", "", z)) as.numeric(z) else z }))
      
      
      sample.prop <- parseSampleProp(dirname(pf))
      ngs.prop <- parseNGSProp(dirname(pf))
      
      
      baseInfo <- data.frame(cov=sample.prop$cov, tumorLoad=sample.prop$tumorLoad, 
                             pe.ins.mean=ngs.prop$pe.ins.mean, pe.ins.sd=ngs.prop$pe.ins.sd, read.length=ngs.prop$read.length)
      
      
      # read in RDS-object
      perfObj <- readRDS(pf)
      
      #cat(sprintf("perfObj in %s has class %s\n", pf, class(perfObj)))
    

      perfDF <- switch(what,
                       screen={cbind(summScreenPerf(perfObj), summScreenPerf(perfObj, SVlengthCondition = "SVlen50_120"), summScreenPerf(perfObj, SVlengthCondition = "SVlen150_250"))},
                       learner=if (inherits(perfObj, "BenchmarkResult")) mlr::getBMRPerformances(perfObj, as.df = TRUE) else NULL,
                       stop("Error: read in evaluation for  ", what, " is not implemented here!"))
      
      
      # return
      if (! is.null(perfDF)) cbind(baseInfo, optDF,perfDF) else NULL
    } else NULL
    
  }
  
  return(allPerformance)
}



#' Plot the screening performance.
#' 
#' The performance of all SVs is used (no subgrouping for certain window of SV lengths).
#' @param useCow flag if we want to use the cowplot package
#' @param useRegion flag if whole region should be used for screen performance evaluation (or margin)
#' @param myMargin numeric. Size of margin used for evaluation (only used if \code{useRegion=FALSE})  myMargin=25L, 
#' @export
plotScreenPerformance <- function(scrPerf, minTumorLoad=20L, useRegion=TRUE, doROC=FALSE,
                                  useCow=TRUE){
                                  #imgFilename=NULL, path=NULL){
  
  if (isTRUE(useRegion)){
    senCol <- "senRegion.SVall"
    specCol <- "specRegion.SVall"
  } else {
    senCol <- "senMargin.SVall"
    specCol <- "specMargin.SVall"
  }
  
  
  dataPrep <- scrPerf %>% 
    #dplyr::rename_(Sensitivity=~senMargin.SVall, Specificity=~specMargin.SVall) %>% 
    dplyr::rename_(Sensitivity=senCol, Specificity=specCol) %>% 
    dplyr::filter_(~tumorLoad > minTumorLoad) %>% 
    dplyr::select_(~cov, ~tumorLoad, ~Sensitivity, ~Specificity) 
  
  
  # further filtering: for region metric we need only one entry
  if (isTRUE(useRegion)){
    dataPrep <- unique(dataPrep)
      #dplyr::select_(~ - margin) %>%
  } 
    
  
  plotObj <- if (isTRUE(doROC)){
    # # ROC curve (with varying tumorload)
    ggplot(data = dataPrep, aes(x=1-Specificity, y=Sensitivity, shape=factor(cov))) +
      geom_point(aes(col=tumorLoad)) + #aes(size=1.1)
      geom_smooth(aes(linetype=factor(cov)), se = FALSE, span = 1.5, colour = "grey") + #method="lm", formula = y ~ splines::ns(x, 2L) ) +
      #scale_color_hue(name="Coverage") +
      scale_colour_gradient() + 
      labs(title="Screen Performance", x="1-Specificity", y = "Sensitivity") +
      if (isTRUE(useCow) && isTRUE(requireNamespace("cowplot"))){
        # cowplot
        background_grid() 
      } else {
        theme(title = element_text(size = rel(1.5)),
              strip.text = element_text(size = rel(1.2)),
              axis.text = element_text(size = rel(1.2)))
      }
  } else {
    # SEN and SPEC plots
    dataPrep %>% 
      tidyr::gather_(key_col = "measure", value_col = "value", gather_cols = c("Sensitivity", "Specificity")) %>% 
      ggplot(aes(x=tumorLoad, y=value, col=factor(cov))) + 
      geom_point() + #aes(size=1.1)
      stat_smooth(se = FALSE, method="loess", span=1) + #method="lm", formula = y ~ splines::ns(x, 2L) ) +
      scale_color_hue(name="Coverage") + 
      facet_wrap(~measure) + 
      labs(title="Screen Performance", x="Tumor Load", y = "") + 
      if (isTRUE(useCow) && isTRUE(requireNamespace("cowplot"))){
        # cowplot
        background_grid() 
      } else {
        theme(title = element_text(size = rel(1.5)),
              strip.text = element_text(size = rel(1.2)),
              axis.text = element_text(size = rel(1.2)))
      }
  }
  
  
  
  # if (! is.null(imgFilename) && is.character(imgFilename) ){
  #   #ggsave(filename = imgFilename, path = path)
  #   cowplot::save_plot(filename = file.path(path, imgFilename), plot = plotObj)
  # }
  
  plotObj
}



#' Plot learner performance dataframe
#' @param allPerf dataframe with performance values
#' @param useCow flag if we want to use the cowplot package
#' @export
plotLearnerPerformance <- function(allPerf, useCow=TRUE) { #}, imgFilename=NULL, path=NULL){
  
  plotObj <- allPerf %>%
    #dplyr::filter_(~tumorLoad > 20) %>%
    dplyr::select_((~-mcc)) %>% 
    dplyr::rename_(Sensitivity=~tpr, Accuracy=~acc, MCC=~my.mcc) %>% 
    dplyr::mutate_(learner.id=~stringr::str_replace_all(learner.id, pattern = c("classif." = "", ".filtered.tuned" = ""))) %>% 
    tidyr::gather_(key_col = "measure", value_col = "value", gather_cols = c("Sensitivity", "Accuracy", "MCC")) %>% 
    ggplot(aes(x=tumorLoad, y=value, col=factor(cov))) + 
    geom_point() + 
    stat_smooth(se = FALSE, method="lm", formula = y ~ splines::ns(x,2L) ) +  #method = "loess") +
    coord_cartesian(ylim=c(0,1)) + 
    scale_color_hue(name="Coverage") + 
    facet_grid(measure~learner.id) + 
    labs(title="Algorithm (Screen + Confirmation) Performance", x="Tumor Load", y = "") + 
    if (isTRUE(useCow) && isTRUE(requireNamespace("cowplot"))){
      # cowplot
      background_grid() 
    } else {
      #plain ggplot2 (w/o cowplot)
      theme(title = element_text(size = rel(1.2)),
            legend.text = element_text(size = rel(1.2)),
            strip.text = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)))
    }
  
  # if (! is.null(imgFilename) && is.character(imgFilename) ){
  #   #ggsave(filename = imgFilename, path = path)
  #   cowplot::save_plot(filename = file.path(path, imgFilename), plot = plotObj)
  # }
  
  plotObj
  
}
