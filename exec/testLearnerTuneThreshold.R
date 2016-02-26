

library(mlr)
set.seed(1909)

svmod::setupParallel(cpus = 13L)

tuneMode <- "grid"

logreg_Learner <- mlr::makeLearner("classif.logreg", config=list(on.learner.error="warn"))
logreg_Learner.f <- mlr::makeFilterWrapper(logreg_Learner, fw.method = "rf.importance") %>% 
  mlr::setPredictType(predict.type = "prob")
#train(logreg_Learner.f)
tuning.ps.logreg <-  ParamHelpers::makeParamSet(
  ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
)


rf_Learner <- mlr::makeLearner("classif.randomForest", config=list(on.learner.error="warn"))
rf_Learner.f <- mlr::makeFilterWrapper(rf_Learner, fw.method = "rf.importance") %>% 
  mlr::setPredictType(predict.type = "prob")
tuning.ps.rf <-  ParamHelpers::makeParamSet(
  ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
  , ParamHelpers::makeIntegerParam("ntree", lower = 5, upper = 10, trafo = function(x) 2^x)
  , ParamHelpers::makeIntegerParam("mtry", lower = 5, upper = 15)
)


svm2_Learner <- mlr::makeLearner("classif.svm", type="nu-classification", config=list(on.learner.error="warn"), predict.type = "prob")
svm2_Learner.f <- mlr::makeFilterWrapper(svm2_Learner, fw.method = "rf.importance") %>% 
  mlr::setPredictType(predict.type = "prob")
tuning.ps.svm2 <-  ParamHelpers::makeParamSet(
  ParamHelpers::makeNumericParam("fw.perc", lower=0.66, upper = 1)
  , ParamHelpers::makeNumericParam("nu", lower=0.05, upper = 0.75)
  #, ParamHelpers::makeNumericParam("cost", lower=0.5, upper = 15)
  , ParamHelpers::makeDiscreteParam("kernel", values = c("linear", "radial", "polynomial"))
  , ParamHelpers::makeNumericParam("gamma", lower=-8, upper = 8, requires = quote(kernel != 'linear'), trafo = function(x) 2^x)
  , ParamHelpers::makeIntegerParam("degree", lower = 2L, upper = 3L, requires = quote(kernel=="polynomial"))
)


tuneCtrl.tuneTh <- switch(tuneMode,
                          random = mlr::makeTuneControlRandom(maxit = 50L, tune.threshold = TRUE),
                          gensa = mlr::makeTuneControlGenSA(maxit = 50L, tune.threshold = TRUE),
                          irace = mlr::makeTuneControlIrace(maxExperiments = 100L, tune.threshold = TRUE),
                          grid = ,
                          mlr::makeTuneControlGrid(resolution = 2L, tune.threshold = TRUE) #c(fw.perc=2L, cost=4L, gamma=3L))
)


rDesc <- makeResampleDesc("Holdout")
k <- 3L
resampDesc.inner <-  mlr::makeResampleDesc("Subsample", split=.75, iters = k, stratify = TRUE)
logreg_Learner.f.t <- mlr::makeTuneWrapper(logreg_Learner.f, rDesc, measures = list(my.mcc),
                                           par.set = tuning.ps.logreg, control = tuneCtrl.tuneTh, show.info=verbose)
rf_Learner.f.t <- mlr::makeTuneWrapper(rf_Learner.f, resampling = resampDesc.inner, measures = list(my.mcc),
                                       par.set = tuning.ps.rf, control = tuneCtrl.tuneTh, show.info=verbose)

svm2_Learner.f.t <- mlr::makeTuneWrapper(svm2_Learner.f, resampling = resampDesc.inner, measures = list(my.mcc),
                                         par.set = tuning.ps.svm2, control = tuneCtrl.tuneTh, show.info=verbose)


trained <- mlr::train(logreg_Learner.f.t, task=somSV_Task, subset = trainInd)
tune_rf <- tuneParams(rf_Learner.f, task = somSV_Task, resampling = rDesc, measures = list(my.mcc), par.set = tuning.ps.rf, control = tuneCtrl.tuneTh)
trained_rf <- mlr::train(rf_Learner.f.t, task=somSV_Task, subset = trainInd)
trained_svm <- mlr::train(svm2_Learner.f.t, task = somSV_Task, subset = trainInd)


# performance

trainedL <- trained_svm #trained_rf

predObj0 <- predict(trainedL, task = somSV_Task, subset = trainInd)
predObj1 <- predict(trainedL, task = somSV_Task, subset = testInd)
#plotLearnerPrediction(svm2_Learner.f.t, task=somSV_Task) #takes time again!
mlr::getConfMatrix(predObj0); mlr::getConfMatrix(predObj1)
table( mlr::getPredictionTruth(predObj1), mlr::getPredictionResponse(predObj1) )
mlr::performance(predObj0, measures = list(my.mcc, mlr::mcc, mlr::acc, mlr::tpr, mlr::tnr))
mlr::performance(predObj1, measures = list(my.mcc, mlr::mcc, mlr::acc, mlr::tpr, mlr::tnr))
