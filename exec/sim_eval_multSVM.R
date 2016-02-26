# mkuhn, 2015-03-02
# Quick evaluation of multiclass SVM.

library(mlr)

#featDat <- 

featDatY <- cbind(Y=SVtype, featDat)


diffSV_Task <- makeClassifTask(id = "diffSV", data = featDatY, target = "Y")


multSVM_Learner <- makeLearner("classif.svm") #, predict.type = "response")

mlr::train(multSVM_Learner, diffSV_Task)

# tuning parameters
tuning.ps <-  makeParamSet(
  makeDiscreteParam("nu", values = c(0.05, 0.1, 0.25, 0.5)),
  makeDiscreteParam("gamma", values = 2^(-3:-2))
)

tuneCtrl <-  makeTuneControlGrid()

resampDesc.inner <-  makeResampleDesc("CV", iters = 3L)
resampDesc.outer <-  makeResampleDesc("Holdout")

res <- tuneParams(multSVM_Learner, task = diffSV_Task, resampling = resampDesc.inner, par.set = tuning.ps,
                 control = tuneCtrl, measures = list(acc, setAggregation(acc, test.sd)), show.info = FALSE)

multSVM_Learner.t <- makeTuneWrapper(multSVM_Learner, resampDesc.inner, measures = list(acc), par.set = tuning.ps, control = tuneCtrl)
r <- resample(multSVM_Learner.t, diffSV_Task, resampling = resampDesc.outer, show.info = FALSE, measures = list(acc, mmce))
r$measures.test
