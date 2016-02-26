



library(mlr)
po <- readRDS("~/SVdata/performance/cov60_TL90/PEm200sd25_R100/perf_learner_margin:25_tune:random_MO_v0.6.1_r20.rds")

# lsit of benchmark results
po
class(po)
sapply(po, class)

tuneRes <- mlr::getBMRTuneResults(po)[["somaticSV"]][["classif.svm.filtered.tuned"]][[1]]

# get tuned learner
svmTuned <- mlr::getBMRLearners(po)[["classif.svm.filtered.tuned"]] %>% 
  setHyperPars(par.vals = tuneRes$x)





# if po is list of bmrResults
# bmark result
bmrKnn <- po[[3]]
class(bmrKnn)

tuneRes <- mlr::getBMRTuneResults(bmrKnn)[[1]][[1]][[1]]
mlr::getBMRAggrPerformances(bmrKnn)[[1]][[1]]
mlr::getBMRPerformances(bmrKnn)[[1]][[1]]
getHyperPars(mlr::getBMRLearners(bmrKnn)[[1]])
getParamSet(mlr::getBMRLearners(bmrKnn)[[1]])
mlr::getBMRLearners(bmrKnn)[[1]] %>% 
  setHyperPars(par.vals = tuneRes$x)


bmrKnn




