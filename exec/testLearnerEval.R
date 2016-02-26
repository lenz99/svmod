# mkuhn, 2016-01-07
# Quick harvest of performance evaluations

library(mlr)
perfDir <- "~/SVdata/performance/cov60_TL90/PEm200sd25_R100/"

p_rand <- readRDS(file.path(perfDir, "perf_learner_margin:25_tune:random_SV2_r300.rds"))

sapply(p_rand, class)
perf_l <- p_rand[[2]]
class(perf_l)

getConfMatrix(getBMRPredictions(perf_l)[[1]][[1]])
getBMRTuneResults(perf_l, as.df = TRUE)
getBMRFilteredFeatures(perf_l, as.df= TRUE)

perf_rand <- lapply(p_rand, function(x) if (inherits(x, "BenchmarkResult")) getBMRPerformances(x, as.df = TRUE) else NULL)
do.call(rbind, perf_rand)
tune_rand <- lapply(p_rand, function(x) if (inherits(x, "BenchmarkResult")) getBMRTuneResults(x, as.df = TRUE) else NULL)
do.call(rbind, tune_rand)


# older SV-code (2016-01-06)
p_rand06 <- readRDS("perf_learner_margin:25_tune:random.rds_20160106") #perf_learner_margin:25_tune:random.rds")

sapply(p_rand06, class)

perf_rand06 <- lapply(p_rand06, function(x) if (inherits(x, "BenchmarkResult")) getBMRPerformances(x, as.df = TRUE) else NULL)
do.call(rbind, perf_rand06)
tune_rand06 <- lapply(p_rand06, function(x) if (inherits(x, "BenchmarkResult")) getBMRTuneResults(x, as.df = TRUE) else NULL)
do.call(rbind, tune_rand)



# use tuned parameters for new SV-code features
kknnTuneRes06 <- getBMRTuneResults(bmr = p_rand06[[3L]])[[1]][[1]][[1]]  #tune_rand06 fw=0.7312252 k=3 epanechnikov
kknnL <- mlr::makeLearner("classif.kknn", scale=TRUE, config=list(on.learner.error="warn")) %>% 
  mlr::makeFilterWrapper(fw.method = "cforest.importance") %>% 
  mlr::setPredictType(predict.type = "prob") %>% 
  setHyperPars(par.vals=kknnTuneRes06$x)

# resample( with fixed holdout test data)
# now, train this learner
#train(knnL, task = , subset = )
#predict() on test data
# and then evaluate this prediction