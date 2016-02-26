#DO_PARALLEL <- TRUE #FALSE

#' which mode of parallelism?
#' parallelization mode: multicore or sockets on localhost?
PCLUST_OBJ <- NULL  # parallel::makePSOCKcluster(names = 3L)



#' Set the parallel computation.
#' 
#' This function can be called at startup or during work with the package.
#' @param myLevel character string on which level to do parallelization. If \code{NA} any regisitered level will run in paralllel.
#' @export
setupParallel <- function(cpus=2L, myLevel=NA){
  stopifnot( is.numeric(cpus), cpus > 0L )
  stopifnot( exists("PCLUST_OBJ") )
  
  if (is.null(myLevel)) myLevel <- NA
  
  
  # multicore on localhost
  if ( is.null(PCLUST_OBJ) ){
    
    # check current settings
    pm.settings <- parallelMap::parallelGetOptions()[["settings"]]
    
    # Check for identical settings: only if not identical I need to restart parallel mode
    if (! isTRUE(pm.settings$mode == 'multicore' && pm.settings$cpus == cpus &&
                 (is.na(pm.settings$level) && is.na(myLevel) || pm.settings$level == myLevel)) ){
      
      
      try(parallelMap::parallelStop(), silent = TRUE)
      
      if (! is.na(myLevel)){
        logging::loginfo("Set parallel mode to multicore with %d cpus and level %s", cpus, myLevel)
        parallelMap::parallelStartMulticore(cpus = cpus, level = myLevel)
      }else {
        logging::loginfo("Set parallel mode to multicore with %d cpus", cpus)
        parallelMap::parallelStartMulticore(cpus = cpus)
      }
    }# fi
    
    # for foreach-looping
    doParallel::registerDoParallel(cores = cpus)
  } else {
    #socket on localhost
    loginfo("Using parallelism via sockets (on localhost)")
    parallelMap::parallelStartSocket(cpus = cpus, level = myLevel)
    doParallel::registerDoParallel(PCLUST_OBJ)
  }
  
  return(invisible(NULL))
}


# mkuhn, 2015-04-24: 
# fix parallel package
.onLoad <- function(libname, pkgname){
  packageStartupMessage("svmod package loaded - register parallelMap levels of mlr and mlr-options") # and register doParallel for foreach-loop")
  
  
  # setting mlr-options
  op <- options()
  op.mlr <- list(
    mlr.show.info=TRUE,
    mlr.on.learner.error="stop",
    mlr.on.learner.warning="warn",
    mlr.on.par.out.of.bounds="stop",
    mlr.on.par.without.desc="stop",
    mlr.show.learner.output=TRUE)
  
  toset <- !(names(op.mlr) %in% names(op))
  if ( any(toset) ) options(op.mlr[toset])
  
  
  # register the mlr-levels so that parallelization happens in the mlr-calls
  parallelMap::parallelRegisterLevels(package="mlr", levels=c("benchmark", "resample", "selectFeatures", "tuneParams"))
  
  
  # default core number
  nbrAvailableCores <- parallel::detectCores()
  nbrCoresToUse.default <- min(nbrAvailableCores, ceiling(sqrt(nbrAvailableCores)+1L))
  
  # defaults to multicore with 3 cpus
  setupParallel(cpus = nbrCoresToUse.default)
  
  invisible(NULL)
}


.onUnload <- function(libpath){
  cat("Clean up parallelMap parallel setup [and doParallel implicit cluster].\n")
  
  if (! is.null(PCLUST_OBJ)) parallel::stopCluster(cl = PCLUST_OBJ)
  doParallel::stopImplicitCluster()
  parallelMap::parallelStop()
}

