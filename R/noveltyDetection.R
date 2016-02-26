
#' Wrapper function to one-class SVM.
#' 
#' Given a feature data matrix \code{X} a one-class SVM is learnt on the data.
#' A one-class SVM for novelty detection predicts \code{TRUE} for values belonging to population and \code{FALSE} for novelty points.
#' 
#' The parameter \eqn{\nu} is an upper bound on the fraction of outliers in the population, i.e. smaller nu = fewer outliers allowed = more complicted boundary to adapt to data, higher variance
#'
#' @author mkuhn, 20140625 
#' @param X matrix of feature data. One row per observational unit. One column per feature.
#' @param kernel character indicating the kernel to be used in the SVM
#' @param gamma parameter for the kernel
#' @param nu numeric parameter value for nu-SVM. In one-class SVM it is an upper bound for the proportion of novelties in the data.
#' @export
applyOneClassSVM <- function(X, nu=0.5, kernel="radial", gamma=0.1, ...){
  stopifnot( is.numeric(X) && is.matrix(X) )
  
  svmMod <- svm(x=X, type="one-classification", scale=TRUE, kernel=kernel, nu=nu, gamma=gamma, ...) #probability=TRUE #not supported for one-class SVM  (yet) 
}


#' Evaluates a one-class SVM model.
#' 
#' A one-class SVM for novelty detection predicts \code{TRUE} for values belonging to population and \code{FALSE} for novelty points (=outliers).
#' If for each observation it is known if it is considered as an outlier (for instance, in a simulation)
#' then we can build the confusion matrix and estimate how well the SVM performed on the data.
#' Either \code{X=} or \code{newdata=} needs to be given. If both are given then \code{newdata=} is used.
#' 
#' @param svmMod a one-class SVM object
#' @param X the training data of \code{svmMod}, i.e. the feature matrix on which the SVM was trained
#' @param newdata if not null the evaluation of SVM is done on this new data set as validation data. Defaults to NULL, i.e. by default the evaluation uses the training data.
#' @param y 0-1 coded vector stating if observation is novelty/outlier. If newdata is given it applies to the validation data, otherwise it refers to the training data \code{X}
#' @param do.plot character flag. 'no'=no plot, '2d'=2-dim plot, '3d'=3-dim plot
#' @param tunedFor optional character: what was this one-class SVM tuned for?
#' @export
evaluateOneClassSVM <- function(svmMod, X=NULL, newdata=NULL, y, do.plot=c("none", "2d", "3d"), tunedFor=""){
  do.plot <- match.arg(do.plot)

  stopifnot( inherits(svmMod, "svm") && svmMod$type == 2 )
  
  stopifnot( is.numeric(y) && is.vector(y) && all( y == 0 | y == 1) )
  stopifnot(! (is.null(X) && is.null(newdata)) ) # one of X or newdata is given
  stopifnot( is.null(X) || NROW(X) == length(resid(svmMod)) )
  stopifnot( is.null(X) ||  NROW(X) == length(y) )
  stopifnot( is.null(newdata) || NROW(newdata) == length(y))
  
  if ( sum(y) > 0.5 * length(y) ) logwarn("There are more 'outliers' (y=1) than normal observations (y=0)!")
  
  print(svmMod)
  
  isTrainingData <- is.null(newdata)
  
  if ( isTRUE(isTrainingData) ) {
    loginfo("Evaluating SVM model on training data (n=%d) with %d novelties.", NROW(X), sum(y))
    evalDat <- X
  } else {
    loginfo("Evaluating SVM model on validation data (n=%d) with %d novelties.", NROW(newdata), sum(y))
    evalDat <- newdata
  }
  
  # svmPredY: prediction of SVM coded as 0=IN and 1=OUTL
  svmPredOutl <- as.numeric(! predict(svmMod, evalDat))
  confTab <- table(Outl=y, svmOutl=svmPredOutl)
  print(confTab)
  cat("\n Row percent:\n")
  print(round(prop.table(confTab, 1), 3))
  
  #acc.overall <- (confTab[1,2] + confTab[2,1]) / sum(confTab)
  acc.overall <- sum(diag(confTab)) / sum(confTab)
  
  # precision (=PPV): proportion of true hits in predicted hits
  prec <- confTab[2,2] / sum(confTab[,2])
  
  # recall (=SEN)
  recall <- confTab[2,2] / sum(confTab[2,])
  
  # SPEC
  spec <- confTab[1,1] / sum(confTab[1,])
  
  #F1-statistic (ignores true negative rate, though)
  f1 <- 2 * prec * recall / (prec + recall)
  
  cat(sprintf("- Overall accuracy: %.2f\n", acc.overall))
  cat(sprintf("- Precision (PPV): %.2f\n", prec))
  cat(sprintf("- Recall (SEN): %.2f\n", recall))
  cat(sprintf("- Specificity (SPEC): %.2f\n", spec))
  cat(sprintf("- F1-score: %.2f\n", f1))
  
  perfMeasures <- list(accuracy=acc.overall, precision=prec, recall=recall, specificity=spec, f1=f1)
  
  
  if ( ! do.plot %in% c("3d", "2d")){
    logdebug("No plot requested for SVM evaluation.")
    mdsCoords <- NULL
  } else {
    # check if any feature is scaled
    if ( all(svmMod$scaled) ){
      evalDatDist <- dist(scale(evalDat, center=svmMod$x.scale[[1]], scale=svmMod$x.scale[[2]]))
    } else {
      logwarn("Data to calculate distance is *NOT* scaled as not all features are scaled in SVM. [There were %d scaled features.]", sum(which(svmMod$scaled)))
      evalDatDist <- dist(evalDat)
    }
    
    titleStr <- if (isTRUE(isTrainingData)) "Model on Training Data" else "Model on Validation Data"
    if (nzchar(tunedFor)) titleStr <- paste0(titleStr,"\nTuned for ", tunedFor)
    
    # visualize the result of SVM:
    # SVM predicts FALSE for novelty points
    if ( do.plot == "3d" && requireNamespace("scatterplot3d", quietly = TRUE) ){
      mdsCoords <- cmdscale(evalDatDist, k = 3)
      scatterplot3d::scatterplot3d(mdsCoords, color=1+y, pch=ifelse(predict(svmMod, evalDat), 1, 19)) #pch=ifelse(1:NROW(featureDat) %in% svmMod$index, 19, 1)) #Support Vectors
      title(main = titleStr)
      legend("topleft", legend = c("Normal", "DSV"), fill = 1:2, xpd=TRUE, inset = c(0, -.1), cex = .8)
      legend("bottomright", legend = c("Normal", "Novelty"), title = "Prediction", pch=c(1, 19), col="grey", cex = .8,
             xpd=TRUE, inset = c(0, -.15))
    } else if (do.plot != 'none'){ #== '2d'){
      mdsCoords <- cmdscale(evalDatDist, k = 2L)
      plot(mdsCoords, col=1+y, pch=ifelse(predict(svmMod, evalDat), 1, 19)) #pch=ifelse(1:NROW(featureDat) %in% svmMod$index, 19, 1)) #Support Vectors
      title(main = titleStr)
      legend("topleft", legend = c("Normal", "DSV"), fill = 1:2, cex = .8)
      legend("bottomright", legend = c("Normal", "Novelty"), title = "Prediction", pch=c(1, 19), col="grey", cex = .8)
    }  
  }#fi do.plot
  
  return(list(conf=confTab, perf=perfMeasures, mdsCoords=mdsCoords))
  
}


#' Tune a one-class SVM given a training data set and status (IN/OUT).
#' 
#' This function allows to tune hyper-parameters of a one-class SVM by means of cross-validation or bootstrapping or 
#' by external validation (hold out) data. The hyper-parameters are varied on a specified grid and 
#' one cross-validation partition or bootstrap split is used for all hyper-parameters.
#' Error is calculated on the external validation data or on the out-of-bag data in CV or bootstrap.
#' It is OK to select the best parameters based on these errors but in order to generalize the error (and not over-optimistic low error)
#' it is best to use external validation or an outer CV loop for the selected model (i.e. selected optimal parameter values).
#' 
#' The CV is using the same proportion of IN/OUT as the full training data [mkuhn, 20150302]
#' 
#' Question is if bootstrapping is disturbing the one-class SVM, because we have so many duplicates in the bootstrapped training data.
#' 
#' The function is based on the \code{\link[e1071]{tune}}-function.
#' @param train.x formula or matrix with predictors
#' @param train.status logical vector that should be coded like what \code{\link[e1071]{predict}} does at one-class SVM, i.e. IN=\code{TRUE}, OUT=\code{FALSE}.
#' @param ranges list with parameter ranges for tuning.
#' @param validation.x a separate validation data set. If specified this validation set will be used to evaluate the models (with varying parameters from the grid). When non-NULL, sampling is set to 'fix'
#' @param validation.status logical vector that should be coded like what \code{\link[e1071]{predict}} does at one-class SVM, IN=\code{TRUE}, OUT=\code{FALSE}.
#' @param tune.measure specifies the evaluation-measure of the confusion table (outlier: TRUE/FALSE vs predicted outlier: TRUE/FALSE) in the validation/test rounds. The error of the confusion table is defined as \code{1-measure}.
#' @export
ocSVMTune <- function(train.x, train.status, data=list(), validation.x = NULL, validation.status = NULL, ranges = NULL,
                      predict.func=predict, tune.measure=c("accuracy", "precision", "recall", "specificity", "f1"),
                      tunecontrol = tune.control(), ...){
  call <- match.call()
  tune.measure <- match.arg(tune.measure)
  
  ##resp <- function(formula, data) { model.response(model.frame(formula, data)) }
  
  # rows: prediction, columns: status
  classAgreement <- function(tab, measure=c("accuracy", "precision", "recall", "specificity", "f1")) {
    measure <- match.arg(measure)
    n <- sum(tab)

    coln <- colnames(tab)
    rown <- rownames(tab)
    stopifnot( identical(c("outlStatus", "outlPred"), names(dimnames(tab))) )
    stopifnot( !is.null(dimnames(tab)) && length(intersect(coln,rown))>0 )
    stopifnot( identical(c("FALSE", "TRUE"), rown) )
    
    # catch cases where the confusion matrix is incomplete (not 2x2)
    prec   <- if (! "TRUE" %in% coln) NA else if (! "TRUE" %in% rown) 0 else tab["TRUE","TRUE"] / sum(tab[,"TRUE"]) #PPV
    recall <- if (! "TRUE" %in% coln) 0 else if (! "TRUE" %in% rown) NA else tab["TRUE","TRUE"] / sum(tab["TRUE",]) #SEN
    spec <- if (! "FALSE" %in% coln) 0 else if (! "FALSE" %in% rown) NA else tab["FALSE","FALSE"] / sum(tab["FALSE",])
    
    p0 <- switch(measure,
                 accuracy=sum(diag(tab))/n, #diag on table (even with only one column, k x 1) seems to work
                 precision=prec,
                 recall=recall,
                 specificity=spec,
                 f1= 2 * prec * recall / (prec + recall)
    )
    # debug line
    #cat("p0=",p0)
    p0
  }
  

  
  
  useFormula <- inherits(train.x, "formula") #is.null(train.y)
  if (useFormula && (is.null(data) || length(data) == 0)) data <- model.frame(train.x)
  if (is.vector(train.x)) train.x <- t(t(train.x)) # vector as nx1-matrix
  ## if (is.data.frame(train.y)) train.y <- as.matrix(train.y)
  if (!is.null(validation.x)){ 
    #mkuhn, 20150302: if validation set is given I use 100% of training data for learning.
    tunecontrol$sampling <- "fix"
    tunecontrol$fix <- 1 #also boot.size=1 #?!
  }
  
  # PREPARE SAMPLING (training/test data splits)
  n <- nrow(if (useFormula) data else train.x)
  perm.ind <- sample(n)
  
  # stratified sampling that respects the proportion of IN and OUT in the full training data
  ind.norm <- which(  train.status)
  ind.outl <- which(! train.status)
  n.norm <- length(ind.norm)
  n.outl <- length(ind.outl)
  stopifnot( n.norm + n.outl == n )
  
  perm.ind.norm <- sample(ind.norm)
  perm.ind.outl <- sample(ind.outl)
  
  if (tunecontrol$sampling == "cross") {
    if (tunecontrol$cross > n) 
      stop(sQuote("cross"), " must not exceed sampling size!")
    if (tunecontrol$cross == 1) 
      stop(sQuote("cross"), " must be greater than 1!")
  }
  
  # perm.ind is the initial permutation of the training data.
  # Resampling is stratified for IN/OUT [mkuhn, 20150302]
  # train.ind is a 
  # + list of a single subset of the permutation [fix], 
  # + list of subsets of the permutation [cross]
  # + list of permutations [boot] 
  train.ind <- if (tunecontrol$sampling == "cross") {
    #tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x]) # get k=tunecontrol$cross subsets from the initial permuation perm.ind
    cv.ind.norm <- tapply(1:n.norm, cut(1:n.norm, breaks = tunecontrol$cross), function(x) perm.ind.norm[-x])
    cv.ind.outl <- tapply(1:n.outl, cut(1:n.outl, breaks = tunecontrol$cross), function(x) perm.ind.outl[-x])
    cv.ind <- vector("list", length=tunecontrol$cross)
    for (i in 1:tunecontrol$cross){
      cv.ind[[i]] <- c(cv.ind.norm[[i]], cv.ind.outl[[i]])
    }
    cv.ind
  } else if (tunecontrol$sampling == "fix"){ 
    #list(perm.ind[1:trunc(n * tunecontrol$fix)])
    fix.ind.norm <- perm.ind.norm[1:trunc(n.norm * tunecontrol$fix)]
    fix.ind.outl <- perm.ind.outl[1:trunc(n.outl * tunecontrol$fix)]
    list(c(fix.ind.norm, fix.ind.outl))
  } else {
    #lapply(1:tunecontrol$nboot, function(x) sample(n, n * tunecontrol$boot.size, replace = TRUE))  #bootstrap with own sample()-command with replacement.
    # bootstrap sample of the normal observations and the outlying observations
    lapply(1:tunecontrol$nboot, function(x) 
      c(sample(perm.ind.norm, size = length(perm.ind.norm) * tunecontrol$boot.size, replace = TRUE),
        sample(perm.ind.outl, size = length(perm.ind.outl) * tunecontrol$boot.size, replace = TRUE) )
    )
  }
  
  #mkuhn, 20140701: debug output
  cat( sprintf("Length of each sampling from training data: %s. [MK]\n", paste(sapply(train.ind, length), collapse = ", ")) )
  
  # Start with parameter-grid
  parameters <- if (is.null(ranges)) data.frame(dummyparameter = 0) else expand.grid(ranges)
  p <- nrow(parameters)
  
  # mkuhn, 20150302: with random one assess only some random parameter combinations
  if (!is.logical(tunecontrol$random)) {
    if (tunecontrol$random < 1) 
      stop("random must be a strictly positive integer")
    if (tunecontrol$random > p) tunecontrol$random <- p
    parameters <- parameters[sample(1:p, tunecontrol$random), ]
  }
  
  
  # mkuhn, 20140701: init vectors with correct length
  model.variances <- model.errors <- numeric(p) #c()
  for (para.set in 1:p) {
    sampling.errors <- numeric(length(train.ind)) #c()
    for (sample in 1:length(train.ind)) {
      repeat.errors <- numeric(tunecontrol$nrepeat) #c()
      for (reps in 1:tunecontrol$nrepeat) { # repeats for random algorithms?
        
        pars <- if (is.null(ranges)) NULL else lapply(parameters[para.set, , drop = FALSE], unlist)
        
        
        ### TRAINING data:
        currentTrain.ind <- train.ind[[sample]]
        
        #train model on training data
        model <- if (useFormula) 
          do.call(svm, c(list(train.x, data = data, subset = currentTrain.ind), 
                         pars, list(...)))
        else do.call(svm, c(list(train.x[currentTrain.ind, ], y = NULL), ##y = train.y[currentTrain.ind]),
                            pars, list(...)) )
        
        # mkuhn, 20140701: debug output of model
        #cat(sprintf("Trained SVM-model of type %d on sample of training data, %d rows. [MK]\n", model$type, NROW(train.x[currentTrain.ind,])))
        #cat("number of unique indices in training sample: ",length(unique(currentTrain.ind)), " [MK]\n")
  
        
        
        ### TEST data
        
        # PREDICTION
        # model prediction on full validation-data OR remaining test-data
        pred <- predict.func(model, if (!is.null(validation.x)) 
          validation.x
          else if (useFormula) 
            data[-currentTrain.ind, , drop = FALSE]
          else if (inherits(train.x, "matrix.csr")) 
            train.x[-currentTrain.ind, ]
          else train.x[-currentTrain.ind, , drop = FALSE])
        
        ## mkuhn, 20140701: debug output for pred
        #cat(sprintf("Dimension of test data: %s. [MK]\n", paste(dim(train.x[-currentTrain.ind, , drop = FALSE]), collapse="-")))
        
        # mkuhn, 20140701: debug output of pred
        #print(table(pred))
        
        
        # OBSERVED
        # observed, i.e. true STATUS in full validation-data or remaining test
        true.y <- if (!is.null(validation.status))
          validation.status
        else train.status[-currentTrain.ind]
        
        
        # mkuhn, 20140701: debug output of true.y
        #print(table(true.y))
        
        
        # ERROR
        # calculate the error in this repetition:
        # with own error-function (if available) or defaults to
        # accuracy [classification] or mean squared error (MSE) [regression]
        stopifnot( is.logical(true.y) )
        repeat.errors[reps] <- if (!is.null(tunecontrol$error.fun)) 
          tunecontrol$error.fun(true.y, pred)
        else if ((is.logical(true.y) || is.factor(true.y)) && 
                   (is.logical(pred) || is.factor(pred) || is.character(pred))) {
          confTab <- table(outlStatus=(! true.y), outlPred=(! pred))
          #print(confTab)
          1 - classAgreement(confTab, measure=tune.measure)
        } 
        # mean squared error (MSE) for continuous predictors
        #else if (is.numeric(true.y) && is.numeric(pred)) crossprod(pred - true.y)/length(pred)
        else stop("Dependent STATUS variable is not logical: TRUE=IN, FALSE=OUT.")
      }#rof reps
      
      # if necessary aggregate the errors in repetition (i.e. runs with same parameters and same training/test sample split)
      sampling.errors[sample] <- if ( tunecontrol$nrepeat == 1 || all(is.na(repeat.errors)) ) repeat.errors[1] else tunecontrol$repeat.aggregate(repeat.errors, na.rm=TRUE)  #default: min
    }# rof sample of CV or bootstrap
    model.errors[para.set] <- tunecontrol$sampling.aggregate(sampling.errors, na.rm=TRUE) #default: mean
    model.variances[para.set] <- tunecontrol$sampling.dispersion(sampling.errors, na.rm=TRUE) #default: sd
  }# rof para.set
  
  best <- which.min(model.errors)
  pars <- if (is.null(ranges))  NULL else lapply(parameters[best, , drop = FALSE], unlist)
  
  #return value
  structure(list(best.parameters = parameters[best, , drop = FALSE], 
                 best.performance = model.errors[best],
                 method = "svm", #if (!is.character(method)) deparse(substitute(method)) else method, 
                 nparcomb = nrow(parameters), train.ind = train.ind,
                 sampling = switch(tunecontrol$sampling, 
                                   fix = "fixed training/validation set", bootstrap = "bootstrapping", 
                                   cross = if (tunecontrol$cross == n) "leave-one-out" else paste(tunecontrol$cross, "-fold cross validation", sep = "")),
                 tune.measure = tune.measure, # mkuhn: added this
                 performances = if (tunecontrol$performances) cbind(parameters, error = model.errors, dispersion = model.variances), 
                 best.model = if (tunecontrol$best.model) {
                   modeltmp <- if (useFormula) do.call(svm, c(list(train.x, 
                                                                      data = data), pars, list(...))) else do.call(svm, 
                                                                                                                   c(list(x = train.x, y = NULL), pars, list(...)))
                   call[[1]] <- as.symbol("best.tune")
                   modeltmp$call <- call
                   modeltmp
                 }), class = "tune")
}





#' Tune function from package e1071 as reference with my comments.
#' 
#' Allows to tune hyper-parameter in the model by cross-validation of bootstrap samples.
#' 
#' I added comments to clarify the inner workings of the tuning.
#' In particular for one-class SVMs it apparently uses on validation data as response variable all TRUE.
#' @param method function (possibly as character)
tune_commented <- function (method, train.x, train.y = NULL, data = list(), validation.x = NULL, 
                    validation.y = NULL, ranges = NULL, predict.func = predict, 
                    tunecontrol = tune.control(), ...) 
{
  call <- match.call()
  
  resp <- function(formula, data) { model.response(model.frame(formula, data)) }
  
  classAgreement <- function(tab) {
    n <- sum(tab)
    if (!is.null(dimnames(tab))) {
      lev <- intersect(colnames(tab), rownames(tab))
      p0 <- sum(diag(tab[lev, lev]))/n
    }
    else {
      m <- min(dim(tab))
      p0 <- sum(diag(tab[1:m, 1:m]))/n
    }
    p0
  }
  
  # mkuhn, 20140701: allow also cross-validation with validation data
  if (tunecontrol$sampling == "cross") {
    validation.x <- validation.y <- NULL
    if (! is.null(validation.x) || ! is.null(validation.y)) cat("Set validation data to NULL. [MK]\n")
  }

  #mkuhn, 20140701: allow matrix for train.x and no response.
  #Original was: useFormula <- is.null(train.y)
  useFormula <- inherits(train.x, "formula") #is.null(train.y)
  if (useFormula && (is.null(data) || length(data) == 0)) data <- model.frame(train.x)
  if (is.vector(train.x)) train.x <- t(t(train.x)) # vector as nx1-matrix
  if (is.data.frame(train.y)) train.y <- as.matrix(train.y)
  if (!is.null(validation.x)) tunecontrol$fix <- 1
  
  # PREPARE SAMPLING (training/test data splits)
  n <- nrow(if (useFormula) data else train.x)
  perm.ind <- sample(n)
  if (tunecontrol$sampling == "cross") {
    if (tunecontrol$cross > n) 
      stop(sQuote("cross"), " must not exceed sampling size!")
    if (tunecontrol$cross == 1) 
      stop(sQuote("cross"), " must be greater than 1!")
  }
  
  # perm.ind is the initial permutation of the training data
  #
  # train.ind is a 
  # + list of a single subset of the permutation [fix], 
  # + list of permutations [boot] 
  # + list of subsets of the permutation [cross]
  train.ind <- if (tunecontrol$sampling == "cross") 
    tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x]) # get k=tunecontrol$cross subsets from the initial permuation perm.ind
  else if (tunecontrol$sampling == "fix") 
    list(perm.ind[1:trunc(n * tunecontrol$fix)])
  else lapply(1:tunecontrol$nboot, function(x) sample(n, n * tunecontrol$boot.size, replace = TRUE))  #bootstrap with own sample()-command with replacement.
  
  #mkuhn, 20140701: debug output
  cat( sprintf("Length of samplings from training data: %s. [MK]\n", paste(sapply(train.ind, length), collapse = ", ")) )
  
  # Start with parameter-grid
  parameters <- if (is.null(ranges)) data.frame(dummyparameter = 0) else expand.grid(ranges)
  p <- nrow(parameters)
  
  if (!is.logical(tunecontrol$random)) {
    if (tunecontrol$random < 1) 
      stop("random must be a strictly positive integer")
    if (tunecontrol$random > p) tunecontrol$random <- p
    parameters <- parameters[sample(1:p, tunecontrol$random), ]
  }
  
  # mkuhn, 20140701: init vectors with correct length
  model.variances <- model.errors <- numeric(p) #c()
  for (para.set in 1:p) {
    sampling.errors <- numeric(length(train.ind)) #c()
    for (sample in 1:length(train.ind)) {
      repeat.errors <- numeric(tunecontrol$nrepeat) #c()
      for (reps in 1:tunecontrol$nrepeat) { # repeats for random algorithms?
        
        pars <- if (is.null(ranges)) NULL else lapply(parameters[para.set, , drop = FALSE], unlist)
        
        
        ### TRAINING data:
        
        #train model on training data
        model <- if (useFormula) 
          do.call(method, c(list(train.x, data = data, 
                                 subset = train.ind[[sample]]), pars, list(...)))
        else do.call(method, c(list(train.x[train.ind[[sample]], ], y = train.y[train.ind[[sample]]]), pars, 
                               list(...)))
   
        # mkuhn, 20140701: debug output of model
        cat(sprintf("Trained SVM-model of type %d on sample of training data, %d rows. [MK]\n", model$type, NROW(train.x[train.ind[[sample]],])))
        cat("number of unique indices in training sample: ",length(unique(train.ind[[sample]])), " [MK]\n")
        #ZZZ is bootstrapping distrubing the one-class SVM, because we have so many duplicates?
        
        ### TEST data
        
        # PREDICTION
        # model prediction on test- or validation-data
        pred <- predict.func(model, if (!is.null(validation.x)) 
          validation.x
          else if (useFormula) 
            data[-train.ind[[sample]], , drop = FALSE]
          else if (inherits(train.x, "matrix.csr")) 
            train.x[-train.ind[[sample]], ]
          else train.x[-train.ind[[sample]], , drop = FALSE])
        
        # mkuhn, 20140701: debug output for pred
        cat(sprintf("Dimension of test data: %s. [MK]\n", paste(dim(train.x[-train.ind[[sample]], , drop = FALSE]), collapse="-")))
        # mkuhn, 20140701: debug output of pred
        print(table(pred))
        
        
        # OBSERVED
        # observed (=true) y in test- or validation-data
        
        # mkuhn, 20140701: if only validation.y is given (and not validation.x) it will be used for evaluation in case of one-class SVM.
        # The vector validation.y should be coded like what predict() does at one-class SVM, IN=TRUE, OUT=FALSE.
        #
        # Originally was:
        #         true.y <- if (!is.null(validation.y)) 
        #           validation.y
        # mkuhn, 20140701: added first if branch
        true.y <- if ( !is.null(validation.y) && is.null(validation.x) && NROW(validation.y) == n ) {
          cat("Use validation.y as secret output of training data for one-class SVM [MK]\n")
          validation.y[-train.ind[[sample]] ]  # abuse validation.y here: if only validation.y is given it serves as a secret response vector for one-class SVM
        } else if (!is.null(validation.y)) 
          validation.y
        else if (useFormula) { # validation.y is NULL in this else-branch here
          if (!is.null(validation.x)) 
            resp(train.x, validation.x)
          else resp(train.x, data[-train.ind[[sample]], ])
        } # else: no validation.y and no formula. Hence, in this branch train.y should be a (1-column?!) matrix
        else train.y[-train.ind[[sample]]]
        
        # mkuhn, 20140701: debug output of true.y
        cat("True Y [MK]\n")
        print(table(true.y))
        
        # mkuhn, 20140701: all TRUE as fall-back (in case of novelty detection!)
        if (is.null(true.y)){
          cat("Setting response for parameter ", paste(pars), " for rep ", reps, " in sample ", sample, " to all TRUE [MK]\n")
          true.y <- rep(TRUE, length(pred))
        }
        
        # ERROR
        # calculate the error in this repetition: with own error-function (if available) or defaults to accuracy [classification] or mean squared error (MSE) [regression]
        repeat.errors[reps] <- if (!is.null(tunecontrol$error.fun)) 
          tunecontrol$error.fun(true.y, pred)
        else if ((is.logical(true.y) || is.factor(true.y)) && 
                   (is.logical(pred) || is.factor(pred) || is.character(pred))) {
          confTab <- table(pred, true.y)
          print(confTab)
          1 - classAgreement(confTab)
        } else if (is.numeric(true.y) && is.numeric(pred)) 
          crossprod(pred - true.y)/length(pred)
        else stop("Dependent variable has wrong type!")
      }#rof reps
      # aggregating the errors in repetion (with same parameters and same training/test sample split)
      sampling.errors[sample] <- tunecontrol$repeat.aggregate(repeat.errors)  #default: min
    }# rof sample
    model.errors[para.set] <- tunecontrol$sampling.aggregate(sampling.errors) #default: mean
    model.variances[para.set] <- tunecontrol$sampling.dispersion(sampling.errors) #default: sd
  }# rof para.set
  
  best <- which.min(model.errors)
  pars <- if (is.null(ranges))  NULL else lapply(parameters[best, , drop = FALSE], unlist)
  
  #return value
  structure(list(best.parameters = parameters[best, , drop = FALSE], 
                 best.performance = model.errors[best], method = if (!is.character(method)) deparse(substitute(method)) else method, 
                 nparcomb = nrow(parameters), train.ind = train.ind,
                 sampling = switch(tunecontrol$sampling, 
                                   fix = "fixed training/validation set", bootstrap = "bootstrapping", 
                                   cross = if (tunecontrol$cross == n) "leave-one-out" else paste(tunecontrol$cross, "-fold cross validation", sep = "")),
                 performances = if (tunecontrol$performances) cbind(parameters, error = model.errors, dispersion = model.variances), 
                 best.model = if (tunecontrol$best.model) {
                   modeltmp <- if (useFormula) do.call(method, c(list(train.x, 
                                                                      data = data), pars, list(...))) else do.call(method, 
                                                                                                                   c(list(x = train.x, y = train.y), pars, list(...)))
                   call[[1]] <- as.symbol("best.tune")
                   modeltmp$call <- call
                   modeltmp
                 }), class = "tune")
}

