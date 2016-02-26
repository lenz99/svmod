# mkuhn, 2015-03-03
# make oneclass and multiclass SVMs available to the mlr-framework as 2-class classifiers



#' ocsvm2: two-class classifier implemented as one-class SVM in mlr
#' 
#' A one-class SVM is used to do a 2-class classification.
#' This classifier is made available to the mlr-framework.
#' @export
makeRLearner.classif.ocsvm2 = function() {
  mlr::makeRLearnerClassif(
    cl = "classif.ocsvm2",
    package = "e1071",
    par.set = ParamHelpers::makeParamSet(
      ParamHelpers::makeNumericLearnerParam(id = "nu", default = 0.5, lower=0), #, requires = quote(type=="nu-classification")),
      ParamHelpers::makeNumericVectorLearnerParam("class.weights", len = NA_integer_, lower = 0),
      ParamHelpers::makeDiscreteLearnerParam(id = "kernel", default = "radial", values = c("linear", "polynomial", "radial", "sigmoid")),
      ParamHelpers::makeIntegerLearnerParam(id = "degree", default = 3L, lower = 1L, requires = quote(kernel=="polynomial")),
      ParamHelpers::makeNumericLearnerParam(id = "coef0", default = 0, requires = quote(kernel=="polynomial" || kernel=="sigmoid")),
      ParamHelpers::makeNumericLearnerParam(id = "gamma", lower = 0, requires = quote(kernel!="linear")),
      ParamHelpers::makeNumericLearnerParam(id = "tolerance", default = 0.001, lower = 0),
      ParamHelpers::makeLogicalLearnerParam(id = "shrinking", default = TRUE),
      ParamHelpers::makeNumericLearnerParam(id = "cachesize", default = 40L)
    ),
    properties = c("twoclass", "numerics", "factors"),  # "prob" not supported for one-class SVM (yet!?)
    name = "oneclass Support Vector Machines (libsvm) as two-class classifier",
    short.name = "ocsvm2",
    note = "The classification runs via a one-class SVM as twoclass classifier"
  )
}


#' Trains a one-class SVM on the given data using all features except the response variable.
#' @return a list object of length 2 that contains the fitted one-class SVM and the status variable
#' @export
trainLearner.classif.ocsvm2 = function(.learner, .task, .subset, .weights = NULL,  ...) {
  
  #f <- getTaskFormula(.task)
  #e1071::svm(f, data = myData, probability = .learner$predict.type == "prob", ...)
  
  myData <- mlr::getTaskData(.task, .subset)
  
  # drop response column
  targetName <- mlr::getTaskTargetNames(.task)
  myStatus <- myData[, targetName]
  myData <- myData[, ! names(myData) %in% targetName]
  
  return( list(ocsvm2=e1071::svm(~ ., data = myData, probability = FALSE, ...), status=myStatus) )
}


#' Prediction method for the one-class SVM classifier
#' 
#' It utilizes the one-class SVM predict method that returns \code{TRUE} for normal cases and \code{FALSE} for outliers.
#' @return Returns a factor with levels IN and OUT.
#' @export
predictLearner.classif.ocsvm2 = function(.learner, .model, .newdata, ...) {
  
  myModelObj <- .model$learner.model
  myOCSVM <- myModelObj[["ocsvm2"]]
  myStatus <- myModelObj[["status"]]
  stopifnot( identical(STATUS_LEVELS, levels(myStatus)) )
  
  if(.learner$predict.type == "response") {
    pred.isNormal <- predict(myOCSVM, newdata = .newdata, ...) #TRUE=1=IN and FALSE=0=OUT
    
    return( factor(STATUS_LEVELS[2L-as.numeric(pred.isNormal)]) )
    
  } else {
    warning("Probabilities are not (yet?!) supported in libsvm's oneclass SVMs")
    #attr(predict(myOCSVM, newdata = .newdata, probability = TRUE, ...), "probabilities")
    return(invisible(NULL))
  }
}



#' mcsvm2: two-class classifier in mlr implemented as multi-class SVM
#' 
#' A multi-class SVM is used to do a 2-class classification.
#' This classifier is made available to the mlr-framework.
#' @export
makeRLearner.classif.mcsvm2 = function() {
  mlr::makeRLearnerClassif(
    cl = "classif.mcsvm2",
    package = "e1071",
    par.set = ParamHelpers::makeParamSet(
      ParamHelpers::makeNumericLearnerParam(id = "nu", default = 0.5, lower=0), #, requires = expression(type=="nu-classification")),
      ParamHelpers::makeNumericVectorLearnerParam("class.weights", len = NA_integer_, lower = 0),
      ParamHelpers::makeDiscreteLearnerParam(id = "kernel", default = "radial", values = c("linear", "polynomial", "radial", "sigmoid")),
      ParamHelpers::makeIntegerLearnerParam(id = "degree", default = 3L, lower = 1L, requires = expression(kernel=="polynomial")),
      ParamHelpers::makeNumericLearnerParam(id = "coef0", default = 0, requires = expression(kernel=="polynomial" || kernel=="sigmoid")),
      ParamHelpers::makeNumericLearnerParam(id = "gamma", lower = 0, requires = expression(kernel!="linear")),
      ParamHelpers::makeNumericLearnerParam(id = "tolerance", default = 0.001, lower = 0),
      ParamHelpers::makeLogicalLearnerParam(id = "shrinking", default = TRUE),
      ParamHelpers::makeNumericLearnerParam(id = "cachesize", default = 40L)
    ),
    properties = c("twoclass", "numerics", "factors"), #"prob" mkuhn, 20150425: leave out for now
    name = "two-class classifier via multi-class Support Vector Machines (libsvm)",
    short.name = "mcsvm2",
    note = "The classification runs via a multi-class SVM as twoclass classifier"
  )
}


#' Trains a 2-class classifier based on a multi-class SVM on the given data using all features
#' @param SVtype the actual multi-class label
#' @return a fitted svm model
#' @export
trainLearner.classif.mcsvm2 = function(.learner, .task, .subset, .weights = NULL, ymultLabel="SVtype", ...) {
  # mlr definition:
  #f <- getTaskFormula(.task)
  #e1071::svm(f, data = myData, probability = .learner$predict.type == "prob", ...)
  
  
  stopifnot( ymultLabel %in% mlr::getTaskFeatureNames(.task) )
  
  taskData <- mlr::getTaskData(.task, .subset)
  yVar2 <- mlr::getTaskTargetNames(.task)
  
  # remove original 2-class response variable
  myTaskData <- dplyr::select_(taskData, ~-yVar2)
  
  e1071::svm(formula(paste(ymultLabel, "~ .")), data = myTaskData, probability = .learner$predict.type == "prob", ...)
  
}


#' Prediction method for the 2-class SVM classifier implemented as multi-class SVM
#' 
#' It utilizes the multi-class SVM predict method.
#' @return Returns a factor with levels IN and OUT.
#' @export
predictLearner.classif.mcsvm2 = function(.learner, .model, .newdata, ...) {
  #ZZZ continue here! how can I do again the predictions on the two-class variable!
  myModelObj <- .model$learner.model
  myOCSVM <- myModelObj[["ocsvm2"]]
  myStatus <- myModelObj[["status"]]
  stopifnot( identical(STATUS_LEVELS, levels(myStatus)) )
  
  if(.learner$predict.type == "response") {
    pred.isNormal <- predict(myOCSVM, newdata = .newdata, ...) #TRUE=1=IN and FALSE=0=OUT
    
    return( factor(STATUS_LEVELS[2L-as.numeric(pred.isNormal)]) )
    
  } else {
    warning("Probabilities are not (yet?!) supported in libsvm's oneclass SVMs")
    #attr(predict(myOCSVM, newdata = .newdata, probability = TRUE, ...), "probabilities")
    return(invisible(NULL))
  }
}



# #' classif.ocsvm1: really as oneclass classifier (does not really work, though.)
# makeRLearner.classif.ocsvm1 = function() {
#   mlr::makeRLearnerClassif(
#     cl = "classif.ocsvm1",
#     package = "e1071",
#     par.set = makeParamSet(
#       #makeDiscreteLearnerParam(id = "type", default = "C-classification", values = c("C-classification", "nu-classification")),
#       #makeNumericLearnerParam(id = "cost",  default = 1, lower = 0, requires = expression(type=="C-classification")),
#       makeNumericLearnerParam(id = "nu", default = 0.5, lower=0), #, requires = expression(type=="nu-classification")),
#       makeNumericVectorLearnerParam("class.weights", len = NA_integer_, lower = 0),
#       makeDiscreteLearnerParam(id = "kernel", default = "radial", values = c("linear", "polynomial", "radial", "sigmoid")),
#       makeIntegerLearnerParam(id = "degree", default = 3L, lower = 1L, requires = expression(kernel=="polynomial")),
#       makeNumericLearnerParam(id = "coef0", default = 0, requires = expression(kernel=="polynomial" || kernel=="sigmoid")),
#       makeNumericLearnerParam(id = "gamma", lower = 0, requires = expression(kernel!="linear")),
#       makeNumericLearnerParam(id = "tolerance", default = 0.001, lower = 0),
#       makeLogicalLearnerParam(id = "shrinking", default = TRUE),
#       makeNumericLearnerParam(id = "cachesize", default = 40L)
#     ),
#     properties = c("oneclass", "numerics", "factors"),  # "prob" not supported for one-class SVM (yet!?)
#     name = "oneclass Support Vector Machines (libsvm)",
#     short.name = "ocsvm1",
#     note = "My test wrapper in mlr for evaluation of one-class SVM directly as real oneclass novelty detection system"
#   )
# }
# 
# 
# trainLearner.classif.ocsvm1 = function(.learner, .task, .subset, .weights = NULL,  ...) {
#   myData <- mlr::getTaskData(.task, .subset)
#   myStatus <- myData$status
#   myData <- myData[, ! names(myData) %in% c("status", mlr::getTaskTargetNames(.task))]
#   
#   return( list(ocsvm1=e1071::svm(~ ., data = myData, probability = .learner$predict.type == "prob", ...),
#               status=myStatus) )
# }
# 
# predictLearner.classif.ocsvm1 = function(.learner, .model, .newdata, ...) {
#   
#   myModelObj <- .model$learner.model
#   myOCSVM <- myModelObj[["ocsvm1"]]
#   myStatus <- myModelObj[["status"]]
#   
#   if(.learner$predict.type == "response") {
#     pred.isNormal <- predict(myOCSVM, newdata = .newdata, ...)
#     
#     is.na(pred.isNormal) <- ! pred.isNormal
#     # factor with only one single level, outliers are NA
#     return( factor(pred.isNormal, labels = "IN") )
#   } else {
#     warning("Probabilities are not (yet?!) supported in oneclass SVMs in libsvm")
#     #attr(predict(myOCSVM, newdata = .newdata, probability = TRUE, ...), "probabilities")
#     return(invisible(NULL))
#   }
# }