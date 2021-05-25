#' Caret Wrapper for Training Model Ensembles
#'
#' Convenient wrapper for training multiple caret models and making predictions on test data.
#'
#' @param trainLRs samples by features training data. *without labels*
#' @param testLRs  samples by features test data. *without class labels*
#' @param ytrain  training class labels
#' @param y_test test class labels
#' @param testIDs row meta-data for test samples
#' @param models model to be trained. Should be valid caret models names.
#' @param Index_ NULL; User supplied vector of folds for cross validation partitions.
#' @param numRepeats Number of repeats for training cross-validation
#' @param numFolds Number of folds for model training
#' @param cvMethod cross-validation method; should be a valid caret CV method
#' @param ntrees number of trees for random forest models
#' @param mtry_ mtry for random forest models
#' @param num_comp number of components for pls models
#' @param ranger_imp variable importance for random forest models
#' @param bagModels should bagging be used to train the model ensemble
#' @param sampleSize sample size for bagging samples
#' @param seed random seed
#'
#' @return A list containing:\tabular{ll}{
#'    \code{performance} \tab training performance for each model  \cr
#'    \tab \cr
#'    \code{predictionMatrix} \tab prediction probs for test samples for each model  \cr
#'    \tab \cr
#'    \code{models} \tab caret model object for each model  \cr
#'    \tab \cr
#' }
#'
#' @seealso \code{\link[caret]{trainControl}}
#' @export
#'
trainML_Models <-
  function(trainLRs,testLRs,
           ytrain,y_test,
           testIDs=NULL,
           models = c("ranger","svmRadial","gbm","pls"),
           Index_ = NULL,
           numRepeats = 3,
           numFolds = 5,
           cvMethod = "repeatedcv",ntrees = 500,mtry_ = 1,num_comp = 2,ranger_imp = "none",
           bagModels = F,sampleSize,
           seed = 08272008){

    finModels = list()
    train_control <- caret::trainControl(method=cvMethod,
                                         repeats = numRepeats,index = Index_,
                                         number=numFolds,seeds = NULL,
                                         classProbs = TRUE,
                                         savePredictions = T,
                                         allowParallel = TRUE,
                                         summaryFunction = caret::multiClassSummary
    )



    performance_ = data.frame()
    i  = 1
    nsamps = nrow(testLRs)
    prediction_matrix_ =data.frame()
    set.seed(seed)
    models_ = paste(models,1:length(models),sep = "")

    for(mdl in models){

      if(bagModels==T){
        bm = data.frame(Status = ytrain,trainLRs)
        bm_i = sample(1:nrow(trainLRs),size = sampleSize,replace = F)
        trainLRs1 = bm[bm_i,-1]
        rownames(trainLRs1) = paste("ID_",bm_i)
        ytrain1 = bm[bm_i,1]
      }else{
        ytrain1 = ytrain
        trainLRs1 = trainLRs
        rownames(trainLRs1) = paste0("ID_",1:nrow(trainLRs))
      }

      suppressMessages(suppressWarnings({


        if(mdl%in%c("ranger")){
          glm.mdl1 = caret::train(x = trainLRs1 ,
                                  y = ytrain1,
                                  metric = "ROC",
                                  max.depth = 0,
                                  method = "ranger",num.trees=ntrees,
                                  importance = ranger_imp,
                                  tuneGrid = expand.grid(min.node.size=1,splitrule = c("gini"), mtry = c(2,round(sqrt(ncol(trainLRs1))))),
                                  trControl = train_control
          )
        }else if(mdl%in%c("xgbTree")){
          glm.mdl1=caret::train(x = trainLRs1 ,
                                y = ytrain1,
                                metric = "ROC",
                                method = "xgbTree",tuneGrid = expand.grid(nrounds = c(50,150),
                                                                          max_depth = c(6),
                                                                          eta = c(0.1),
                                                                          gamma = c(0),
                                                                          colsample_bytree = c(.1,.2),
                                                                          min_child_weight = 1,
                                                                          subsample = c(.8)),
                                trControl = train_control)

        } else if(mdl%in%c("rangerE")){
          glm.mdl1=caret::train(x = trainLRs1 ,
                                y = ytrain1,
                                metric = "ROC",
                                max.depth = 0,
                                method = "ranger",num.trees=ntrees,
                                importance = ranger_imp,
                                tuneGrid = expand.grid(min.node.size=1,
                                                       splitrule = c("extratrees"),
                                                       mtry = c(1,2,round(sqrt(ncol(trainLRs1))))),
                                trControl = train_control
          )
        } else if(mdl%in%c("gbm")){
          glm.mdl1 = caret::train(x = trainLRs1 ,
                                  y = ytrain1,
                                  metric = "ROC",
                                  method = "gbm",
                                  verbose = F,
                                  trControl = train_control
          )}else if(mdl%in%c("knn")){
            glm.mdl1 = caret::train(x = trainLRs1 ,
                                    y = ytrain1,
                                    metric = "ROC",
                                    method = "knn",
                                    tuneGrid = expand.grid(k = round(sqrt(nrow(trainLRs1)))),
                                    trControl = train_control
            )}else if(mdl%in%c("pls")){
              glm.mdl1 = caret::train(x = trainLRs1 ,
                                      y = ytrain1,
                                      metric = "ROC",
                                      method = "pls",
                                      tuneGrid = expand.grid(ncomp = num_comp),
                                      trControl = train_control
              )}else if(mdl%in%c("rf")){
                glm.mdl1 = caret::train(x = trainLRs1 ,
                                        y = ytrain1,
                                        metric = "ROC",
                                        method = "rf",
                                        tuneGrid = expand.grid(mtry = round(sqrt(ncol(trainLRs1)))),
                                        trControl = train_control
                )}else {
                  glm.mdl1 = caret::train(x = trainLRs1 ,
                                          y = ytrain1,
                                          method=mdl,
                                          metric = "ROC",
                                          trControl = train_control
                  )
                }

        preds = caret::predict.train(glm.mdl1,testLRs, type= "prob")
        preds$model = models_[i]
        if(!is.null(testIDs)){
          preds = cbind(testIDs,preds)
        }
        prediction_matrix_ = rbind(prediction_matrix_,preds)
        ph1 = cbind(data.frame(( data.frame( caret::getTrainPerf(glm.mdl1) ) )))
        ph1$method = models_[i]
        #colnames(ph1) = paste(models[i],colnames(ph1),sep = "_")
        performance_ = rbind(performance_,ph1)

      }))

      message(models_[i])
      finModels[[models_[i]]] = glm.mdl1


      i = i+1

    }



    return(list(performance = performance_,
                predictionMatrix = prediction_matrix_,
                models = finModels)
    )

  }
