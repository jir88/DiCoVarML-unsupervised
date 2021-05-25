

#' Train and Compute Performance for the kth-Fold
#'
#' General wrapper to train and compute performance for an entire fold (Train (model training) / Test (Held-out performance estimate)). This function compute the performance of each
#' model in the ensemble. Performance can be returned for each model and for optional average prob ensemble and stacked ensembles
#'
#' @importFrom foreach %dopar%
#'
#' @param train_x samples by features matrix for model training
#' @param train_y training class labels
#' @param test_x test samples by feature matrix
#' @param test_y test class labels
#' @param test_ids row metadata for samples. Streamlines downstream analysis and further covariate concatenation
#' @param num_folds number of cross-validation folds for training models and hyper-parameter selection
#' @param num_repeats number of repeats for model train cross-validation
#' @param ranger_ntrees number of trees for random forest models (if applicable)
#' @param ranger_mtry usual random forest mtry parameter (if applicable)
#' @param cvMethod cross validation method. Should be a valid caret cross validation method
#' @param train_stackedModel Should a stacked ensemble be trained?
#' @param avgModelPrediction Should an average probability model be applied. Ensembles models together and in general improves performance.
#' @param metaLearner metalearner for stacked ensemble. Usually a neural net or glm (binary). Can be any valid caret model.
#' @param ensembleModels which models to train. Models that make up the ensemble(if applicable). Must be valid caret models.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{performance} \tab data.frame with performance metrics on test data for each model and ensemble(if applicable). Metrics include: train-test AUC, MCC, Accuracy, etc.    \cr
#'    \tab \cr
#'    \code{roc_plotdata} \tab tidy data to plot roc curves for each model \cr
#'    \tab \cr
#'    \code{models} \tab caret model objects for each model trained  \cr
#'    \tab \cr
#'    }
#'
#' @export
#'
processModel =
  function(train_x,
           train_y ,
           test_x ,
           test_y ,
           test_ids,
           num_folds = 10,num_repeats  = 5,ranger_ntrees = 750,ranger_mtry = 2,cvMethod = "repeatedcv",
           train_stackedModel = F,avgModelPrediction = T,metaLearner = "mlpML",
           ensembleModels = c("gbm","xgbTree","ranger","rangerE","pls")

  ) {


    Class = NULL
    ID = NULL
    Model = NULL
    Resample = NULL
    Status = NULL
    jjj = NULL
    model = NULL

    ## classes
    train_y = stringr::str_replace_all(factor((train_y)),pattern = "-","_")
    test_y = stringr::str_replace_all(factor((test_y)),pattern = "-","_")
    train_y = stringr::str_replace_all(factor((train_y)),pattern = "/","_")
    test_y = stringr::str_replace_all(factor((test_y)),pattern = "/","_")
    classes = (unique(train_y))


    ## Train Model(s) ####
    message("Train Models")
    mdls.dcv = trainML_Models(trainLRs =  data.frame(train_x),
                              testLRs = data.frame(test_x),
                              ytrain= train_y,
                              y_test = test_y,
                              cvMethod = cvMethod,
                              mtry_ = ranger_mtry,
                              numFolds = num_folds,
                              numRepeats = num_repeats,
                              testIDs = test_ids,
                              ntrees = ranger_ntrees,
                              models = ensembleModels)

    ## Extract Model Names and Prediction Matrix
    mdlNames = mdls.dcv$performance$method
    pd = mdls.dcv$predictionMatrix
    ## Data storage
    perf = data.frame()
    ensMatrix = data.frame()
    train_ensprobs = data.frame()
    rocPlots = list()
    model_list = list()

    ## Process Performance #####
    message("Process Performance")
    pb <- progress::progress_bar$new(
      format = " Process Performance [:bar] :percent eta: :eta",
      total = length(mdlNames), clear = FALSE, width= 60)
    for(m in mdlNames){
      ##get Model
      mdl = mdls.dcv$models[[m]]
      bestTune = data.frame( mdl$bestTune)
      suppressMessages(suppressWarnings({
        pred.matrix = dplyr::left_join(bestTune,mdl$pred)
      }))

      ## Get Train Scores
      et = pred.matrix[,as.character(classes)];ytl = pred.matrix$obs
      et = selEnergyPermR::fastImputeZeroes(et)

      ## Optimize weights
      optf = function(x){
        et1 =  log(sweep(et,2,STATS = x ,FUN = "/"))
        votes = sapply(1:nrow(et1), function(x) colnames(et1)[which.max(et1[x,])] )
        AC1 = caret::confusionMatrix(factor(votes,levels = classes),factor(ytl,levels = classes))
        xx = as.matrix(tidyr::spread(data.frame(AC1$table),"Reference","Freq")[,-1])
        mltools::mcc(confusionM = xx)
      }
      nrr = 1000
      weights = compositions::clo(table(train_y))
      cc = foreach::foreach(jjj =1:nrr,.combine = rbind)%dopar%{
        set.seed(jjj)
        prs = as.numeric(compositions::rDirichlet.acomp(1,alpha = weights))
        data.frame(t(prs),MCC = optf(prs))
      }
      weights = as.numeric(cc[which.max(cc$MCC),1:length(classes)])
      if(optf(compositions::clo(table(train_y)))>optf(weights)){
        weights = (compositions::clo(table(train_y)))
      }
      train_mcc = optf(weights)

      ## Train weight adjust scores
      et1 = et[,as.character(classes)];ytl = pred.matrix$Status
      et1 = selEnergyPermR::fastImputeZeroes(et1)
      et1 =  log(sweep(et1,2,STATS = weights ,FUN = "/"))
      train_ensprobs = rbind(train_ensprobs,data.frame(Model = m,Resample = pred.matrix$Resample,
                                                       Status = pred.matrix$obs,
                                                       ID = pred.matrix$rowIndex,et1))

      ## Use weights to adjust probabilities
      et = pd %>%
        dplyr::filter(model==m)
      et1 = et[,as.character(classes)];ytl = et$Status
      et1 = selEnergyPermR::fastImputeZeroes(et1)
      et1 =  log(sweep(et1,2,STATS = weights ,FUN = "/"))
      votes = sapply(1:nrow(et1), function(x) colnames(et1)[which.max(et1[x,])] )
      mroc = pROC::multiclass.roc(et$Status,et1)
      auc_ = pROC::auc(mroc)
      AC1 = caret::confusionMatrix(factor(votes,levels = classes),factor(et$Status,levels = classes))
      xx = as.matrix(tidyr::spread(data.frame(AC1$table),"Reference","Freq")[,-1])
      rocPlots[[m]] = rocPlot(mroc)
      ## store prediction
      et1 = data.frame(Model = m,et[,1:2],et1)
      ensMatrix = rbind(ensMatrix,et1)

      ## Get Confusion Matrix Performance
      if(length(classes)>2){
        pp = t(rowMeans( t(AC1$byClass)))
      }else{
        pp = t(AC1$byClass)
      }

      model_list[[m]] = mdls.dcv$models[[m]]

      trainAUC = mdls.dcv$performance$TrainAUC[mdlNames==m]
      pb$tick()
      perf = rbind(perf,data.frame(Model = m,
                                   train_AUC = trainAUC,
                                   train_MCC = train_mcc,
                                   AUC = auc_,
                                   MCC = mltools::mcc(confusionM = xx),
                                   Accuracy = t(AC1$overall)[1],
                                   pp)
      )
    }


    ## Get Avg Ensemble Performance ####
    train_ensprobs1 = train_ensprobs[,-2:-1] %>%
      dplyr::group_by(ID,Status) %>%
      dplyr::summarise_all(.funs = mean)
    ph = train_ensprobs1[,as.character(classes)]
    votes = sapply(1:nrow(train_ensprobs1), function(x) colnames(ph)[which.max(ph[x,])] )
    mroc = pROC::multiclass.roc(train_ensprobs1$Status,data.frame(ph))
    AC1 = caret::confusionMatrix(factor(votes,levels = classes),factor(train_ensprobs1$Status,levels = classes))
    xx = as.matrix(tidyr::spread(data.frame(AC1$table),"Reference","Freq")[,-1])
    train_mcc = mltools::mcc(confusionM = xx)
    trainAUC = pROC::auc(mroc)
    ensMatrix1 = ensMatrix[,-1] %>%
      dplyr::group_by(ID,Status) %>%
      dplyr::summarise_all(.funs = mean)
    probs.ens = ensMatrix1[,as.character(classes)]
    votes = sapply(1:nrow(probs.ens), function(x) colnames(probs.ens)[which.max(probs.ens[x,])] )
    mroc = pROC::multiclass.roc(ensMatrix1$Status,data.frame(probs.ens))
    rocPlots[["avgEnsemble"]] = rocPlot(mroc)
    auc_ = pROC::auc(mroc)
    AC1 = caret::confusionMatrix(factor(votes,levels = classes),factor(ensMatrix1$Status,levels = classes))
    xx = as.matrix(tidyr::spread(data.frame(AC1$table),"Reference","Freq")[,-1])
    ## Get Confusion Matrix Performance
    if(length(classes)>2){
      pp = t(rowMeans( t(AC1$byClass)))
    }else{
      pp = t(AC1$byClass)
    }
    perf = rbind(perf,data.frame(Model = "avg_ensemble",
                                 train_AUC = trainAUC,
                                 train_MCC = train_mcc,
                                 AUC = auc_,
                                 MCC = mltools::mcc(confusionM = xx),
                                 Accuracy = t(AC1$overall)[1],
                                 pp
    )
    )


    ## Train Stacked Model ####
    if(train_stackedModel){
      message("Training Stacked Model")
      ## Training data stacked model
      stacked.df = tidyr::gather(train_ensprobs,key = "Class",value = "Score",as.character(classes))
      stacked.df$Resample = stringr::str_split(stacked.df$Resample,pattern = "\\.",simplify = T)[,2]
      ## avg Model PRedictions Ensemble
      if(avgModelPrediction){
        stacked.df = stacked.df %>%
          dplyr::select(-Resample) %>%
          dplyr::group_by(ID,Status,Model,Class) %>%
          dplyr::summarise_all(.funs = mean)
      }
      stacked.df = tidyr::unite(stacked.df,"model_class",c(Model,Class))
      mod_class = unique(stacked.df$model_class)
      stacked.df = tidyr::spread(stacked.df,"model_class","Score")


      et = data.frame(stacked.df[,mod_class])
      stacked.model = trainML_Models(trainLRs =  et,testLRs = et,
                                     ytrain = stacked.df$Status,
                                     cvMethod = cvMethod,
                                     mtry_ = 1,numFolds = num_folds,
                                     numRepeats = num_repeats,
                                     y_test = stacked.df$Status,testIDs = NULL,
                                     models = metaLearner)
      ## Compute train MCC
      ph = stacked.model$predictionMatrix[,as.character(classes)]
      votes = sapply(1:nrow(stacked.model$predictionMatrix), function(x) colnames(ph)[which.max(ph[x,])] )
      AC1 = caret::confusionMatrix(factor(votes,levels = classes),factor(stacked.df$Status,levels = classes))
      xx = as.matrix(tidyr::spread(data.frame(AC1$table),"Reference","Freq")[,-1])
      train_mcc = mltools::mcc(confusionM = xx)

      message("Make Predictions")
      ## Get Real PRedictions
      stacked.df = tidyr::gather(ensMatrix,key = "Class",value = "Score",as.character(classes))
      stacked.df = tidyr::unite(stacked.df,"model_class",c(Model,Class))
      stacked.df = tidyr::spread(stacked.df,"model_class","Score")
      ## make predicion
      preds = caret::predict.train(stacked.model$models[[1]],stacked.df[,mod_class],type = "prob")
      votes = caret::predict.train(stacked.model$models[[1]],stacked.df[,mod_class])
      mroc = pROC::multiclass.roc(stacked.df$Status,data.frame(preds))
      rocPlots[["stackedModel"]] = rocPlot(mroc)
      auc_ = pROC::auc(mroc)
      AC1 = caret::confusionMatrix(factor(votes,levels = classes),factor(stacked.df$Status,levels = classes))
      xx = as.matrix(tidyr::spread(data.frame(AC1$table),"Reference","Freq")[,-1])
      ## Get Confusion Matrix Performance
      if(length(classes)>2){
        pp = t(rowMeans( t(AC1$byClass)))
      }else{
        pp = t(AC1$byClass)
      }
      perf = rbind(perf,data.frame(Model = "stacked_ensemble",
                                   train_AUC = stacked.model$performance$TrainAUC,
                                   train_MCC = train_mcc,
                                   AUC = auc_,
                                   MCC = mltools::mcc(confusionM = xx),
                                   Accuracy = t(AC1$overall)[1],
                                   pp
      )
      )
    }


    return(list(performance = perf,roc_plotdata = rocPlots,models = model_list))
  }

