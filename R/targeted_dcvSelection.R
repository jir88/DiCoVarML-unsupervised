#' Targeted DCV Feature Selection
#'
#' Performs DCV feature selection to build a model with a specific number of parts.
#'
#' @param trainx train data partition samples x parts/taxa/etc.
#' @param testx test data partition samples x parts/taxa/etc.
#' @param y_label train labels
#' @param y_test test labels
#' @param dcv a dcv matrix can optionally be provided
#' @param lrs.train An optional supplied logratio matrix for training data; must be supplied if dcv scores are pre computed
#' @param lrs.test An optional supplied logratio matrix for test data; must be supplied if dcv scores are pre computed
#' @param minConnected should the ratio network be minimally connected or threshold with dcv_score>0
#' @param ensemble ensemble model definition, or FALSE to skip fitting ensemble models
#' @param tarFeatures number of parts/taxa/etc. to retain (independent from number of ratios to retain)
#' @param imp_factor factor to multiplicative impute count zeros
#' @param select_randomFeatures should random parts/taxa/etc. be selected; useful for benchmarking
#' @param use_rfe should a random forest recursive elimination model be trained
#' @param alpha_ glmnet alpha parameter for fitting the final ridge regression model
#' @param ts.id test set id matrix. Can be found in output list from partitioning functions i.e.  kfoldDataPartition(),lodo_partition(),etc. Can be manually entered nrows = nrows(testx)
#' @param seed Random Seed control for reproducibility
#' @param max_sparsity max sparsity of parts/taxa/etc. For example max_sparsity=0.10 would mean to only retain parts/taxa/etc present in at least 10\% of all samples.
#' @param useRidgeWeights should ensemble model first be weighted by ridge regression coefficients?
#' @param scaledata should train data be scaled (test data scaled based on train)
#'
#' @return a list containing the models and diagnostic outputs
#'
#' @importFrom rlang .data
#' @export
#'
#'
targeted_dcvSelection = function(trainx,
                                 testx,
                                 y_label,
                                 y_test,
                                 dcv=NULL,lrs.train=NULL,lrs.test = NULL,
                                 minConnected = F,
                                 ensemble  = c("ranger","pls","svmRadial","glmnet","rangerE"),
                                 tarFeatures = 5,
                                 imp_factor = 1e-7,
                                 select_randomFeatures = F,use_rfe = T,alpha_ = 0,
                                 ts.id, seed = 08272008, max_sparsity = .9, useRidgeWeights = T,
                                 scaledata = T){

  result = data.frame()
  # compTime = 0

  prob_matrices = list()

  classes = as.character(unique(y_label))

  # Select Ratios -----

  baseDims = ncol(trainx) # number of available features


  # if we're actually trying to select good features
  if(!select_randomFeatures){

      ## Compute DCV
      message("Perform DCV Log Ratio Selection - (Step 1 of 4)")

      # if no pre-computed DiCoVar matrix was supplied, calculate one
      if(is.null(dcv)){
        # compTime = system.time({

          ## impute zeroes for both training and test data
          trainx = data.frame(selEnergyPermR::fastImputeZeroes(trainx,impFactor = imp_factor))
          testx = data.frame(selEnergyPermR::fastImputeZeroes(testx,impFactor = imp_factor))

          ## compute log ratios
          lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = y_label,trainx))
          lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = y_test,testx))

          # compute scores on full training set without CV averaging
          cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                              includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                              rankOrder = F)

        # })
        # pull out the scoring matrix
        dcv = cc.dcv$lrs
      }

      ## Compute Node Strength
      dcv_strength = DiCoVarML::computeDCVStrength(list(dcv  =dcv))


      ## get subcomposition - only features in the DCV log-ratio network
      train_subcomp = subset(trainx,select = dcv_strength$Node[1:tarFeatures])
      test_subcomp = subset(testx,select = dcv_strength$Node[1:tarFeatures])


      ## get pairwise logratios from subcomposition
      lrs.train = selEnergyPermR::calcLogRatio(df = data.frame(Status = y_label,
                                                               train_subcomp))[,-1]
      lrs.test = selEnergyPermR::calcLogRatio(df = data.frame(Status = y_test,
                                                              test_subcomp))[,-1]

      ## recompute dcvScores on the subcomposition
      if(minConnected){ # minimally-connected DCV network
        # calculate DCV scores and rank them, also computing feature counts
        cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = data.frame(Status = y_label,lrs.train),
                                            includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                            rankOrder = T)
        dcv = cc.dcv$lrs
        # keep only a pre-specified number of distinct features
        w = min(which(dcv$nDistinct==tarFeatures))
      }else{ # thresholded DCV network with all positive scores
        # calculate DCV scores, don't bother to rank
        cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = data.frame(Status = y_label,lrs.train),
                                            includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                            rankOrder = F)
        dcv = cc.dcv$lrs
        w = max(which(dcv$rowmean>0))

      }


      ### Select key ratios
      glm.train = subset(lrs.train,select = dcv$Ratio[1:w])
      glm.test = subset(lrs.test,select = dcv$Ratio[1:w])
      dcv = dcv[1:w,]


      ## Verify ratio subset contains all target parts
      nn = stringr::str_split(colnames(glm.train),pattern = "___",simplify = T)
      # how many distinct features are actually present?
      empTarFeats = dplyr::n_distinct(c(nn[,1],nn[,2]))

      # if we got too few features with threshold method, use minimal network instead
      if(empTarFeats<tarFeatures){
          # calculate DCV scores and rank them
          cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = data.frame(Status = y_label,lrs.train),
                                              includeInfoGain = T, nfolds = 1, numRepeats = 1,
                                              rankOrder = T)
          dcv = cc.dcv$lrs
          # grab the correct number of features
          w = min(which(dcv$nDistinct==tarFeatures))
          glm.train = subset(lrs.train,select = dcv$Ratio[1:w])
          glm.test = subset(lrs.test,select = dcv$Ratio[1:w])
          dcv = dcv[1:w,]

      }


    }else{ # we're selecting random features, probably for null model benchmarks
      cn = colnames(trainx)
      # compTime = 0
      rand_feats = sample(cn,replace = F)[1:tarFeatures]
      train_subcomp = subset(trainx,select = rand_feats)
      test_subcomp = subset(testx,select = rand_feats)
      ## get pair logratio from subcomp
      glm.train = selEnergyPermR::calcLogRatio(df = data.frame(Status = y_label,train_subcomp))[,-1]
      glm.test = selEnergyPermR::calcLogRatio(df = data.frame(Status = y_test,test_subcomp))[,-1]
    }




    ## Compute Number of Parts (features)
    cn = colnames(glm.train)
    uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
    n_parts  = dplyr::n_distinct(uniqueParts)


    ## scale data to unit standard deviation
    if(scaledata ){
      # calculate scaling factors on training set
      pp = caret::preProcess(glm.train,method = "scale")
      glm.train <- stats::predict(pp, glm.train)
      # use same scaling factors on test set
      glm.test     <- stats::predict(pp, glm.test)
    }


    # Train Ridge Model ----
    message("Train Ridge Regression Model - (Step 2 of 4)")

    # 2 groups, or >2 groups?
    type_family = dplyr::if_else(length(classes)>2,"multinomial","binomial")
    # 10-fold CV of ridge regression on standardized training data
    # compTime2 = system.time({
      cv.clrlasso <- glmnet::cv.glmnet(as.matrix(glm.train),y_label, standardize=F, alpha=alpha_,family=type_family)
    # })
    if(type_family=="binomial"){
      # pull optimal model coefficients
      features = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))
      # drop the intercept (why???)
      features = features[-1,]
      # drop any features with coefficient == 0
      features = features[abs(features)>0]
      # report number of features we ended up with
      length(features)
      # same thing, but keep zero coefficients
      c = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))[-1,]
    }else{
      # multinomial coefficient extraction
      features = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))
      feat.df = data.frame()
      for(o in 1:length(features)){
        ph = as.matrix(features[[o]])
        feat = ph[-1,]
        keep = feat[abs(feat)>0]
        feat.df = rbind(feat.df,data.frame(Ratio = names(keep), coef = as.numeric(keep)))
      }
      feat.df = feat.df %>%
        dplyr::group_by(.data$Ratio) %>%
        dplyr::summarise(coef = sum(.data$coef)) %>%
        dplyr::filter(.data$coef != 0)

    }

    ## make predictions using the ridge regression model

    # get fitted probabilities for the test set
    p = stats::predict(cv.clrlasso, newx = as.matrix(glm.test), s = "lambda.min",type = "response")
    # record model and (scaled) data used to train and test it
    model_ = list(mdl = cv.clrlasso,data = list(train = glm.train,test = glm.test))
    if(type_family=="binomial"){
      # construct ROC curve with ridge model probabilities
      mroc = pROC::roc(y_test,p[, 1])
      # calculate and report AUC for ridge model
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
      # record class probabilities from ridge model
      pmat = data.frame(p,1-p);colnames(pmat) = classes

    }else{
      ## multiclass ROC and AUC calculations
      mroc = pROC::multiclass.roc(y_test,p[,,1])
      mroc.dcvlasso = pROC::auc(mroc);mroc.dcvlasso
      pmat = data.frame(p[,,1])#;colnames(pmat) = classes

    }
    ## Save Performance of ridge regression model
    perf = data.frame(Approach = "DCV-ridgeRegression",
                      # AUC of ridge model on the test split
                      AUC = as.numeric(mroc.dcvlasso),
                      # number of features in model
                      number_parts = n_parts,
                      # number of logratios in model
                      number_ratios = ncol(glm.train),
                      comp_time = 1,#compTime[3]+compTime2[3],
                      # number of features before feature selection
                      base_dims = baseDims)
    result = rbind(result,perf)
    # save predicted probabilities
    prob_matrices[["ridgeRegression"]] = data.frame(Status = y_test,pmat)




    # Train Ensemble with Ridge Features ----

    # optionally scale features using ridge regression weights before
    # fitting an ensemble model

    # use ridge weights
    if(type_family=="binomial"){ # if 2 groups
      # grab features used in the ridge regression model
      train_data2 = subset(glm.train,select = unique(names(features)))
      test_data2 = subset(glm.test,select = unique(names(features)))

      # scale features by ridge weights/coefficients
      if(useRidgeWeights){
        train_data2 = sweep(train_data2,MARGIN = 2,STATS = c,FUN = "*")
        test_data2 = sweep(test_data2,MARGIN = 2,STATS = c,FUN = "*")
      }

    }else{ # if >2 groups
      # grab features used in the ridge regression model
      train_data2 = subset(glm.train,select = feat.df$Ratio)
      test_data2 = subset(glm.test,select = feat.df$Ratio)

      # scale features by ridge weights/coefficients
      if(useRidgeWeights){
        train_data2 = sweep(train_data2,MARGIN = 2,STATS = feat.df$coef,FUN = "*")
        test_data2 = sweep(test_data2,MARGIN = 2,STATS = feat.df$coef,FUN = "*")
      }

    }

    # copy the (possibly ridge weighted) data
    weight.train = train_data2
    weight.test = test_data2

    ## Train Ensemble Model

    # if user has supplied ensemble models to use
    if(is.character(ensemble) & sum(is.na(ensemble)) == 0) {
      message("Train Ensemble Model - (Step 3 of 4)")
      ph = trainML_Models(trainLRs = train_data2,
                          testLRs = test_data2,
                          ytrain = y_label,
                          y_test = y_test,
                          testIDs = ts.id,
                          models = ensemble )
      # get test set predictions for each model
      pmat = ph$predictionMatrix
      # average predictions across all models
      pmat = pmat %>%
        dplyr::group_by(.data$ID, .data$Status) %>%
        dplyr::select(-.data$model) %>%
        dplyr::summarise_all(.funs = mean)
      pmat = data.frame(pmat)
      # calculate ROC on averaged predictions
      mroc = pROC::multiclass.roc(pmat$Status,pmat[,classes]);mroc
      ## Compute Number of Parts (again)
      cn = colnames(train_data2)
      n_ratios = length(cn)
      uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
      n_parts  = dplyr::n_distinct(uniqueParts)
      ## Save Performance
      perf = data.frame(Approach = "DCV-ridgeEnsemble",
                        AUC = as.numeric(pROC::auc(mroc)),
                        number_parts = n_parts,
                        number_ratios = ncol(train_data2),
                        comp_time = 1, #compTime[3]+compTime2[3],
                        base_dims = baseDims)
      # add performance to existing ridge model performance data frame
      result = rbind(result,perf)
      # save averaged predictions from ensemble
      prob_matrices[["ridgeEnsemble"]] = pmat
    } else {
      message("Skipping Ensemble Model - (Step 3 of 4)")
    }


    message("Train RFE Model - (Step 4 of 4)")
    # Train RFE model ----
    rfe_Features = NULL
    if(use_rfe){ # if user wants a random forest elimination analysis
      suppressMessages(suppressWarnings({

        type = "Dense"
        ## Feature Selection: rf-RFE

        # start by removing highly correlated logratios
        # correlation threshold
        tc = .99
        c1.cor = stats::cor(glm.train,method = "spearman")
        # find highly-correlated logratios
        c.fc = data.frame(Ratio = caret::findCorrelation(c1.cor,cutoff = tc,names = T))
        # drop these ratios
        keep =!colnames(glm.train) %in% c.fc$Ratio
        glm.train = glm.train[,keep]
        glm.test = glm.test[,keep]

        # do the RFE analysis (similar to rfcv function in randomForest package)
        # compTime2 = system.time({
          pp = rfeSelection.ByMetric(train_ratio = glm.train,
                                     test_ratio = glm.test,
                                     ytrain =y_label,
                                     ntrees = 750,
                                     sets = 10,
                                     impMeasure = "impurity_corrected",
                                     kfold = 5,
                                     minPercentFeatReturn = .3)
        # })

        # get the reduced set of ratios
        train_data2 = pp$reducedTrainRatios
        test_data2 = pp$reducedTestRatio
        message("number of features = ",ncol(train_data2))

        ## Train Ensemble Model
        # if user has supplied ensemble models to use
        if(is.character(ensemble) & sum(is.na(ensemble)) == 0) {
          ph = trainML_Models(trainLRs = train_data2,
                              testLRs = test_data2,
                              ytrain = y_label,
                              y_test = y_test,
                              testIDs = ts.id,
                              models = ensemble)

          ## Compute Performance

          pmat = ph$predictionMatrix
          pmat = pmat %>%
            dplyr::group_by(.data$ID, .data$Status) %>%
            dplyr::select(-.data$model) %>%
            dplyr::summarise_all(.funs = mean)
          pmat = data.frame(pmat)
          mroc = pROC::multiclass.roc(pmat$Status,pmat[,classes]);mroc


          ## Compute Number of Part
          cn = colnames(train_data2)
          n_ratios = length(cn)
          uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
          n_parts  = dplyr::n_distinct(uniqueParts)

          ## Save Performance
          perf = data.frame(Approach = "DCV-rfRFE",AUC = as.numeric(pROC::auc(mroc)),
                            number_parts = n_parts,number_ratios = ncol(train_data2),
                            comp_time = 1, #compTime[3]+compTime2[3],
                            base_dims = baseDims)
          result = rbind(result,perf)
          prob_matrices[["rfRFE"]] = pmat
        }
      }))
      rfe_Features = list(train = train_data2,test =test_data2)
    }




    # Return results ----
    if(is.character(ensemble) & sum(is.na(ensemble)) == 0) {
      all_model_preds <- ph$predictionMatrix
      ensemble_pmat <- pmat
    } else {
      all_model_preds <- NULL
      ensemble_pmat <- NULL
    }
    return(list(Performance = result, # summary stats for fitted models
                # test sample predictions from individual ensemble components
                all_model_preds = all_model_preds,
                # the ridge regression model and its input data
                glm_model = model_,
                # optionally ridge-weighted features used for ensemble fitting
                weighted_features = list(train = weight.train, test = weight.test),
                # parts in the DCV matrix, on the simplex, zero imputed
                part_matrices = list(train = train_subcomp, test = test_subcomp),
                # coefficients from glm_model, without intercept
                ridge_coefficients = c,
                # class probabilities from each fitted model type
                probMatrices  = prob_matrices,
                # DCV row means for final logratios
                final_dcv = dcv,
                # prediction probabilities for test samples for each ensemble component
                ensemble_pmat = pmat,
                # features selected by RFE
                rfe_features = rfe_Features,
                # fitted test set probabilities from the ridge model
                ridge_pmat = p ))




}



