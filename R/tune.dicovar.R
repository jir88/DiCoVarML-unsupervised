#' Tune Differential Compositional Variation Machine Learning (DiCoVarML) Models
#'
#' A function that takes raw data, transforms it into the logratio space, and
#' tests DiCoVarML model performance using subsets of the logratios. Model fitting
#' will use the foreach package to take advantage of multicore processing, if
#' a cluster has been registered.
#'
#' @param X matrix-like object with one row per sample and one column per part
#' @param Y vector containing response variable for each sample
#' @param sample_ids character vector of unique labels or names for each sample
#' @param dv_labels numeric vector identifying each sample as part of the discovery dataset (1) or validation dataset (2)
#' @param scale_data Should data be scaled to unit standard deviation?
#' @param performRFE Should random forest feature elimination be used when fitting DCV models?
#' @param useRidgeWeight Should features be scaled using ridge regression weights before DCV ensemble fitting?
#' @param min_connected Should logratio network be constructed with the minimum number of edges required to include a given number of parts, or should all edges with DCV scores >0 be retained?
#' @param ensemble character vector of models to use in the DCV ensemble
#' @param max_sparsity Parts with a higher proportion of zero/missing values will be dropped before zero imputation.
#' @param test.parts integer vector giving the number of parts to be tested
#' @param k_fold estimate model performance using k-fold cross validation
#' @param repeats the number of times to repeat k-fold cross validation
#' @param seed random seed to use for reproducible results, or NULL to use the current seed
#'
#' @return A list containing:\tabular{ll}{
#'    \code{cv_folds} \tab The cross-validation folds used for performance estimation. \cr
#'    \tab \cr
#'    \code{train_auc} \tab The training set ROC AUC values for each model type and target part count. \cr
#'    \code{inner_perf} \tab ROC AUCs from the inner cross-validation. \cr
#'    \code{target_max} \tab the number of parts which yields the highest AUC \cr
#'    \code{target_1se} \tab the smallest number of parts where AUC is within 1 standard error of optimum \cr
#'    \tab \cr
#'}
#' @export
#'
#' @seealso \code{\link[diffCompVarRcpp]{dcvScores}}
#'
#' @importFrom foreach %:%
#'
tune.dicovar <- function(X,
                         Y,
                         sample_ids,
                         dv_labels,
                         scale_data = TRUE,
                         performRFE = FALSE,
                         useRidgeWeight = FALSE,
                         min_connected = FALSE,
                         ensemble = c("ranger","pls","svmRadial","glmnet","rangerE"),
                         max_sparsity = .9,
                         test.parts = c(4, 8, 12),
                         k_fold = 2,
                         repeats = 5,
                         seed = NULL
                         ) {
  # get levels from response variable
  classes = as.character(levels(as.factor(Y)))

  # set up ROC AUC function depending on number of classes
  if(length(classes) == 2) {
    frm <- stats::as.formula(paste0("Status ~ ", paste(classes[-1], collapse = "+")))
    # df - all data
    # mdl - model for which we're calculating ROC AUC
    roc_auc <- function(df, mdl) {
      df_ <- dplyr::filter(df, model == mdl)
      return(as.numeric(pROC::auc(formula = frm, data = df_, levels = classes,
                                  direction = "<")))
    }
  } else if(length(classes) > 2) {
    frm <- stats::as.formula(paste0("Status ~ ", paste(classes, collapse = "+")))
    # df - all data
    # mdl - model for which we're calculating ROC AUC
    roc_auc <- function(df, mdl) {
      df_ <- dplyr::filter(df, model == mdl)

      return(as.numeric(pROC::auc(pROC::multiclass.roc(formula = frm, data = df_,
                                                       levels = classes,
                                                       direction = ">"))))
    }
  } else {
    simpleError("Data needs at least 2 classes, but only found 1!")
  }

  # get part of the current random seed if none has been specified
  if(!is.numeric(seed) || is.na(seed)) {
    seed <- .Random.seed[3]
  }

  dcv_samples <- data.frame(Status = Y, X)
  rownames(dcv_samples) <- sample_ids

  # 'leave-one-dataset-out' partition to split out validation samples
  allData <- lodo_partition(dcv_samples, dataset_labels = dv_labels, seed = seed)

  # pre-process the data
  ttData <- extractTrainTestSplit(foldDataList = allData,
                                  fold = 2,## leave out true test data which is fold 2
                                  permLabels = FALSE,
                                  maxSparisty = max_sparsity,
                                  extractTelAbunance = FALSE)
  ##get discovery/validation partitions
  train.data = ttData$train_Data
  test.data = ttData$test_data
  ## Compute Total Parts
  number_parts = ncol(train.data);
  nsamps = nrow(train.data)
  # table(ttData$y_test)

  # DiCoVar ----------

  ## Cross-Validation ----

  # First, we do 2-fold cross validation, repeated 5 times, to pick the optimum
  # number of features to include in the final model.

  # create a list of CV folds
  # set RNG seed
  set.seed(seed)
  # RNG seeds for each repeat
  rep_seeds <- stats::runif(n = repeats, min = -.Machine$integer.max,
                     max = .Machine$integer.max)
  cv_folds <- foreach::foreach(sd1 = 1:repeats, rsd = rep_seeds, .packages = c("DiCoVarML", "selEnergyPermR", "caret")) %dopar% {
    set.seed(rsd)
    overll_folds = caret::createFolds(ttData$y_train, k = k_fold, list = F)
    innerfold_data = lodo_partition(data.frame(Status = ttData$y_train,
                                               ttData$train_Data),
                                    dataset_labels = overll_folds,
                                    rsd)
    list(sd1 = sd1, rng_seed = rsd, k_fold = k_fold, innerfold_data = innerfold_data)
  }

  # calculate logratios and DCV scores for each set of CV folds
  lr_dcv_folds <- foreach::foreach(fld = cv_folds, .combine = c) %:%
    foreach::foreach(f = 1:fld$k_fold, .packages = c("DiCoVarML", "selEnergyPermR", "foreach"), .combine = c) %dopar% {
      message(paste("Starting rep", fld$sd1, "fold", f))
      ## Partition inner fold
      innerFold = DiCoVarML::extractTrainTestSplit(foldDataList = fld$innerfold_data,
                                                   fold = f,
                                                   maxSparisty = max_sparsity,
                                                   extractTelAbunance = F)
      ## impute zeroes
      message(paste0("R", fld$sd1, "F", f, " - impute zeroes"))
      trainx = data.frame(selEnergyPermR::fastImputeZeroes(innerFold$train_Data,
                                           impFactor = innerFold$imp_factor))
      testx = data.frame(selEnergyPermR::fastImputeZeroes(innerFold$test_data,
                                          impFactor = innerFold$imp_factor))

      ## compute log ratios
      message(paste0("R", fld$sd1, "F", f, " - compute logratios"))
      lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = innerFold$y_train,
                                                          trainx))
      lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = innerFold$y_test,
                                                         testx))
      # calculate DCV scores
      message(paste0("R", fld$sd1, "F", f, " - calculate DCV scores"))
      cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.train,
                                          includeInfoGain = T, nfolds = 1,
                                          numRepeats = 1, rankOrder = F,
                                          seed_ = fld$rng_seed)

      # estimate performance with different numbers of features
      tmp <- foreach::foreach(tar_Features = test.parts, .packages = c("DiCoVarML", "magrittr")) %do% {
        tar_dcvInner = targeted_dcvSelection(trainx = trainx,
                                             minConnected = min_connected,
                                             useRidgeWeights = useRidgeWeight,
                                             use_rfe = performRFE,
                                             scaledata = scale_data,
                                             testx = testx,
                                             dcv = cc.dcv$lrs,
                                             lrs.train = lrs.train,
                                             lrs.test = lrs.test,
                                             y_label = innerFold$y_train,
                                             seed = fld$rng_seed,
                                             ensemble = ensemble,
                                             y_test = innerFold$y_test,
                                             tarFeatures = tar_Features,
                                             ts.id = innerFold$test_ids,
                                             max_sparsity = max_sparsity
        )
        # extract model performance
        perf = data.frame(Seed = fld$sd1,Fold = f,tar_Features ,tar_dcvInner$Performance)
        # calculate model AUCs on test data
        pmat = tar_dcvInner$all_model_preds %>%
          dplyr::group_by(model) %>%
          dplyr::summarise(train_auc = roc_auc(., model)) %>%
          data.frame()

        # add in model AUCs on training data
        pmat = rbind(pmat,
                     data.frame(model = tar_dcvInner$Performance$Approach,
                                train_auc = as.numeric(tar_dcvInner$Performance$AUC)))
        pmat$targetFeatures = tar_Features
        # report progress
        message(paste0("Finished fitting DCV model with ", tar_Features, " features"))
        # return list of results
        list(Inner_Perf = perf, Train_AUC = pmat,
             Ridge_Model = tar_dcvInner$glm_model$mdl,
             DCV_Scores = tar_dcvInner$final_dcv)
      } # end target features loop
      # return inner results
      message(paste0("Finished rep ", fld$sd1, ", fold ", f))
      tmp
    } # end nested 5x2-fold CV

  # invert the list of results so we can combine individual parts
  df <- purrr::transpose(lr_dcv_folds)
  # can't row-bind glmnet models, so pull those out
  dcv_perf_ridge_mdls <- df$Ridge_Model
  df$Ridge_Model <- NULL

  dcv_perf_results <- lapply(df, function(l) do.call(rbind, l, fill = TRUE))

  message(paste("NAs in DCV results:", lapply(X = dcv_perf_results, FUN =  function(x) { return(sum(is.na(x))) })))

  ## lambda.1se equivalent ----

  # how many performance estimates did we calculate?
  n_folds <- k_fold*repeats

  # the smallest number of features within 1 std. error of best AUC
  target_1se <- dcv_perf_results$Train_AUC %>%
    dplyr::filter(stringr::str_starts(model, "DCV")) %>%
    dplyr::group_by(model,targetFeatures) %>%
    dplyr::summarise(Mean_AUC = mean(train_auc),
              SE_AUC = sqrt(mean((train_auc - Mean_AUC)^2)/(n_folds - 1))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Upper_SE = Mean_AUC + SE_AUC,
           Max_AUC = max(Mean_AUC),
           In_Range = Upper_SE >= Max_AUC) %>%
    dplyr::filter(In_Range) %>%
    { min(.$targetFeatures) }

  ## optimum number of parts ----

  # the number of features that gave best AUC in DCV models
  train_auc1 = dcv_perf_results$Train_AUC %>%
    dplyr::group_by(model,targetFeatures) %>%
    dplyr::summarise_all(.funs = mean) %>%
    dplyr::filter(stringr::str_detect(model,"DCV"))
  target_max = train_auc1$targetFeatures[which.max(train_auc1$train_auc)]

  message(paste("Optimum model contains", target_max, "features."))
  message(paste("Simplified 1-SE model contains", target_1se, "features."))

  return(list(tt_data = ttData,
              cv_folds = cv_folds,
              train_auc = dcv_perf_results$Train_AUC,
              inner_perf = dcv_perf_results$Inner_Perf,
              inner_dcv_scores = dcv_perf_results$DCV_Scores,
              inner_ridge_models = dcv_perf_ridge_mdls,
              target_max = target_max,
              target_1se = target_1se))
}

