#' Unsupervised Differential Compositional Variation (DiCoVar) Analysis
#'
#' A function that takes raw data and constructs a synthetic copy by resampling
#' each part. These data sets are transformed into the logratio space, and
#' DCV scores used to select logratios which differ in the real data versus the
#' resampled data.
#'
#' @param X matrix-like object with one row per sample and one column per part
#' @param max_sparsity Parts with a higher proportion of zero/missing values will be dropped before zero imputation.
#' @param scale_data Should data be scaled to unit standard deviation before ridge regression?
#' @param alpha amount of LASSO regression to mix into the final ridge regression model
#' @param seed random seed to use for reproducible results, or NULL to use the current seed
#'
#' @return A list containing:\tabular{ll}{
#'    \code{Performance} \tab How well a glmnet model can distinguish real from fake data. \cr
#'    \code{glm_model} \tab The glmnet model itself. \cr
#'    \code{part_matrix} \tab The closed, zero-imputed part matrix. Includes fake data. \cr
#'    \code{final_dcv} \tab The DCV row means for logratios included in the model. \cr
#'    \code{ridge_pmat} \tab Fitted class probabilities from the glmnet model. \cr
#'    \tab \cr
#'}
#' @export
#'
#' @seealso \code{\link[diffCompVarRcpp]{dcvScores}}
#'
unsupervised.dicovar <- function(X, max_sparsity = 0.9, scale_data = TRUE,
                                 alpha = 0, seed = NULL) {
  # get part of the current random seed if none has been specified
  if(!is.numeric(seed) || is.na(seed)) {
    seed <- .Random.seed[3]
  }
  set.seed(seed)

  # Generate synthetic data ----
  message("Generating synthetic data...")

  # this method samples new data from the columns of X
  synth_data <- sapply(X, function(x) { sample(x, replace = TRUE) })

  # combine fake and real data
  combo_data <- rbind(as.data.frame(X), as.data.frame(synth_data))
  # add labels
  combo_data <- cbind(Status = rep(c("real", "fake"), each = nrow(X)),
                      combo_data)

  # impute zeroes and close data ----
  message("Imputing zeroes...")

  # drop sparse parts and calculate imputation factor
  cd_unsparse <- selEnergyPermR::processCompData(combo_data, minPrevalence = max_sparsity)
  # drop any samples that now have zero sum (no counts at all)
  cd_nz <- cd_unsparse$processedData[rowSums(cd_unsparse$processedData[, -1]) != 0, ]
  # impute zeroes and close data
  combo_x = data.frame(selEnergyPermR::fastImputeZeroes(cd_nz[, -1],
                                                       impFactor = cd_unsparse$impFactor))
  combo_y <- as.vector(cd_unsparse$processedData[, 1])

  # compute log ratios ----
  message("Computing logratios...")

  lrs.combo = selEnergyPermR::calcLogRatio(data.frame(Status = combo_y,
                                                      combo_x))
  # calculate DCV scores ----
  message(paste0("Calculating DCV scores..."))

  cc.dcv = diffCompVarRcpp::dcvScores(logRatioMatrix = lrs.combo,
                                      includeInfoGain = T, nfolds = 1,
                                      numRepeats = 1, rankOrder = F,
                                      seed_ = seed)

  ## Compute Node Strength
  # dcv_strength = DiCoVarML::computeDCVStrength(list(dcv = cc.dcv$lrs))

  # which logratios have mean scores >0?
  pos_dcv_ratios <- cc.dcv$lrs$Ratio[cc.dcv$lrs$rowmean > 0]
  cc.dcv$lrs

  ### Select key ratios from logratio matrix
  glm.combo = subset(lrs.combo, select = pos_dcv_ratios)
  # dcv = dcv[1:w,]

  ## Compute Number of Parts (features)
  cn = colnames(glm.combo)
  uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
  n_parts  = dplyr::n_distinct(uniqueParts)

  ## scale data to unit standard deviation
  if(scale_data){
    # calculate scaling factors on training set
    pp = caret::preProcess(glm.combo,method = "scale")
    glm.combo <- stats::predict(pp, glm.combo)
  }


  # Train Ridge Model ----
  message("Train Ridge Regression Model...")

  # 10-fold CV of ridge regression on standardized training data
  message(colnames(as.matrix(glm.combo)))
  cv.clrlasso <- glmnet::cv.glmnet(x = as.matrix(glm.combo),
                                   y = combo_y,
                                   nfolds = 10,
                                   standardize=F, alpha=alpha,family="binomial")

  # pull optimal model coefficients
  features = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))
  # drop any features with coefficient == 0
  features = features[abs(features)>0]
  # report number of features we ended up with
  message(paste0("Ridge model has ", length(features), " features."))
  # same thing, but keep zero coefficients
  c = as.matrix(stats::coef(cv.clrlasso, s = "lambda.min"))

  ## make predictions using the ridge regression model ----

  # get fitted probabilities for the test set
  p = stats::predict(cv.clrlasso, newx = as.matrix(glm.combo), s = "lambda.min",type = "response")
  # record model and (scaled) data used to train and test it
  model_ = list(mdl = cv.clrlasso, data = list(x = glm.combo, y = combo_y))

  # construct ROC curve with ridge model probabilities
  mroc = pROC::roc(response = combo_y, predictor = as.vector(p),
                   auc = TRUE, ci = TRUE)
  # calculate and report AUC for ridge model
  mroc.dcvlasso = mroc$auc
  mroc$ci
  # record class probabilities from ridge model
  pmat = data.frame(p,1-p);colnames(pmat) = c("real", "fake")

  ## Save Performance of ridge regression model
  perf = data.frame(Approach = "DCV-ridgeRegression",
                    AUC = as.numeric(mroc.dcvlasso),
                    # number of features in model
                    number_parts = n_parts,
                    # number of logratios in model
                    number_ratios = ncol(glm.combo),
                    # number of features before feature selection
                    base_parts = ncol(X))

  # Return results ----
  return(list(Performance = perf, # summary stats for fitted models
              # the ridge regression model and its input data
              glm_model = model_,
              # parts in the DCV matrix, on the simplex, zero imputed
              part_matrix = combo_x,
              # DCV row means for final logratios
              final_dcv = cc.dcv$lrs,
              # fitted class probabilities from the ridge model
              ridge_pmat = pmat))
}
