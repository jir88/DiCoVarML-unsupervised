#' Compute Differential Compositional Variation (DCV)
#'
#' A function that impute-zeros, closes data and computes DCV for a set of log-ratios. DCV is a univariate 5-component log-ratio scoring algorithm.
#' Helps to identify important log-ratio feature for classification task.
#'
#' @param train_data a set of log-ratio features to be scored
#' @param impute_factor zero-imputation impute factor for multiplicative replacement #'
#' @param y_train class labels
#' @param num_repeats number of cross-validation repeats for DCV. Default = 1. Can be computationally expensive to repeat multiple times.
#' @param num_folds number of cross-validation folds. Reduces over-fitting.
#' @param rank_orderDCV should the log-ratios be order and number of distinct part be computed. Can increase computational time.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{dcv} \tab Aggregated (across-folds) DCV scores \cr
#'    \tab \cr
#'    \code{raw_dcv} \tab raw DCV non aggregated DCV scores with scores from each component \cr
#'    \tab \cr
#'}
#' @export
#'
#' @seealso \code{\link[diffCompVarRcpp]{dcvScores}}
#'
computeDCV = function(train_data,
                      y_train,
                      num_repeats = 1,
                      num_folds = 5,
                      impute_factor=1e-7,
                      rank_orderDCV = F){

  Status = NULL

  message("Compute Logratios")
  suppressMessages(suppressWarnings({
    ## CLose data and impute zeroes
    trainData1 = selEnergyPermR::fastImputeZeroes(train_data,impFactor = impute_factor)

    ## Compute logratio;s train set computed with
    lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = y_train,trainData1))

  }))

  classes = as.character(unique(y_train))
  message("Compute DCV Scores")
  ## compute DCV
  comb = combinat::combn2(classes)
  dcv.all = data.frame()
  keyFeatures = data.frame()
  for(c in 1:nrow(comb)){
    cx = c(comb[c,])
    cc =  lrs.train %>%
      dplyr::filter(Status %in% cx )
    cc$Status = factor(as.character(cc$Status))
    cc.dcv  = diffCompVarRcpp::dcvScores(logRatioMatrix = cc,includeInfoGain = T,nfolds = num_folds,
                                         numRepeats = num_repeats,
                                         rankOrder = rank_orderDCV )
    dcv.all = rbind(dcv.all,data.frame(Cx = c,Comp1 = cx[1],Comp2 = cx[2],cc.dcv$lrs))
  }
  dcv.all$Cx = factor(dcv.all$Cx)
  dcv.spread = tidyr::spread(dcv.all[,-3:-2],"Cx","rowmean")
  dcv.spread[dcv.spread<0]=0
  xx = apply(data.frame(dcv.spread[,-1]), 2, scale)
  xx = rowMeans(xx)
  dcv = data.frame(Ratio = dcv.spread$Ratio,rowmean = xx)
  return(list(dcv = dcv,raw_dcv = dcv.all))

}
