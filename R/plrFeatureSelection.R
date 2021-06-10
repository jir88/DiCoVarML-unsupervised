#' All Pairwise Log Ratio (PLR) Feature Selection using BORUTA / LASSO
#'
#'Performs feature selection on all PLR using BORUTA or LASSO. Can also return all PLR without feature selection.
#'
#' @importFrom magrittr %>%
#'
#' @param train_data a samples by log ratio matrix for featire discovery on training data
#' @param y_train training set class labels
#' @param test_data a samples by log ratio matrix to subset and return. Note key log-ratios are discovered using training set.
#' @param featureSelectionMethod 1 - Boruta; 2- LASSO ; 3 - None
#' @param impute_factor impute factor multiplicative replacement of zeroes
#' @param num_borutaRuns number of runs of the boruta algorithm
#' @param glm_alpha glmnet alpha parameter (0,1] where 1-LASSO and (0,1) - elasticnet
#' @param glm_family which family to use for logistics regression. Should be a valid glmnet family. default = 'binomial'
#'
#' @return A list containing:\tabular{ll}{
#'    \code{train_Data} \tab samples by features (derived from featureSelectionMethod) training data  \cr
#'    \tab \cr
#'    \code{test_data} \tab samples by features (derived from featureSelectionMethod) test data \cr
#'    \tab \cr
#'    }
#' @export
#'

plrFeatureSelection =
  function(train_data,
           y_train,
           test_data,
           featureSelectionMethod = 1,
           impute_factor  = 1e-7,
           num_borutaRuns = 100,
           glm_family='binomial',glm_alpha=1){

    Decision = NULL

    ## CLose data and impute zeroes
    trainData1 = selEnergyPermR::fastImputeZeroes(train_data,impFactor = impute_factor)
    testData1 = selEnergyPermR::fastImputeZeroes(test_data,impFactor = impute_factor)

    ## Get Log-ratios
    lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = y_train,trainData1))[,-1]
    lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = "test",testData1))[,-1]

    switch (featureSelectionMethod,
            {
              ### Boruta
              b = Boruta::Boruta(x = lrs.train,y = y_train,maxRuns = num_borutaRuns,
                                 doTrace = 2,getImp = Boruta::getImpExtraGini)
              dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
              keep = dec %>%
                dplyr::filter(Decision!="Rejected")
              kr =as.character(keep$Ratio)
              ## select final ratios
              if(length(kr)>2){
                trainData2 = subset(lrs.train,select = c(kr))
                testData2 = subset(lrs.test,select = c(kr))
              }else{
                trainData2 = lrs.train
                testData2 =  lrs.test
              }

            },
            {
              ## LASSO
              cv.lasso <- glmnet::cv.glmnet(as.matrix(lrs.train), y_train, family=glm_family, alpha=glm_alpha,
                                            #nfolds = numFold_glm,nlambda = num_lambda,
                                            parallel=TRUE, type.measure='auc')
              # Select Features
              df_coef <- (as.matrix(stats::coef(cv.lasso, s=cv.lasso$lambda.min)))
              # See all contributing variables
              impVar = df_coef[df_coef[, 1] != 0, ]
              impVar_names = names(impVar[-1])
              ## select final ratios
              if(length(impVar_names)>2){
                trainData2 = subset(lrs.train,select = c(impVar_names))
                testData2 = subset(lrs.test,select = c(impVar_names))
              }else{
                trainData2 = lrs.train
                testData2 =  lrs.test
              }
            },
            {
              ## ALL ALR
              trainData2 = lrs.train
              testData2 = lrs.test
            }
    )

    return(list(train_data =trainData2,test_data = testData2))

  }

