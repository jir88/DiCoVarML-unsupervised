#' Additive- Log Ratio (ALR) Feature Selection using BORUTA / LASSO
#'
#' @importFrom magrittr %>%
#'
#' @param train_data a samples by log ratio matrix for feature discovery on training data
#' @param y_train training set class labels
#' @param test_data a samples by log ratio matrix to subset and return. Note key log-ratios are discovered using training set.
#' @param alr_denom which column to use as ALR denominator
#' @param featureSelectionMethod 1 - Boruta; 2- LASSO ; 3 - None
#' @param impute_factor impute factor multiplicative replacement of zeroes
#' @param num_borutaRuns number of runs of the Boruta algorithm
#' @param glm_alpha glmnet alpha parameter (0,1] where 1-LASSO and (0,1) - elasticnet
#'#'
#' @return A list containing:\tabular{ll}{
#'    \code{train_Data} \tab samples by features (derived from featureSelectionMethod) training data  \cr
#'    \tab \cr
#'    \code{test_data} \tab samples by features (derived from featureSelectionMethod) test data \cr
#'    \tab \cr
#'    }
#'
#' @export
#'
alrFeatureSelection =
  function(train_data,
           y_train,
           test_data,
           alr_denom = 1,
           featureSelectionMethod = 1,
           impute_factor = 1e-7,
           num_borutaRuns = 100,
           glm_alpha=1){

    Decision = NULL

    ## CLose data and impute zeroes
    trainData1 = selEnergyPermR::fastImputeZeroes(train_data,impFactor = impute_factor)
    testData1 = selEnergyPermR::fastImputeZeroes(test_data,impFactor = impute_factor)

    ## Get Log-ratios
    lrs.train = data.frame(compositions::alr(trainData1,ivar = alr_denom))
    lrs.test = data.frame(compositions::alr(testData1,ivar = alr_denom))

    switch (featureSelectionMethod,
            {
              ### Boruta
              b = Boruta::Boruta(x = lrs.train,y = y_train,
                                 doTrace = 2,maxRuns = num_borutaRuns,
                                 getImp = Boruta::getImpExtraGini)
              dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
              keep = dec %>%
                dplyr::filter(Decision!="Rejected")
              kr =as.character(keep$Ratio)
              ## select final ratios
              trainData2 = subset(lrs.train,select = c(kr))
              testData2 = subset(lrs.test,select = c(kr))

            },
            {
              ## LASSO
              cv.lasso <- glmnet::cv.glmnet(as.matrix(lrs.train), y_train, family='binomial', alpha=glm_alpha,
                                            #nfolds = numFold_glm,nlambda = num_lambda,
                                            parallel=TRUE, type.measure='auc')
              # Select Features
              df_coef <- (as.matrix(stats::coef(cv.lasso, s=cv.lasso$lambda.min)))
              # See all contributing variables
              impVar = df_coef[df_coef[, 1] != 0, ]
              impVar_names = names(impVar[-1])
              ## select final ratios
              trainData2 = subset(lrs.train,select = c(impVar_names))
              testData2 = subset(lrs.test,select = c(impVar_names))
            },
            {
              ## ALL ALR
              trainData2 = lrs.train
              testData2 = lrs.test
            }
    )

    return(list(train_data =trainData2,test_data = testData2))

  }

