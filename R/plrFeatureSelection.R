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
              ## penalized regression
              train_control <- caret::trainControl(method="cv",
                                                   repeats = 1,
                                                   number=5,seeds = NULL,
                                                   classProbs = TRUE,
                                                   savePredictions = T,
                                                   allowParallel = TRUE,
                                                   summaryFunction = caret::multiClassSummary
              )

              glm.mdl1 = caret::train(x = as.matrix(lrs.train) ,
                                      y =y_train,
                                      metric = "ROC",
                                      max.depth = 0,
                                      method = "glmnet",
                                      trControl = train_control
              )
              imp  = caret::varImp(glm.mdl1)
              imp  = data.frame(feature = rownames(imp$importance),imp = imp$importance,total = rowSums(imp$importance))
              keep = imp[imp$total>0,]
              keep = keep$feature
              if(length(keep)>2){
                trainData2 = subset(lrs.train,select = c(keep))
                testData2 = subset(lrs.test,select = c(keep))
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

