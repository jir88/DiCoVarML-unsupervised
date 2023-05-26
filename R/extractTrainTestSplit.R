#' Extract Train-Test Split
#'
#'Function to extract partitioned train-test. Extract row metadata, performs zero-imputation and returns labels for train-test data
#'
#' @param foldDataList partition fold data. Usually output from the partition data function
#' @param fold which fold to extract
#' @param maxSparisty percentage of samples with 0 counts for a each feature/column/taxa/etc.
#' @param permLabels should the train data labels be permuted
#' @param extractTelAbunance should data be closed and imputed (bool)? If not, the is returned as-is in the train-test split.
#'
#'
#' @return A list containing:\tabular{ll}{
#'    \code{train_Data} \tab samples by features training data  \cr
#'    \tab \cr
#'    \code{test_data} \tab samples by features test data \cr
#'    \tab \cr
#'    \code{train_ids} \tab row metadata for training samples  \cr
#'    \tab \cr
#'    \code{test_ids} \tab row metadata for test samples  \cr
#'    \tab \cr
#'    \code{y_train} \tab training class/group labels  \cr
#'    \tab \cr
#'    \code{y_test} \tab test class/group labels \cr
#'    \tab \cr
#'    \code{imputeFactor} \tab Multiplicative imputation factor calculated on both training and test data together. \cr
#'    \tab \cr
#'    \code{minCLass} \tab minority class for AUPRC calculations and Class imbalance Problems. *Valid for binary class labels only*  \cr
#'    \tab \cr
#'    \code{majClass} \tab minority class for AUPRC calculations and Class imbalance Problems. *Valid for binary class labels only*   \cr
#'    \tab \cr
#' }
#' @export
#'

extractTrainTestSplit =
  function(foldDataList,fold,
           maxSparisty = .9,permLabels=F,
           extractTelAbunance=T){

    if(extractTelAbunance){
      message("Fold-",fold," Extract Partitions, Estimate Sparisty, Impute Zeroes and Close Data")
    }else{
      message("Fold-",fold," Extract Partitions, Estimate Sparisty, and Return Counts (if Counts Exist) without Imputation")
    }

    suppressMessages(suppressWarnings({
      ### Train Data Estimate Sparsity and Impute Zeroes ####
      trainIDs = foldDataList[[fold]]$xtrain_IDs
      td = foldDataList[[fold]]$xtrain_combinedFolds
      td  = selEnergyPermR::processCompData(td,minPrevalence = maxSparisty)
      impFact = td$impFactor
      minorityClass = td$minClss
      majorityClass = td$majClass
      td = td$processedData
      ### remove samples with total sum after sparity ==0
      bool = !rowSums(td[,-1])==0
      td = td[bool,]
      trainIDs = trainIDs[bool,]
      trainData = td[,-1]
      if(extractTelAbunance){
        trainData = selEnergyPermR::fastImputeZeroes(td[,-1],impFactor = impFact)
      }
      ytrain = factor(td[,1])
      if(permLabels){
        ytrain = sample(ytrain)
      }
      classes = as.character(unique(ytrain))

      ## retrieve test data from kfold partition ####
      testIDS = foldDataList[[fold]]$xtest_IDs
      ts = foldDataList[[fold]]$xtest_kthFold
      ts = data.frame(Status = ts[,1],ts[,-1])
      ts = subset(ts,select = colnames(td))
      ### remove samples with total sum after sparity ==0
      bool = !rowSums(ts[,-1])==0
      ts = ts[bool,]
      testIDS = testIDS[bool,]
      testData = ts[,-1]

      if(extractTelAbunance){
        ##Define Test Data
        testData = selEnergyPermR::fastImputeZeroes(ts[,-1],impFactor = impFact)
      }
      ytest = factor(ts[,1])

    }))


      trainData1 = trainData
      testData1 = testData


    return(dataSplit = list(train_Data = trainData1,test_data = testData1,
                            train_ids = trainIDs,test_ids = testIDS,
                            y_train = ytrain,y_test = ytest,
                            imputeFactor = impFact,
                            minCLass = minorityClass,
                            majClass = majorityClass))
  }
