#' Label Stratified Data Partitioning For K-Fold Cross-Validation
#'
#' Easily creates k train-test splits and stores it partition in a list. Allows for easy retrieval of data partitions and sample identification.
#'
#' @importFrom magrittr %>%
#' @param df  a data.frame with samples x features data.frame with the first column being the class labels
#' @param kfold number of data train-test splits
#' @param seed random seed; enables reproducibility of data partitions
#' @param permuteLabel should the labels be permuted before performing partitions. Useful when determining permutation significance of a classifier.
#'
#' @return a k-element list containing the train-test (xtrain_combinedFolds / xtest_kthFold) splits along with sample IDs (xtrain_IDs / xtest_IDs )
#' @export
#'
#' @examples
#' \dontrun{
#' ## Sample Data Matrix with Binary Classes (A/B)
#' dat = data.frame(CLass =sample(x = c("A","B"),size = 50,replace = T) ,
#'    matrix(runif(50*50),nrow = 50))
#'
#' ## Perform 5-fold Data Paritioning of dat
#' foldData = kfoldDataPartition(df = dat,
#'    kfold = 5,
#'    permuteLabel = F,
#'    seed = 1)
#'
#'  ## get first fold
#'   trainData = foldData[[1]]$xtrain_combinedFolds
#'   trainIDS = foldData[[1]]$xtrain_IDs
#'   testData = foldData[[1]]$xtest_kthFold
#'   testIDs = foldData[[1]]$xtest_IDs
#' }
#'
kfoldDataPartition <-
  function(df,kfold,seed = 08272008,permuteLabel = F){
    set.seed(seed)

    ID = NULL
    fold = NULL

    if(permuteLabel==T){
      df[,1] = sample(df[,1])
    }
    f = caret::createFolds(df[,1],k = kfold,list = FALSE)
    df$fold = f
    df_ = df
    df_ = data.frame(ID = rownames(df_),df_)

    #Partition Data
    foldData = list()
    ph_list = list()

    for (j in 1:kfold){
      xtrain = df%>%
        dplyr::filter(fold != j)%>%
        dplyr::select(-fold)

      xtrain_ID = df_%>%
        dplyr::filter(fold != j)%>%
        dplyr::select(ID,1,2,fold)

      xtest = df%>%
        dplyr::filter(fold == j)%>%
        dplyr::select(-fold)

      xtest_ID = df_%>%
        dplyr::filter(fold == j)%>%
        dplyr::select(ID,1,2,fold)

      ph_list[[1]] = xtrain
      names(ph_list)[1] = "xtrain_combinedFolds"
      ph_list[[2]] = xtest
      names(ph_list)[2] = "xtest_kthFold"
      ph_list[[3]] = xtrain_ID
      names(ph_list)[3] = "xtrain_IDs"
      ph_list[[4]] = xtest_ID
      names(ph_list)[4] = "xtest_IDs"

      foldData[[j]] = ph_list
      names(foldData)[j] = paste("fold_",j,sep = "")
    }

    return(foldData)
  }


#' Leave-One-Dataset-Out (LODO) cross validation partitioning
#'
#' @param df samples by features matrix to be partitioned. First columns must contain class labels.
#' @param dataset_labels vector of dataset label for each sample in df
#' @param seed random seed
#'
#' @return a k-element list containing the train-test (xtrain_combinedFolds / xtest_kthFold) splits along with sample IDs (xtrain_IDs / xtest_IDs )#' @export
#'
lodo_partition = function (df, dataset_labels, seed = 8272008) {

  set.seed(seed)
  ID = NULL
  fold = NULL

  f = as.numeric(as.factor(dataset_labels))
  df$fold = f
  df_ = df
  df_ = data.frame(ID = rownames(df_), df_)
  foldData = list()
  ph_list = list()
  for (j in 1:dplyr::n_distinct(f)) {
    xtrain = df %>% dplyr::filter(fold != j) %>% dplyr::select(-fold)
    xtrain_ID = df_ %>% dplyr::filter(fold != j) %>% dplyr::select(ID,
                                                                   1, 2, fold)
    xtest = df %>% dplyr::filter(fold == j) %>% dplyr::select(-fold)
    xtest_ID = df_ %>% dplyr::filter(fold == j) %>% dplyr::select(ID,
                                                                  1, 2, fold)
    ph_list[[1]] = xtrain
    names(ph_list)[1] = "xtrain_combinedFolds"
    ph_list[[2]] = xtest
    names(ph_list)[2] = "xtest_kthFold"
    ph_list[[3]] = xtrain_ID
    names(ph_list)[3] = "xtrain_IDs"
    ph_list[[4]] = xtest_ID
    names(ph_list)[4] = "xtest_IDs"
    foldData[[j]] = ph_list
    names(foldData)[j] = paste("fold_", j, sep = "")
  }
  return(foldData)
}


#' Leave-All_But-One-Out (LABOO) cross-validation partitioning
#'
#' @param df samples by features matrix to be partitioned. First columns must contain class labels.
#' @param dataset_labels vector of dataset label for each sample in df
#' @param seed random seed
#'
#' @return a k-element list containing the train-test (xtrain_combinedFolds / xtest_kthFold) splits along with sample IDs (xtrain_IDs / xtest_IDs )#' @export
#'
laboo_partition = function (df, dataset_labels, seed = 8272008) {

  set.seed(seed)
  ID = NULL
  fold = NULL

  f = as.numeric(as.factor(dataset_labels))
  df$fold = f
  df_ = df
  df_ = data.frame(ID = rownames(df_), df_)
  foldData = list()
  ph_list = list()
  for (j in 1:dplyr::n_distinct(f)) {
    xtrain = df %>% dplyr::filter(fold == j) %>% dplyr::select(-fold)
    xtrain_ID = df_ %>% dplyr::filter(fold == j) %>% dplyr::select(ID,
                                                                   1, 2, fold)
    xtest = df %>% dplyr::filter(fold != j) %>% dplyr::select(-fold)
    xtest_ID = df_ %>% dplyr::filter(fold != j) %>% dplyr::select(ID,
                                                                  1, 2, fold)
    ph_list[[1]] = xtrain
    names(ph_list)[1] = "xtrain_combinedFolds"
    ph_list[[2]] = xtest
    names(ph_list)[2] = "xtest_kthFold"
    ph_list[[3]] = xtrain_ID
    names(ph_list)[3] = "xtrain_IDs"
    ph_list[[4]] = xtest_ID
    names(ph_list)[4] = "xtest_IDs"
    foldData[[j]] = ph_list
    names(foldData)[j] = paste("fold_", j, sep = "")
  }
  return(foldData)
}





#' Get log-ratios from a list of ratio names
#'
#' @param Ratio Ratio names where numerator and denominator are separated by '___'. Numerator and Denominator names should match column names of raMatrix
#' @param raMatrix a zero-imputed relative abundance or count matrix where the log-ratios will be computed from
#' @param Class class labels for raMAtrix samples.
#'
#' @return a set of logratios corresponding the the Ratio list input provided computed from the raMatrix
#' @export
#'
getLogratioFromList <-
  function(Ratio,raMatrix,Class){
    Ratio = data.frame(Ratio)
    keyRats = tidyr::separate(Ratio,1,into = c("Num","Denom"),sep = "___",remove = F)
    el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
    ad = selEnergyPermR::stageData(raMatrix,labels = Class)
    selEnergyPermR::getLogRatios(ad$allData,el_)
  }

