#' Compute DCV Log Ratio Node Strength
#'
#' Calculates the DCV Strength of each node in a log ratio network
#'
#' @param dcv_mat The DCV log ratio matrix. Should be a matrix a p x 2 matrix with columns being the log ratio names and the second column being the DCV score (rowmean).
#' Usually output from the computeDCV function
#'
#' @return a  numberParts x 2 matrix of node strengths
#'
#' @export
#'
#' @seealso \code{\link{computeDCV}}
computeDCVStrength = function(dcv_mat){
  Str = NULL
  dcvScores = dcv_mat$dcv

  dcvScores$rowmean[dcvScores$rowmean<0] = 0
  trainData.md = caret::preProcess(data.frame(dcvScores$rowmean),
                                   method = "range",rangeBounds = c(0,1) )
  scaledScore = stats::predict(trainData.md, data.frame(dcvScores$rowmean))

  ## Compute Node Strength
  el = data.frame(Ratio = dcvScores$Ratio,Score = scaledScore[,1])
  el = tidyr::separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
  g = igraph::graph_from_edgelist(as.matrix(el[,2:3]))
  igraph::E(g)$weight = el$Score
  nodeStrength = data.frame(Node = names(igraph::strength(g)),
                            Str = igraph::strength(g)) %>%
    dplyr::arrange(dplyr::desc(Str))
  nodeStrength
}




#' Impute and Compute all pairwise Log Ratios (PLR)
#'
#' @param train_data count / relative abundance sample by feature matrix where all PLR(p) are to be computed.
#' @param impute_factor multiplicative replacement factor to replace zeroes
#' @param y_train class labels of data
#'
#' @return a n x p matrix of log ratios
#'
#' @export
#'
computeLRs = function(train_data,impute_factor,y_train){


  message("Compute Logratios")
  suppressMessages(suppressWarnings({
    trainData1 = selEnergyPermR::fastImputeZeroes(train_data,
                                                  impFactor = impute_factor)
    lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = y_train,
                                                        trainData1))
  }))

  lrs.train

}

#' Recursive Feature Elimination (RFE) using random forest feature importance
#'
#' Performs nested RFE using a OOB AUC of a random forest (RF) model. Features are ranked using common RF feature importance.
#' RFE is implemented using the ranger R package where feature importance measures should correspond to valid ranger metrics.
#'
#' @importFrom foreach %do%
#'
#' @param train_ratio training log ratios for feature subset discovery. A sample by log ratio matrix
#' @param test_ratio test log ratios to be subset after feature discovery on training dats
#' @param ytrain training class labels
#' @param metric OOB error metric. Currently c("multiClass.ROC"(default),"AUPRC","AUROC")
#' @param impMeasure RF Feature importance metric.
#' @param sets number of recursive set to consider
#' @param ntrees number of RF trees to grow for each set
#' @param kfold number of folds of nested cross validation
#' @param minPercentFeatReturn min percent of log ratio features to return
#'
#' @return A list containing:\tabular{ll}{
#'    \code{reducedTrainRatios} \tab a n x f (retained log ratios)  matrix with RFE selected log ratios \cr
#'    \tab \cr
#'    \code{reducedTestRatio} \tab a n x f (retained log ratios)  matrix with train set discovered log ratios \cr
#'    \tab \cr
#'    \code{trainPerformance} \tab RFE performance results  \cr
#'    \tab \cr
#'    }
#' @export
#'
#' @seealso \code{\link[ranger]{ranger}}
rfeSelection.ByMetric = function(train_ratio,test_ratio,ytrain,
                                 metric = "multiClass.ROC",kfold = 10,
                                 impMeasure = "impurity_corrected",
                                 sets = 10,ntrees = 2000,
                                 minPercentFeatReturn = .2){

  ## Global Bindings
  r = NULL
  Ratio = NULL
  Imp = NULL
  meanImp = NULL
  AUC  =NULL
  Feats = NULL


  #test matrix to select reduced features from
  train_ratio = data.frame(train_ratio)
  baseRatios = train_ratio
  xtst.pc = data.frame(test_ratio)

  #AUPC
  minorityClass = names(which.min(table(ytrain)))
  majorityClass = unique(ytrain[ytrain!=minorityClass])


  #prevent from training with 0 features
  subsets = round(seq(round(ncol(train_ratio)*minPercentFeatReturn),ncol(train_ratio),length.out = sets))
  if(min(subsets)<=1){
    pos = which(subsets<=1)
    subsets = subsets[-pos]
  }

  ## define RFE sets
  subsets = unique(subsets[sets:1])
  sets = length(subsets)
  subsetAUC = data.frame()
  importanceList = list()

  #random forest recirsive feature selection
  importance.df =c()
  glm.trainAUC = c()
  glm.auprc = c()
  multiClass.roc = c()
  f = caret::createFolds(y = ytrain,k = kfold,list = F)
  rfeCV = foreach::foreach(r=1:kfold,.combine = rbind ) %do% {

    ## input data
    ph.df = data.frame(Status = factor(ytrain),train_ratio)

    suppressMessages(suppressWarnings({


      n=ceiling(sqrt(ncol(train_ratio)))
      bool = f!=r
      rf.ranger = ranger::ranger(formula = Status~.,data = ph.df[bool,],
                                 num.trees = ntrees,
                                 mtry = n,
                                 max.depth = 0,
                                 min.node.size = 1,
                                 write.forest = T,probability = T,
                                 sample.fraction = 1,
                                 importance = impMeasure,
                                 scale.permutation.importance = T)



      if(metric=="multiClass.ROC"){
        #multiclass ROC
        p = stats::predict(rf.ranger,train_ratio[!bool,])
        rocobj = pROC::multiclass.roc(ytrain[!bool],p$predictions)
        multiClass.roc[r] = pROC::auc(rocobj)
      }
      #output
      data.frame(Rep = r,Ratio = names(ranger::importance(rf.ranger)),
                 Imp = ranger::importance(rf.ranger))

    }))

  }


  suppressMessages(suppressWarnings({
    vim = rfeCV %>%
      dplyr::group_by(Ratio) %>%
      dplyr::summarise(meanImp = mean(Imp))

    #select subset ny metric
    if(metric=="AUROC"){
      ph = data.frame(Feats = ncol(train_ratio),AUC = mean(glm.trainAUC),sdAUC = stats::sd(glm.trainAUC))
      subsetAUC = rbind(subsetAUC,ph)
      importanceList[[1]] = vim
    }else if(metric=="AUPRC"){
      ph = data.frame(Feats = ncol(train_ratio),AUC = mean(glm.auprc),sdAUC = stats::sd(glm.auprc))
      subsetAUC = rbind(subsetAUC,ph)
      importanceList[[1]] = vim
    }else if(metric=="multiClass.ROC"){
      ph = data.frame(Feats = ncol(train_ratio),AUC = mean(multiClass.roc),sdAUC = stats::sd(multiClass.roc))
      subsetAUC = rbind(subsetAUC,ph)
      importanceList[[1]] = vim
    }
  }))

  for(i in 2:sets){
    #subset
    sn = dplyr::top_n(importanceList[[i-1]],n = subsets[i],wt = meanImp)
    train_ratio = data.frame(subset(train_ratio, select=as.character(sn$Ratio)))
    importance.df =c()
    glm.trainAUC = c()
    glm.auprc = c()
    multiClass.roc = c()
    f = caret::createFolds(y = ytrain,k = kfold,list = F)
    rfeCV = foreach::foreach(r=1:kfold,.combine = rbind ) %do% {

      ## data
      ph.df = data.frame(Status = factor(ytrain),train_ratio)
      suppressMessages(suppressWarnings({


        n=ceiling(sqrt(ncol(train_ratio)))
        bool = f!=r
        rf.ranger = ranger::ranger(formula = Status~.,data = ph.df[bool,],
                                   num.trees = ntrees,
                                   mtry = n,
                                   max.depth = 0,
                                   min.node.size = 1,
                                   write.forest = T,probability = T,
                                   sample.fraction = 1,
                                   importance = impMeasure,
                                   scale.permutation.importance = T)




        if(metric=="multiClass.ROC"){
          #multiclass ROC
          p = stats::predict(rf.ranger,train_ratio[!bool,])
          rocobj = pROC::multiclass.roc(ytrain[!bool],p$predictions)
          multiClass.roc[r] = pROC::auc(rocobj)
        }
        #output
        data.frame(Rep = r,Ratio = names(ranger::importance(rf.ranger)),
                   Imp = ranger::importance(rf.ranger))

      }))
    }



    suppressMessages(suppressWarnings({
      vim = rfeCV %>%
        dplyr::group_by(Ratio) %>%
        dplyr::summarise(meanImp = mean(Imp))

      #select subset ny metric
      if(metric=="AUROC"){
        ph = data.frame(Feats = ncol(train_ratio),AUC = mean(glm.trainAUC),sdAUC = stats::sd(glm.trainAUC))
        subsetAUC = rbind(subsetAUC,ph)
        importanceList[[i]] = vim
      }else if(metric=="AUPRC"){
        ph = data.frame(Feats = ncol(train_ratio),AUC = mean(glm.auprc),sdAUC = stats::sd(glm.auprc))
        subsetAUC = rbind(subsetAUC,ph)
        importanceList[[i]] = vim
      }else if(metric=="multiClass.ROC"){
        ph = data.frame(Feats = ncol(train_ratio),AUC = mean(multiClass.roc),sdAUC = stats::sd(multiClass.roc))
        subsetAUC = rbind(subsetAUC,ph)
        importanceList[[i]] = vim
      }
    }))

    message(paste("Subset-",subsets[i]," ",i,"/",sets," Calculated............",sep = ""))
  }

  #select best features
  subsetAUC$i = 1:nrow(subsetAUC)
  subsetAUC = subsetAUC %>%
    dplyr::top_n(1,AUC) %>%
    dplyr::arrange(dplyr::desc(Feats))
  i = subsetAUC$i[1]
  varImp_ = importanceList[[i]]
  train_ratio = data.frame(subset(baseRatios, select=as.character(varImp_$Ratio)))
  xtst.pc1 = data.frame(subset(xtst.pc, select=as.character(varImp_$Ratio)))
  subsetAUC$CI = (1.96*(subsetAUC$sdAUC))/sqrt(kfold)
  subsetAUC$lower = subsetAUC$AUC - subsetAUC$CI
  subsetAUC$upper = subsetAUC$AUC + subsetAUC$CI



  return(list(reducedTrainRatios = train_ratio,
              reducedTestRatio = xtst.pc1,
              trainPerformance = subsetAUC)
  )
}






#' Performs DCV Hybrid Feature Selection via Hub Node Significance
#'
#' Selects 2 subset of log ratios containing statistically significant hub nodes (permutation testing). Subsets returned are first in three stages.
#' In stage 1 log-ratios are scored using DCV. Afterwards hub 'parts/taxa/etc.' are identified within the log ratio network and significant (permutation thresholded) log ratios containing
#' at least one hub are retained. Finally in stages 2 & 3 two rounds of feature selection is carried out the final subsets. Theses subset are:
#' \itemize{
#'  \item{"MST-set"}{Non Redundant Subset of logratios pruned after stage 1}
#'  \item{"Dense"}{ All log ratios from stage 1 }
#' }
#'
#' @param xtrain A samples(n) by log ratios(p) matrix to be used for feature selection
#' @param ytrain Class labels vector of length n
#' @param xtest A samples (xn)by log ratios (xp) matrix to subset after features are discovered on the xtrain matrix.
#' @param impute_factor multiplicative factor for imputed zero counts/abundance
#' @param nperms number of permutation for hub node discovery
#' @param nfold_dcv number of cross validation folds for DCV computation. Default is 1(all data). Increasing folds increases computational time exponentially.
#' @param alpha significance threshold for hub node identification. Default is 0.10. Lower alpha returns fewer hubs and vice-versa.
#' @param maxBorutaRuns max number of Boruta runs
#' @param num_treesRFE number of Random Forest trees to be grown during recursive feature elimination (RFE). Default = 2,000
#' @param num_setRFE number of sets in RFE. Default = 20
#' @param verbose should status of RFE and Boruta be displayed
#' @param rfImportance RF feature importance metric to be used with RFE. Requires a valid 'ranger::ranger' importance measure
#' @param k_fold number of cross-validation folds for RFE
#' @param minPercFeat minimum number of feature to be considered in RFE
#' @param useGLM should penalized regression be used instead of hybrid
#'
#' @return A list containing:\tabular{ll}{
#'    \code{MST} \tab a n x f (retained log ratios) derived from MST pruned stage-1 log ratio network \cr
#'    \tab \cr
#'    \code{Dense} \tab a n x f (retained log ratios)  derived from dense stage-1 network  \cr
#'    \tab \cr
#'    }
#' @export
#'
hybrid_dcvfeatureSelection = function(xtrain,ytrain,xtest = NULL,impute_factor = 1e-7,useGLM = F,
                                 nperms = 20,nfold_dcv = 1,alpha = 0.1,
                                 num_treesRFE = 2000,num_setRFE = 20,verbose = 0,
                                 rfImportance = "impurity",k_fold = 10,minPercFeat = 0.2,
                                 maxBorutaRuns = 200){

  ## Global Bindings
  score = NULL
  num = NULL
  denom = NULL
  Ratio = NULL
  rowmean = NULL
  Decision = NULL

  message("Compute Log Ratios")
  suppressMessages(suppressWarnings({
    lrs_train = computeLRs(train_data = xtrain,impute_factor = impute_factor,y_train = ytrain)
    if(is.null(xtest)){
      lrs_test = lrs_train
    }else{
      lrs_test = computeLRs(train_data = xtest,impute_factor = impute_factor,y_train = "test")
    }
  }))

  message("Identify Log Ratio Network Hubs")
  pb <- progress::progress_bar$new(
    format = " - Compute DCV Strength Null Distribution [:bar] :percent eta: :eta",
    total = nperms, clear = FALSE, width= 60)
  strength.df = data.frame()
  dcvScore = data.frame()
  for(s in 1:nperms){
    set.seed(s)
    yt  = sample(ytrain)
    suppressMessages(suppressWarnings({
      dcv_mat = computeDCV(train_data = xtrain,lrs.train = lrs_train,
                            y_train = yt,num_folds = nfold_dcv,
                            impute_factor =impute_factor)
      ph = dcv_mat$dcv
      ph$seed = s
      dcvScore = rbind(dcvScore,ph)
      ss = computeDCVStrength(dcv_mat)
      ss$seed = s
      strength.df = rbind(strength.df,ss)
    }))
    pb$tick()
  }
  strength.df = tidyr::spread(strength.df,"seed","Str")
  strength.df = data.frame(strength.df)
  dcvScore = tidyr::spread(dcvScore,"seed","rowmean")
  dcvScore = data.frame(dcvScore)

  message(" - Compute Empirical DCV Strength")
  suppressMessages(suppressWarnings({
    dcv_mat = computeDCV(train_data = xtrain,lrs.train = lrs_train,
                         y_train = ytrain,num_folds = nfold_dcv,
                         impute_factor = impute_factor)
    trueStrength = computeDCVStrength(dcv_mat)
    pre.df = dplyr::left_join(trueStrength,strength.df)
    permCheck = function(x){
      sum(x[1]<x[2:length(x)])
    }
    trueStrength = computeDCVStrength(dcv_mat)
    trueScore = dcv_mat$dcv
    pre.df = dplyr::left_join(trueStrength,strength.df)
    score.df = dplyr::left_join(trueScore,dcvScore)
    permCheck = function(x){
      sum(x[1]<x[2:length(x)])
    }
    pre.df$score = apply(pre.df[,-1], 1, permCheck)
    score.df$score = apply(score.df[,-1], 1, permCheck)

    th = round(nperms*alpha)
    keep = pre.df %>%
      dplyr::filter(score<th)
    keep_ratios = score.df%>%
      dplyr::filter(score<th)
  }))


  message("Prune Log Ratio Network")
  ratios = data.frame()
  el = data.frame(Ratio = dcv_mat$dcv$Ratio)
  el = tidyr::separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
  for(r in keep$Node){
    ph = el %>%
      dplyr::filter(num==r | denom==r)
    ph = dplyr::left_join(ph,dcv_mat$dcv,by = "Ratio")
    ratios = rbind(ratios,ph)
  }
  ## Prune Log Ratio Network
  ratios = dplyr::distinct(ratios) %>%
    dplyr::filter(Ratio %in% keep_ratios$Ratio)
  ratios = dplyr::distinct(ratios) %>%
    dplyr::filter(rowmean>0)
  ## Select Full Subset
  train_x = subset(lrs_train,select = ratios$Ratio)
  test_x = subset(lrs_test,select = ratios$Ratio)

  if(useGLM){
    message("penalized glm Log Ratio Selection from Dense Network")
    train_x1 = train_x
    test_x1 = test_x
    ## penalized regression
    train_control <- caret::trainControl(method="cv",
                                         repeats = 1,
                                         number=5,seeds = NULL,
                                         classProbs = TRUE,
                                         savePredictions = T,
                                         allowParallel = TRUE,
                                         summaryFunction = caret::multiClassSummary
    )

    glm.mdl1 = caret::train(x = as.matrix(train_x1) ,
                            y =ytrain,
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
      train_data2 = subset(train_x1,select = c(keep))
      test_data2 = subset(test_x1,select = c(keep))
    }else{
      pp = rfeSelection.ByMetric(train_ratio = train_x1,test_ratio = test_x1,
                                 ytrain =ytrain,ntrees = num_treesRFE,sets = num_setRFE,
                                 impMeasure = rfImportance,kfold = k_fold,
                                 minPercentFeatReturn = minPercFeat)
      train_data2 = pp$reducedTrainRatios
      test_data2 = pp$reducedTestRatio
    }


    message("penalized Log Ratio Selection from MST")
    train_x1 = diffCompVarRcpp::mstAll(train_x,ratios)
    test_x1 = diffCompVarRcpp::mstAll(test_x,ratios)
    glm.mdl1 = caret::train(x = as.matrix(train_x1) ,
                            y =ytrain,
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
      train_data1 = subset(train_x1,select = c(keep))
      test_data1 = subset(test_x1,select = c(keep))
    }else{
      pp = rfeSelection.ByMetric(train_ratio = train_x1,test_ratio = test_x1,
                                 ytrain =ytrain,ntrees = num_treesRFE,sets = num_setRFE,
                                 impMeasure = rfImportance,kfold = k_fold,
                                 minPercentFeatReturn = minPercFeat)
      train_data1 = pp$reducedTrainRatios
      test_data1 = pp$reducedTestRatio
    }
  }else{
    message("Hybrid Log Ratio Selection from Dense Network")
    train_x1 = train_x
    test_x1 = test_x
    message(" - Compute Stage-1")
    b = Boruta::Boruta(x = train_x1,y = ytrain,doTrace = verbose,
                       maxRuns = maxBorutaRuns,
                       getImp = Boruta::getImpExtraGini)
    dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
    keep = dec %>%
      dplyr::filter(Decision!="Rejected")
    kr =as.character(keep$Ratio)
    if(length(kr)>2){
      train_data = subset(train_x1,select = c(kr))
      test_data = subset(test_x1,select = c(kr))
    }else{
      train_data = train_x1
      test_data = test_x1
    }

    message(" - Compute Stage-2")
    if(verbose==0){
      suppressMessages(suppressWarnings({
        pp = rfeSelection.ByMetric(train_ratio = train_data,test_ratio = test_data,
                                   ytrain =ytrain,ntrees = num_treesRFE,sets = num_setRFE,
                                   impMeasure = rfImportance,kfold = k_fold,
                                   minPercentFeatReturn = minPercFeat)
      }))
    }else{
      pp = rfeSelection.ByMetric(train_ratio = train_data,test_ratio = test_data,
                                 ytrain =ytrain,ntrees = num_treesRFE,sets = num_setRFE,
                                 impMeasure = rfImportance,kfold = k_fold,
                                 minPercentFeatReturn = minPercFeat)

    }
    train_data2 = pp$reducedTrainRatios
    test_data2 = pp$reducedTestRatio


    message("Hybrid Log Ratio Selection from MST")
    train_x1 = diffCompVarRcpp::mstAll(train_x,ratios)
    test_x1 = diffCompVarRcpp::mstAll(test_x,ratios)
    message(" - Compute Stage-1")
    b = Boruta::Boruta(x = train_x1,y = ytrain,doTrace = verbose,
                       maxRuns = maxBorutaRuns,
                       getImp = Boruta::getImpExtraGini)
    dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
    keep = dec %>%
      dplyr::filter(Decision!="Rejected")
    kr =as.character(keep$Ratio)
    if(length(kr)>2){
      train_data = subset(train_x1,select = c(kr))
      test_data = subset(test_x1,select = c(kr))
    }else{
      train_data = train_x1
      test_data = test_x1
    }
    message(" - Compute Stage-2")
    if(verbose==0){
      suppressMessages(suppressWarnings({
        pp = rfeSelection.ByMetric(train_ratio = train_data,test_ratio = test_data,
                                   ytrain =ytrain,ntrees = num_treesRFE,sets = num_setRFE,
                                   impMeasure = rfImportance,kfold = k_fold,
                                   minPercentFeatReturn = minPercFeat)
      }))
    }else{
      pp = rfeSelection.ByMetric(train_ratio = train_data,test_ratio = test_data,
                                 ytrain =ytrain,ntrees = num_treesRFE,sets = num_setRFE,
                                 impMeasure = rfImportance,kfold = k_fold,
                                 minPercentFeatReturn = minPercFeat)

    }
    train_data1 = pp$reducedTrainRatios
    test_data1 = pp$reducedTestRatio
  }




  return(list(MST = list(train = train_data1,test = test_data1),
              Dense = list(train = train_data2,test = test_data2)
              )
         )
}







#' Performs DCV Ratio Filtering
#'
#' Selects 2 subset of log ratios based dcvScore thresholding. Theses subset are:
#' \itemize{
#'  \item{"MST-set"}{ - Non Redundant Subset of logratios pruned after thresholding}
#'  \item{"Dense"}{ - All log ratios after thresholding }
#' }
#'
#' @param xtrain A samples(n) by log ratios(p) matrix to be used for feature selection
#' @param ytrain Class labels vector of length n
#' @param xtest A samples (xn)by log ratios (xp) matrix to subset after features are discovered on the xtrain matrix.
#' @param impute_factor multiplicative factor for imputed zero counts/abundance
#' @param th_percent dcv score percentile for thresholding i.e. keep ratio where dcvScore > threshold(th_percent)
#' @param nfold_dcv  number of cross validation folds for DCV computation. Default is 1(all data). Increasing folds increases computational time exponentially.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{MST} \tab a n x f (retained log ratios) derived from MST  \cr
#'    \tab \cr
#'    \code{Dense} \tab a n x f (retained log ratios)  derived from dense thresholding \cr
#'    \tab \cr
#'    }
#' @export
#'
dcvRatioFilter = function(xtrain,ytrain,xtest = NULL,
                          impute_factor = 1e-7,
                          nfold_dcv = 1,
                          th_percent=.5){

  ## Global Bindings
  score = NULL
  num = NULL
  denom = NULL
  Ratio = NULL
  rowmean = NULL
  Decision = NULL
  nfold_dcv = NULL

  message("Compute Log Ratios")
  suppressMessages(suppressWarnings({
    lrs_train = computeLRs(train_data = xtrain,impute_factor = impute_factor,y_train = ytrain)
    if(is.null(xtest)){
      lrs_test = lrs_train
    }else{
      lrs_test = computeLRs(train_data = xtest,impute_factor = impute_factor,y_train = "test")
    }
  }))


  message(" - Compute Empirical DCV Strength")
  suppressMessages(suppressWarnings({
    dcv_mat = computeDCV(train_data = xtrain,lrs.train = lrs_train,
                         y_train = ytrain,num_folds = nfold_dcv,
                         impute_factor = impute_factor)

  }))

  trueScore = dcv_mat$dcv
  trueScore[trueScore$rowmean<0,] = 0
  th = stats::quantile(trueScore$rowmean,probs = th_percent)
  ratios =trueScore %>%
    dplyr::filter(rowmean>th)

  ## Select Full Subset
  train_x = subset(lrs_train,select = ratios$Ratio)
  test_x = subset(lrs_test,select = ratios$Ratio)

  train_x1 = diffCompVarRcpp::mstAll(train_x,ratios)
  test_x1 = diffCompVarRcpp::mstAll(test_x,ratios)



  return(list(MST = list(train = train_x1,test = test_x1),
              Dense = list(train = train_x,test = test_x)
  )
  )
}
