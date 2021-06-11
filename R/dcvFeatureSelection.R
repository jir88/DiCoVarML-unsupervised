

#' Differential Compositional Variation (DCV) Guided Feature Selection on all Pairwise Log Ratios
#'
#' Performs DCV derived feature selection using three approaches:
#' \itemize{
#'  \item{"Approach-1"}{ Targeted featutre selection - (Forward Selection)finds best subset given constraints on the number parts/taxa}
#'  \item{"DEPRECATED Approach-2"}{ Untargeted network based selection - (Recursive Selection)finds reduced subsets with no constran on numer of parts}
#'  \item{"DEPRECATED Approach-3"}{ Untargeted DCV guided Boruta - reduces computational time of BORUTA Selection with DCV filtering}
#' }
#'
#'  @importFrom magrittr %>%
#'  @importFrom foreach %dopar%#'
#'
#' @param train_data a samples by log ratio matrix for feature discovery on training data
#' @param y_train training set class labels
#' @param test_data a samples by log ratio matrix to subset and return. Note key log-ratios are discovered using training set.
#' @param impute_factor impute factor multiplicative replacement of zeroes
#' @param dcv_matrix Default = NULL; DCV matrix can be supplied
#' @param dcv_method 1:  Approach-1; 2:  Approach-2; 3:  Approach-3
#' @param percent_topDCV Using DCV scores threshold for percent of top log ratios to return (Approach-3)
#' @param num_sets number of sets for recursive selection (Approach-2)
#' @param tarFeats targeted parts/taxa/etc. to select; (Approach-1)
#' @param minFeats Minimum number of log ratio set to consider (Approach-1/2)
#' @param maxBorutaRuns number of runs of the Boruta algorithm (Approach-3)
#' @param nruns_rfAUC Number of runs for random forest AUC estimate (Approach-2)
#' @param num_trees number of trees random forest (Approach-2)
#' @param upperBound_percent Max percent of parts/taxa to be considered. Upperbound  = upperBound_percent*total number of parts where sets = seq(upperbound,minFeats,length.out = num_sets) (Approach-2)
#' @param rf_importance random forest importance measure (Approach-2)
#'
#' @return A list containing:\tabular{ll}{
#'    \code{train_Data} \tab samples by features (derived from dcv_method) training data  \cr
#'    \tab \cr
#'    \code{test_data} \tab samples by features (derived from dcv_method) test data \cr
#'    \tab \cr
#'    }
#'
#' @export
#'
dcvFeatureSelection = function(train_data,
                               y_train,
                               test_data,
                               impute_factor = 1e-7,
                               dcv_matrix = NULL,
                               dcv_method=2,
                               percent_topDCV = .95,
                               num_sets = 10,
                               tarFeats,
                               minFeats=4,
                               maxBorutaRuns = 100,
                               nruns_rfAUC = 20,
                               num_trees = 1000,upperBound_percent = .5,
                               rf_importance = "impurity_corrected"){

  ## Global Bindings
  Decision = NULL
  Str = NULL
  Ratio = NULL
  Status = NULL
  num = NULL
  denom = NULL
  Score = NULL
  Score.y = NULL

  message("Compute Logratios")
  suppressMessages(suppressWarnings({
    ## CLose data and impute zeroes
    trainData1 = selEnergyPermR::fastImputeZeroes(train_data,impFactor = impute_factor)
    testData1 = selEnergyPermR::fastImputeZeroes(test_data,impFactor = impute_factor)
    ## Compute logratio;s train set computed with
    ## trainData1(adjusted Ratios if covariate matrix exist)
    lrs.train = selEnergyPermR::calcLogRatio(data.frame(Status = y_train,trainData1))
    lrs.test = selEnergyPermR::calcLogRatio(data.frame(Status = "test",testData1))
  }))

  classes = as.character(unique(y_train))

  if(is.null(dcv_matrix)){
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
      cc.dcv  = diffCompVarRcpp::dcvScores(logRatioMatrix = cc,includeInfoGain = T,nfolds = 5,numRepeats = 1,rankOrder = F )
      kf =  diffCompVarRcpp::mstAll(featMatrix = lrs.train,dcvRanking = cc.dcv$lrs)
      keyFeatures = rbind(keyFeatures,data.frame(Cx = c,Comp1 = cx[1],Comp2 = cx[2],Ratios = colnames(kf)))
      dcv.all = rbind(dcv.all,data.frame(Cx = c,Comp1 = cx[1],Comp2 = cx[2],cc.dcv$lrs))
    }
    unKeyFeatures = unique(keyFeatures$Ratios)
    dcv.all$Cx = factor(dcv.all$Cx)
    dcv.spread = tidyr::spread(dcv.all[,-3:-2],"Cx","rowmean")
    dcv.spread[dcv.spread<0]=0
    xx = apply(data.frame(dcv.spread[,-1]), 2, scale)
    xx = rowMeans(xx)
    dcv = data.frame(Ratio = dcv.spread$Ratio,rowmean = xx)
  }else{
    dcv = dcv_matrix
  }



  switch (dcv_method,
          {message("Perform Targeted Feature Selection. Retain ",tarFeats,"-Parts")},
          {message("Perform Optimal Tree-Based Feature Selection")},
          {message("Perform DCV Guided Boruta Feature Selection")}
  )
  switch (dcv_method,
          {  ### 1 - Targeted Feature Selection  ######

            suppressMessages(suppressWarnings({
              ## Scale Scores
              dcv$rowmean[dcv$rowmean<0] = 0
              trainData.md = caret::preProcess(data.frame(dcv$rowmean),
                                               method = "range",rangeBounds = c(0,1) )
              scaledScore = stats::predict(trainData.md, data.frame(dcv$rowmean))

              ## Compute Node Strength
              el = data.frame(Ratio = dcv$Ratio,Score = scaledScore[,1])
              el = tidyr::separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
              g = igraph::graph_from_edgelist(as.matrix(el[,2:3]))
              igraph::E(g)$weight = el$Score
              nodeStrength = data.frame(Node = names(igraph::strength(g)),
                                        Str = igraph::strength(g)) %>%
                dplyr::arrange(dplyr::desc(Str))


              nodeStrength1 = nodeStrength %>%
                dplyr::top_n(n = tarFeats,wt = Str) %>%
                dplyr::arrange(dplyr::desc(Str))
              nodeStrength1 = nodeStrength1[1:tarFeats,]

              featureList = list()
              nodePerf = data.frame()
            }))


            ## Select Features From Network
            ra.train = subset(trainData1,select = nodeStrength1$Node[1:tarFeats])
            ra.test = subset(testData1,select = nodeStrength1$Node[1:tarFeats])
            suppressMessages(suppressWarnings({
              ph.train = selEnergyPermR::calcLogRatio(data.frame(Status = "test",ra.train))[,-1]
              ph.test = selEnergyPermR::calcLogRatio(data.frame(Status = "test",ra.test))[,-1]
            }))
            ac = foreach::foreach(i = 1:nruns_rfAUC,.combine = rbind)%dopar%{
              testh  = ranger::ranger(formula = Status~.,
                                      data = data.frame(Status = y_train,ph.train),
                                      importance = rf_importance,
                                      replace = T,
                                      num.trees = num_trees,
                                      probability = T)
              vi = testh$variable.importance
              testH = pROC::auc(pROC::multiclass.roc(y_train,testh$predictions))
              data.frame(Seed = i,Ratio = names(vi) , rowmean = as.numeric(vi),AUC = as.numeric(testH))
            }
            ac = ac %>%
              dplyr::group_by(Ratio) %>%
              dplyr::summarise_all(.funs = mean)
            testH = unique(ac$AUC);testH
            imp.df = data.frame(ac[,c(-2,-4)])

            ## CONSTRUCT NETWORK
            el = data.frame(Ratio = imp.df$Ratio,Score = imp.df$rowmean)
            el = tidyr::separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
            ## add diagonal
            parts = unique(c(el$num,el$denom))
            el_ = rbind(el,data.frame(Ratio = paste0(parts,"___",parts),num = parts,denom = parts,Score=0))
            dcv_adj = tidyr::spread(el_[,-1],"num","Score",fill = 0)
            rownames(dcv_adj) = dcv_adj[,1]
            dcv_adj = dcv_adj[,-1]
            dcv_adj = knnADJtoSYM((dcv_adj))

            ## determine knn sets
            end = tarFeats-1
            knns = round(seq(1,end,length.out = num_sets))
            knns = unique(knns[knns<=nrow(nodeStrength1)] )
            knn_search = data.frame()
            feature_list = list()
            imp_list = list()
            flp = 1
            pb <- progress::progress_bar$new(
              format = " KNN-Network Selection [:bar] :percent eta: :eta",
              total = length(knns), clear = FALSE, width= 60)
            for(k in knns ){

              suppressMessages(suppressWarnings({
                g = simplexDataAugmentation::knn_graph(dcv_adj,K = k,sim_ = F)
                g = g$Graph
                el.knn = data.frame(igraph::get.edgelist(g))
                colnames(el.knn) = c("num","denom")
                el.knn$Ratio = paste0(el.knn$num,"___",el.knn$denom)
                el.knn = dplyr::left_join(el.knn,el)
                el.knn1 = el.knn %>%
                  dplyr::filter(!is.na(Score))
                el.knn2 = el.knn %>%
                  dplyr::filter(is.na(Score)) %>%
                  dplyr::mutate(Ratio = paste0(denom,"___",num)) %>%
                  dplyr::left_join(el,by = "Ratio") %>%
                  dplyr::filter(!is.na(Score.y))
                keep = c(el.knn1$Ratio,el.knn2$Ratio)
                ## Save features
                trainData2 = subset(ph.train,select = keep)
                testData2 = subset(ph.test,select = keep)


                ## caret fold approach
                # ph = DiCoVarML::trainML_Models(trainLRs =trainData2,
                #                                testLRs = trainData2,
                #                                ytrain = y_train,
                #                                y_test = y_train,
                #                                testIDs = NULL,
                #                                ntrees = num_trees,
                #                                numFolds = 5,
                #                                numRepeats = 5,
                #                                mtry_ = round(sqrt(ncol(trainData2))),
                #                                ranger_imp = rf_importance,
                #                                models = "ranger" )
                # ph.perf = data.frame(K = k,ratios = ncol(trainData2),AUC = ph$performance$TrainAUC)
                # ## variable importance
                # vi = ranger::importance(ph$models[[1]]$finalModel)
                # imp.df = data.frame(Ratio =names(vi),rowmean = as.numeric(vi))

                ## oob randomfrest approach
                ac = foreach::foreach(i = 1:nruns_rfAUC,.combine = rbind)%dopar%{
                  testh  = ranger::ranger(formula = Status~.,
                                          data = data.frame(Status = y_train,trainData2),
                                          importance = rf_importance,
                                          replace = T,
                                          num.trees = num_trees,
                                          probability = T)
                  vi = testh$variable.importance
                  testH = pROC::auc(pROC::multiclass.roc(y_train,testh$predictions))
                  data.frame(Seed = i,Ratio = names(vi) , rowmean = as.numeric(vi),AUC = as.numeric(testH))
                }
                ac = ac %>%
                  dplyr::group_by(Ratio) %>%
                  dplyr::summarise_all(.funs = mean)
                testH = unique(ac$AUC);testH
                imp.df = data.frame(ac[,c(-2,-4)])
                ph.perf = data.frame(K = k,ratios = igraph::ecount(g),AUC = testH)

                ## Store results
                imp_list[[flp]] = imp.df
                feature_list[[flp]] = list(train = trainData2,test = testData2)
                flp = flp+1
                knn_search = rbind(knn_search,ph.perf)
              }))

              pb$tick()

            }

            print(knn_search)
            kk = which.max(knn_search$AUC)[1]
            trainData2 = feature_list[[kk]]$train
            testData2 = feature_list[[kk]]$test
          },

          {## 2 - Optimal Tree-Based Feature Selection


            suppressMessages(suppressWarnings({

              ## Scale Scores
              scr = dcv
              trainData.md = caret::preProcess(data.frame(scr$rowmean),method = "range",rangeBounds = c(0,1) )
              scaledScore = stats::predict(trainData.md, data.frame(scr$rowmean))

              ##node strength
              el = data.frame(Ratio = scr$Ratio,Score = scaledScore[,1])
              el = tidyr::separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
              g = igraph::graph_from_edgelist(as.matrix(el[,2:3]))
              igraph::E(g)$weight = el$Score
              nodeStrength = data.frame(Node = names(igraph::strength(g)),Str = igraph::strength(g)) %>%
                dplyr::arrange(dplyr::desc(Str))

              sets = round(seq(nrow(nodeStrength)*upperBound_percent,minFeats,length.out = num_sets))
              sets = unique(sets)
              ## Select All Features
              nodeStrength1 = nodeStrength %>%
                dplyr::top_n(n = sets[1],wt = Str)
              featureList = list()
              nodePerf = data.frame()
              imp.df = scr
              unFeats = imp.df$Ratio

            }))

            ## Setup progress bar
            pb <- progress::progress_bar$new(
              format = " Recursively selecting features [:bar] :percent eta: :eta",
              total = length(sets), clear = FALSE, width= 60)

            ## Estimate Performance
            for(i in 1:length(sets)){
              ph = nodeStrength1[1:sets[i],]

              suppressMessages(suppressWarnings({
                imp.df$rowmean[imp.df$rowmean<0] = 0
                el = data.frame(Ratio = imp.df$Ratio,Score = imp.df$rowmean)
                el = tidyr::separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
                g = igraph::graph_from_edgelist(as.matrix(el[,2:3]))
                igraph::E(g)$weight = el$Score
                ss = data.frame(s = igraph::strength(g))
                dcv_adj = igraph::as_adjacency_matrix(graph = g,sparse = F,
                                                      attr = "weight")
                dcv_adj = knnADJtoSYM(dcv_adj)

                bool = colnames(dcv_adj)%in%ph$Node
                dcv_adj = dcv_adj[bool,bool]

                end = tarFeats-1
                knns = round(seq(1,end,length.out = num_sets))
                knns = unique(knns[knns<=nrow(nodeStrength1)] )
                knn_search = data.frame()
                feature_list = list()
                imp_list = list()
                flp = 1
                for(k in knns ){
                  g = simplexDataAugmentation::knn_graph(dcv_adj,K = k,sim_ = F)
                  g = g$Graph

                  # plot(g, main = "test",
                  #      layout = igraph::layout.fruchterman.reingold,
                  #      vertex.size = abs(log((igraph::strength(g)+1)))+3,
                  #      vertex.label = NA,
                  #      vertex.label.cex = 0.75,
                  #      edge.curved = 0.2,
                  #      edge.width = igraph::E(g)$weight * 0.5,
                  #      #edge.color = cols,
                  #      #edge.arrow.size = 0.5,
                  #      #edge.arrow.width = 1
                  # )


                  el.knn = data.frame(igraph::get.edgelist(g))
                  colnames(el.knn) = c("num","denom")
                  el.knn$Ratio = paste0(el.knn$num,"___",el.knn$denom)
                  el.knn = dplyr::left_join(el.knn,el)
                  trainData2 = getLogratioFromList(Ratio = el.knn$Ratio,raMatrix = trainData1,Class = y_train)
                  testData2 = getLogratioFromList(Ratio = el.knn$Ratio,raMatrix = testData1,Class = "test")



                  ac = foreach::foreach(i = 1:nruns_rfAUC,.combine = rbind)%dopar%{
                    testh  = ranger::ranger(formula = Status~.,
                                            data = data.frame(Status = y_train,trainData2),
                                            importance = rf_importance,
                                            replace = T,
                                            num.trees = num_trees,
                                            probability = T)
                    vi = testh$variable.importance
                    testH = pROC::auc(pROC::multiclass.roc(y_train,testh$predictions))
                    data.frame(Seed = i,Ratio = names(vi) , rowmean = as.numeric(vi),AUC = as.numeric(testH))
                  }
                  ac = ac %>%
                    dplyr::group_by(Ratio) %>%
                    dplyr::summarise_all(.funs = mean)
                  testH = unique(ac$AUC);testH
                  imp.df = data.frame(ac[,c(-2,-4)])

                  imp_list[[flp]] = imp.df
                  feature_list[[flp]] = list(train = trainData2,test = testData2)
                  flp = flp+1
                  ph.perf = data.frame(K = k,ratios = igraph::ecount(g),AUC = testH)
                  knn_search = rbind(knn_search,ph.perf)

                }

                kk = which.max(knn_search$AUC)[1]
                trainData2 = feature_list[[kk]]$train
                testData2 = feature_list[[kk]]$test
                cn = colnames(trainData2)
                numberRatios = length(cn)
                uniqueParts = unique(as.vector(stringr::str_split(cn,"___",2,simplify = T)))
                numUniqPart  = dplyr::n_distinct(uniqueParts)
                testH = knn_search$AUC[kk]
                nodePerf = rbind(nodePerf,data.frame(Set = i,
                                                     numParts = numUniqPart,
                                                     numRatios = numberRatios,
                                                     AUC = as.numeric(testH)))
                imp.df = imp_list[[kk]]



                #imp.df$rowmean[imp.df$rowmean<0] = 0
                trainData.md = caret::preProcess(data.frame(imp.df$rowmean),method = "range",rangeBounds = c(0,1) )
                scaledScore = stats::predict(trainData.md, data.frame(imp.df$rowmean))
                ## Compute node strength
                el = data.frame(Ratio = imp.df$Ratio,Score = scaledScore[,1])
                el = tidyr::separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
                g = igraph::graph_from_edgelist(as.matrix(el[,2:3]))
                igraph::E(g)$weight = el$Score
                nodeImp = data.frame(Node = names(igraph::strength(g)),Str = igraph::strength(g)) %>%
                  dplyr::arrange(dplyr::desc(Str))

                featureList[[i]] = list(trainData2,testData2)

                if(i<length(sets)){
                  nodeStrength1 = nodeImp[1:sets[i+1],]
                }

              }))

              pb$tick()

            }

            trainData2 = featureList[[which.max(nodePerf$AUC)]][[1]]
            testData2 = featureList[[which.max(nodePerf$AUC)]][[2]]

          },

          {## 3 - All DCV Ranked Boruta Selection


            th = stats::quantile(dcv$rowmean,probs = percent_topDCV)
            unKeyFeatures = dcv$Ratio[dcv$rowmean>th]
            trainData2 = subset(lrs.train,select = unKeyFeatures)
            testData2 = subset(lrs.test,select = unKeyFeatures)

            ### Boruta
            b = Boruta::Boruta(x = trainData2,y = y_train,doTrace = 0,
                               maxRuns = maxBorutaRuns,
                               getImp = Boruta::getImpExtraGini)
            dec = data.frame(Ratio = names(b$finalDecision),Decision = b$finalDecision)
            keep = dec %>%
              dplyr::filter(Decision!="Rejected")
            kr =as.character(keep$Ratio)
            ## select final ratios
            trainData2 = subset(trainData2,select = c(kr))
            testData2 = subset(testData2,select = c(kr))

            ## Estimate Perfomance
            ac = foreach::foreach(i = 1:nruns_rfAUC,.combine = c)%dopar%{
              testh  = ranger::ranger(formula = Status~.,
                                      data = data.frame(Status = y_train,trainData2),
                                      importance = "none",
                                      num.trees = num_trees,probability = T)
              pROC::auc(pROC::multiclass.roc(y_train,testh$predictions))

            }
            mean(ac)
          }
  )


  return(list(train_data = trainData2,test_data = testData2))


}
