#' Build DiCoVar network from logratio coefficients
#'
#' Takes a set of logratio coefficients from a DiCoVar logistic regression model
#' and assembles a logratio network.
#'
#' @param mdl Logistic regression model generated by DiCoVar
#' @param s Lambda penalty parameter to use when retrieving model coefficients
#' @param case_lab Label used for cases when fitting regression model. For
#'    multinomial models, this label MUST match one of the model's class labels.
#' @param ctrl_lab Label used for controls when fitting regression model. For
#'    multinomial models, the 'control' group is everything except members of the
#'    'case' group. Best practice would be to label this "other" or similar.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{keyRats} \tab Data frame with model ratios and coefficients. \cr
#'    \code{graph} \tab An igraph object containing the logratio network. \cr
#'    \code{part_weights} \tab Data frame of weights for individual features. \cr
#'    \tab \cr
#'}
#' @importFrom rlang .data
#' @export
#'
build_logratio_network <- function(mdl, s = "lambda.min",
                                   case_lab = "case", ctrl_lab = "control") {
  ##Retrieve Coefficients
  if(mdl$glmnet.fit[1] == "multnet") {
    # check case label -- must match one of the possible model classes
    if(!(case_lab %in% mdl$glmnet.fit$classnames)) {
      stop(paste0("Case label '", case_lab, "' does not match any class (",
                  paste(mdl$glmnet.fit$classnames, collapse = ", "),
                  ") in model!"))
    }

    # pull coefficients at desired penalty
    cc <- glmnet::coef.glmnet(object = mdl, s = s)
    # result is a list of model coefficients for each possible y value
    # we grab the one matching the case_lab parameter and use that
    cc <- as.matrix(cc[[case_lab]])
    # matrix should have one column of coefficients, which we rename here
    colnames(cc) <- c(case_lab)
    # convert to tibble
    cc <- tibble::as_tibble(cc, rownames = "Feature")
  } else { # binomial model of some kind
    # check control/case labels -- still works OK if they're wrong
    if(mdl$glmnet.fit$classnames[2] != case_lab) {
      warning(paste0("Case label '", case_lab, "' does not match second class (",
                  mdl$glmnet.fit$classnames[2], ") in model!"))
    }
    if(mdl$glmnet.fit$classnames[1] != ctrl_lab) {
      warning(paste0("Control label '", ctrl_lab, "' does not match first class (",
                     mdl$glmnet.fit$classnames[1], ") in model!"))
    }

    # pull coefficients at desired penalty and convert to non-sparse matrix
    cc <- as.matrix(glmnet::coef.glmnet(object = mdl, s = s))
    # matrix should have one column of coefficients, which we rename here
    colnames(cc) <- c(case_lab)
    # convert to tibble
    cc <- tibble::as_tibble(cc, rownames = "Feature")
  }
  # extract ridge regression coefficients
  cc1 <- cc[[case_lab]][-1]
  # Imp = absolute coefficient, raw = original coefficient
  imp.df = data.frame(Ratio = cc$Feature[-1],Imp = abs(as.numeric(cc1)),raw = as.numeric(cc1))
  keyRats = tidyr::separate(imp.df,1,into = c("Num","Denom"),sep = "___",remove = F)

  ## stack ratio for consistent interpretation
  keyRats2 = keyRats
  keyRats2$Num = keyRats$Denom
  keyRats2$Denom= keyRats$Num
  keyRats2$Ratio= paste0(keyRats$Denom,"___",keyRats$Num)
  keyRats2$raw = -keyRats$raw
  # mark original ratio orientation
  keyRats$Orientation = "original"
  keyRats2$Orientation = "inverted"
  keyRats2 = rbind(keyRats,keyRats2)
  ### keep negative edges (more abundance more likely in control outcome)
  keyRats = keyRats2 %>%
    dplyr::filter(raw<0)

  ## Define weight such that:  weight * log(a/b) = weight * log(a) - weight * log(b)
  weights.df = data.frame(Part = keyRats$Num,Coef = keyRats$raw)
  weights.df = rbind(weights.df,data.frame(Part = keyRats$Denom,Coef = -1*keyRats$raw))
  weights.df = weights.df %>%
    dplyr::group_by(.data$Part) %>%
    dplyr::summarise_all(.funs = sum)
  weights.df$col = dplyr::if_else(weights.df$Coef>0, case_lab, ctrl_lab)
  weights.df$col = factor(weights.df$col,levels = c(ctrl_lab, case_lab))

  # build logratio network

  el_= data.frame(keyRats$Num,keyRats$Denom,keyRats$Ratio)
  g = igraph::graph_from_edgelist(as.matrix(el_[,1:2]), directed = T)
  # drop self-loops and combine multiple edges between the same nodes
  g = igraph::simplify(g, remove.loops = TRUE,
                       edge.attr.comb = list(weight = "sum", name = "concat", "ignore"))
  el_act = data.frame(igraph::get.edgelist(g))
  el_act$Ratio = paste0(el_act$X1,"___",el_act$X2)
  # el_act came from el_, came from keyRats
  # imp.df direct from ridge regression
  # some ratios have been inverted
  # so we should actually just use keyRats here
  el_act <- dplyr::left_join(el_act, keyRats, by = "Ratio")
  igraph::E(g)$width <- el_act$Imp *10
  # apply edge weights to graph in Gephi format
  igraph::E(g)$weight <- igraph::E(g)$width

  # set vertex colors and sizes
  vertices = data.frame(Part = igraph::V(g)$name,
                        Label = igraph::V(g)$name)
  vertices = dplyr::left_join(vertices, weights.df, by = "Part")
  # apply sizes to graph
  igraph::V(g)$size <- abs(vertices$Coef)*10 + 1
  # set custom vertex colors based on part coefficients
  # these are 'stolen' from ggsci package to avoid an extra dependency
  v_color = dplyr::if_else(vertices$Coef>0,
                           "#ED000099",
                           "#00468B99")
  igraph::V(g)$color <- v_color
  # embed custom vertex colors in Gephi format
  v_c_rgb <- grDevices::col2rgb(v_color)
  igraph::V(g)$r <- v_c_rgb["red", ]
  igraph::V(g)$g <- v_c_rgb["green", ]
  igraph::V(g)$b <- v_c_rgb["blue", ]

  # make vertex names readable
  igraph::V(g)$name <- igraph::V(g)$name %>%
    stringr::str_sub(start = 2) %>%
    stringr::str_replace(pattern = "(\\..*?)\\.", replacement = "\\1@")
  # copy vertex names so Gephi will show them
  igraph::V(g)$label <- igraph::V(g)$name

  return(list(keyRats = keyRats, graph = g, part_weights = weights.df))
}
