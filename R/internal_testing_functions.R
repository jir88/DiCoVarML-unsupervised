# Functions here are for internal testing purposes only, not for external use.

#' Toy Dilution Test
#'
#' Generates a toy data set consisting of 6 groups. Groups 1-3 contain 3 random
#' normal variables with means 3, 5, and 7 and standard deviation 1. The groups
#' are scaled by random normal scaling factors with mean factors of 0.5x, 1x,
#' and 2x respectively, with standard deviation 0.01x. Groups 4-6 will contain
#' the same variables, but with means 7, 5, and 3, respectively. Scaling factors
#' will be the same.
#'
#' This pathological data is then run through the unsupervised DiCoVar algorithm.
#'
#' Hypothesis: incorrect permutation techniques will make the truly uninformative
#' variables _look_ important because dilution information will leak into the
#' synthetic data set.
#'
#' Desired result: method identifies two clusters, groups 1-3 and groups 4-6.
#'
#' Incorrect result: method identifies 3+ clusters, separated in part by scaling
#'
#' @param n the number of samples to generate in each group
#' @param replace Should fake data be sampled from real data with replacement?
#' @param dcv_folds Partition the data into folds and average DCV scores from each fold
#' @param dcv_repeats How many times should DCV repeat the partitioning process?
#' @param scale_data Should data be scaled to unit standard deviation before ridge regression?
#' @param alpha amount of LASSO regression to mix into the final ridge regression model
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom magrittr %>%
toy_unsupervised_dilution_test <- function(n = 301,
                                           add_junk = TRUE,
                                           replace = TRUE,
                                           dcv_folds = 1,
                                           dcv_repeats = 1,
                                           scale_data = TRUE,
                                           alpha = 0) {
  gen_samples <- function(n, a, b, c, var_sd, dil, dil_sd) {
    # generate the variables
    df <- cbind(a = stats::rnorm(n, mean = a, sd = var_sd),
                b = stats::rnorm(n, mean = b, sd = var_sd),
                c = stats::rnorm(n, mean = c, sd = var_sd))
    # generate dilution factors for each sample
    dil_factors <- stats::rnorm(n, mean = dil, sd = dil_sd)
    # apply dilution factors
    df <- sweep(df, MARGIN = 1, STATS = dil_factors, FUN = "*")
    return(df)
  }

  a = 19; b = 7; c = 11

  td1 <- rbind(G1 = gen_samples(n, a, b, c, var_sd = 1, dil = 0.5, dil_sd = 0.01),
               G2 = gen_samples(n, a, b, c, var_sd = 1, dil = 1.0, dil_sd = 0.01),
               G3 = gen_samples(n, a, b, c, var_sd = 1, dil = 2.0, dil_sd = 0.01)) %>%
    tibble::as_tibble()
  td2 <- rbind(G4 = gen_samples(n, c, b, a, var_sd = 1, dil = 0.5, dil_sd = 0.01),
               G5 = gen_samples(n, c, b, a, var_sd = 1, dil = 1.0, dil_sd = 0.01),
               G6 = gen_samples(n, c, b, a, var_sd = 1, dil = 2.0, dil_sd = 0.01)) %>%
    tibble::as_tibble()

  toy_data <- dplyr::bind_rows(td1, td2) %>%
    dplyr::mutate(Status = rep(1:6, each = n)) %>%
    dplyr::relocate(Status)

  if(add_junk) {
    toy_data <- toy_data %>%
    # add some junk features
    dplyr::mutate(d = stats::runif(n*6),
                  e = stats::runif(n*6),
                  f = stats::runif(n*6))
  }

  # PCA plot showing worst-case scenario: all groups resolved, PC1 is scaling
  pca_toy <- mixOmics::pca(X = as.matrix(toy_data[, -1]), ncomp = 3)
  plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = toy_data$Status,
                      legend = TRUE)

  udcv_res <- unsupervised.dicovar(X = as.data.frame(toy_data[, -1]),
                                   replace = replace,
                                   dcv_folds = dcv_folds,
                                   dcv_repeats = dcv_repeats,
                                   scale_data = scale_data,
                                   alpha = alpha)

  graphics::barplot(height = udcv_res$dcv_scores$lrs$rowmean,
          names.arg = udcv_res$dcv_scores$lrs$Ratio,
          cex.names = 0.5)

  # PCA plot showing result of permuting before logratio calculation:
  # real data falls in two clusters, while fake data is all over the place
  pca_toy <- mixOmics::pca(X = udcv_res$glm_model$data$x, ncomp = 2)
  # plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = udcv_res$glm_model$data$y,
                      legend = TRUE)

  # grab the important logratio data
  df <- cbind(Status = udcv_res$glm_model$data$y, udcv_res$glm_model$data$x)
  # drop fake data
  df <- df[df$Status == "real", -1]

  # PCA plot showing only real data, filtered logratios
  pca_toy <- mixOmics::pca(X = df, ncomp = 2)
  # plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = toy_data$Status,
                      legend = TRUE)
}

#' Toy Dilution Test
#'
#' Generates a toy data set consisting of 6 groups. Groups 1-3 contain 3 random
#' normal variables with means 3, 5, and 7 and standard deviation 1. The groups
#' are scaled by random normal scaling factors with mean factors of 0.5x, 1x,
#' and 2x respectively, with standard deviation 0.01x. Groups 4-6 will contain
#' the same variables, but with means 7, 5, and 3, respectively. Scaling factors
#' will be the same.
#'
#' This pathological data is then run through the supervised DiCoVar algorithm.
#'
#' Hypothesis:
#'
#' Desired result:
#'
#' Incorrect result:
#'
#' @param n the number of samples to generate in each group
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom magrittr %>%
toy_supervised_dilution_test <- function(n = 301) {
  gen_samples <- function(n, a, b, c, var_sd, dil, dil_sd) {
    # generate the variables
    df <- cbind(a = stats::rnorm(n, mean = a, sd = var_sd),
                b = stats::rnorm(n, mean = b, sd = var_sd),
                c = stats::rnorm(n, mean = c, sd = var_sd))
    # generate dilution factors for each sample
    dil_factors <- stats::rnorm(n, mean = dil, sd = dil_sd)
    # apply dilution factors
    df <- sweep(df, MARGIN = 1, STATS = dil_factors, FUN = "*")
    return(df)
  }

  a = 19; b = 7; c = 11

  td1 <- rbind(G1 = gen_samples(n, a, b, c, var_sd = 1, dil = 0.5, dil_sd = 0.01),
               G2 = gen_samples(n, a, b, c, var_sd = 1, dil = 1.0, dil_sd = 0.01),
               G3 = gen_samples(n, a, b, c, var_sd = 1, dil = 2.0, dil_sd = 0.01)) %>%
    tibble::as_tibble()
  td2 <- rbind(G4 = gen_samples(n, c, b, a, var_sd = 1, dil = 0.5, dil_sd = 0.01),
               G5 = gen_samples(n, c, b, a, var_sd = 1, dil = 1.0, dil_sd = 0.01),
               G6 = gen_samples(n, c, b, a, var_sd = 1, dil = 2.0, dil_sd = 0.01)) %>%
    tibble::as_tibble()

  toy_data <- dplyr::bind_rows(td1, td2) %>%
    dplyr::mutate(Status = stringr::str_c("Group", rep(1:2, each = n*3))) %>%
    dplyr::relocate(Status) %>%
    # add some junk features
    dplyr::mutate(d = stats::runif(n*6),
                  e = stats::runif(n*6),
                  f = stats::runif(n*6))

  # PCA plot showing worst-case scenario: all groups resolved, PC1 is scaling
  pca_toy <- mixOmics::pca(X = as.matrix(toy_data[, -1]), ncomp = 3)
  plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = toy_data$Status,
                      legend = TRUE)

  # randomly select 20% of data as validation set
  dv_folds <- caret::createDataPartition(y = toy_data$Status, times = 1,
                                  p = 0.8, list = FALSE)
  dv_fold_ids <- purrr::rep_along(toy_data$Status, 2)
  # 1-is train data; 2 is test data
  dv_fold_ids[dv_folds] <- 1


  # set up parallel execution
  n_cores <- 5
  doParallel::registerDoParallel(cl = n_cores, cores = n_cores)

  # udcv_res <- unsupervised.dicovar(X = as.data.frame(toy_data[, -1]),
  #                                  alpha = 0.5)
  dcv_res <- DiCoVarML::tune.dicovar(X = as.data.frame(toy_data[, -1]),
                                     Y = toy_data$Status,
                                     sample_ids = seq_along(toy_data$Status),
                                     dv_labels = dv_fold_ids,
                                     ensemble = c("ranger", "glmnet", "xgbTree"),
                                     test.parts = 3:5)


  # extract DCV scores for each feature ratio across inner CV runs
  inner_dcv_scores <- dcv_res$inner_dcv_scores

  ggplot2::ggplot(inner_dcv_scores, ggplot2::aes(x = Ratio, y = rowmean)) +
    ggplot2::geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.5)

  # get ridge model coefficients
  df <- lapply(X = dcv_res$inner_ridge_models, FUN = function(mdl) {
    coefs <- glmnet::coef.glmnet(mdl, s = "lambda.min")
    coefs <- tibble::as_tibble(as.matrix(coefs), rownames = "Feature")
    return((coefs))
    })

  df2 <- dplyr::bind_rows(df) %>%
    dplyr::filter(Feature != "(Intercept)")

  ggplot2::ggplot(df2, ggplot2::aes(x = Feature, y = lambda.min)) +
    ggplot2::geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.1)

  # PCA plot showing result of permuting before logratio calculation:
  # real data falls in two clusters, while fake data is all over the place
  pca_toy <- mixOmics::pca(X = dcv_res$glm_model$data$x, ncomp = 2)
  plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = udcv_res$glm_model$data$y,
                      legend = TRUE)

  # grab the important logratio data
  df <- cbind(Status = udcv_res$glm_model$data$y, udcv_res$glm_model$data$x)
  # drop fake data
  df <- df[df$Status == "real", -1]

  # PCA plot showing only real data, filtered logratios
  pca_toy <- mixOmics::pca(X = df, ncomp = 2)
  plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = toy_data$Status,
                      legend = TRUE)
}

#' Toy Dilution Test
#'
#' Generates a toy data set consisting of 6 groups. Groups 1-3 contain 3 random
#' normal variables with means 3, 5, and 7 and standard deviation 1. The groups
#' are scaled by random normal scaling factors with mean factors of 0.5x, 1x,
#' and 2x respectively, with standard deviation 0.01x. Groups 4-6 will contain
#' the same variables, but with means 7, 5, and 3, respectively. Scaling factors
#' will be the same.
#'
#' This pathological data is then run through the unsupervised DiCoVar algorithm.
#'
#' Hypothesis: incorrect permutation techniques will make the truly uninformative
#' variables _look_ important because dilution information will leak into the
#' synthetic data set.
#'
#' Desired result: method identifies two clusters, groups 1-3 and groups 4-6.
#'
#' Incorrect result: method identifies 3+ clusters, separated in part by scaling
#'
#' @param n the number of samples to generate in each group
#' @param replace Should fake data be sampled from real data with replacement?
#' @param dcv_folds Partition the data into folds and average DCV scores from each fold
#' @param dcv_repeats How many times should DCV repeat the partitioning process?
#' @param scale_data Should data be scaled to unit standard deviation before ridge regression?
#' @param alpha amount of LASSO regression to mix into the final ridge regression model
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom magrittr %>%
iris_unsupervised_DCV_test <- function(X = datasets::iris[-5],
                                       Y = datasets::iris$Species,
                                       mean_dilution = NA,
                                       sd_dilution = NA,
                                       add_junk = c(),
                                       replace = TRUE,
                                           dcv_folds = 1,
                                           dcv_repeats = 1,
                                           scale_data = TRUE,
                                           alpha = 0) {
  # optionally add some junk variables
  if(!is.null(add_junk)) {
    n_obs <- length(Y)
    junk_data <- sapply(X = add_junk, FUN = function(j) {
      return(stats::rnorm(n = n_obs, mean = j, sd = j/10))
    })
    colnames(junk_data) <- stringr::str_c("J", 1:length(add_junk))
    X <- dplyr::bind_cols(X, junk_data)
  }
  # optionally apply a random dilution factor to each data point
  if(!is.na(mean_dilution)) {
    dil_factor <- stats::rnorm(n = length(Y), mean = mean_dilution, sd = sd_dilution)
    X <- sweep(X, MARGIN = 1, STATS = dil_factor, FUN = "/")
  }

  # PCA plot showing worst-case scenario: all groups resolved, PC1 is scaling
  pca_toy <- mixOmics::pca(X = X, ncomp = 3)
  plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = Y,
                      legend = TRUE)

  udcv_res <- unsupervised.dicovar(X = X,
                                   replace = replace,
                                   dcv_folds = dcv_folds,
                                   dcv_repeats = dcv_repeats,
                                   scale_data = scale_data,
                                   alpha = alpha)

  lrs <- dplyr::arrange(udcv_res$dcv_scores$lrs, (rowmean))
  lrs$Ratio <- factor(lrs$Ratio, levels = lrs$Ratio)

  ggplot2::ggplot(lrs, ggplot2::aes(x = rowmean, y = Ratio)) +
    ggplot2::geom_col()

  # PCA plot showing result of permuting before logratio calculation:
  # real data falls in two clusters, while fake data is all over the place
  pca_toy <- mixOmics::pca(X = udcv_res$glm_model$data$x, ncomp = 2)
  # plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = udcv_res$glm_model$data$y,
                      legend = TRUE)

  # grab the important logratio data
  df <- cbind(Status = udcv_res$glm_model$data$y, udcv_res$glm_model$data$x)
  # drop fake data
  df <- df[df$Status == "real", -1]

  # PCA plot showing only real data, filtered logratios
  pca_toy <- mixOmics::pca(X = df, ncomp = 3)
  # plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = Y,
                      legend = TRUE)

  rf_pre <- randomForest::randomForest(x = X, y = as.factor(Y))
  plot(rf_pre)
  rf_pre

  rf_iris <- randomForest::randomForest(x = df, y = as.factor(Y))
  plot(rf_iris)
  rf_iris
}
