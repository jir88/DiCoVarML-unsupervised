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
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom magrittr %>%
toy_dilution_test <- function(n = 301) {
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

  # PCA plot showing worst-case scenario: all groups resolved, PC1 is scaling
  pca_toy <- mixOmics::pca(X = as.matrix(toy_data[, -1]), ncomp = 3)
  plot(pca_toy)
  mixOmics::plotIndiv(pca_toy, ind.names = FALSE, group = toy_data$Status,
                      legend = TRUE)

  udcv_res <- unsupervised.dicovar(X = as.data.frame(toy_data[, -1]),
                                   alpha = 0.5)

  # PCA plot showing result of permuting before logratio calculation:
  # real data falls in two clusters, while fake data is all over the place
  pca_toy <- mixOmics::pca(X = udcv_res$glm_model$data$x, ncomp = 2)
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
