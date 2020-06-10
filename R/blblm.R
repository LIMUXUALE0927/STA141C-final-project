#' @import purrr
#' @import furrr
#' @import future
#' @import stats
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @aliases NULL
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' linear regression with mini bag bootstrap
#'
#' @param formula formula
#' @param data dataset
#' @param m how many times we split the dataset
#' @param B number of bootstrap performed
#' @param num_worker number of workers when we run parallel
#'
#' @return blblm class object
#'
#' @export
#'
blblm <- function(formula, data, m = 10, B = 5000, num_worker = 1) {
  if (is.data.frame(data)) {
    data_list <- split_data(data, m)
    if (num_worker == 1) {
      estimates <- map(data_list,
                              ~ lm_each_subsample(
                                formula = formula,
                                data = .,
                                n = nrow(data),
                                B = B
                              ))
    } else{
      suppressWarnings(plan(multiprocess, workers = num_worker))
      estimates <- future_map(data_list,
                                     ~ lm_each_subsample(
                                       formula = formula,
                                       data = .,
                                       n = nrow(data),
                                       B = B
                                     ))
    }
  } else{
    if (num_worker == 1) {
      nm = data %>% map( ~ {
        df <- read.csv(., )
        nrow(df)
      }) %>% reduce(`+`)
      estimates = data %>% map( ~ {
        df <- read.csv(., )
        lm_each_subsample(
          formula = formula,
          data = df,
          n = nm,
          B = B
        )
      })
    } else{
      suppressWarnings(plan(multiprocess, workers = num_worker))
      nm = data %>% map( ~ {
        df <- read.csv(., )
        nrow(df)
      }) %>% reduce(`+`)
      estimates = data %>% future_map( ~ {
        df <- read.csv(., )
        lm_each_subsample(
          formula = formula,
          data = df,
          n = nm,
          B = B
        )
      })
    }
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' split data into m parts of approximated equal sizes
#'
#' @param data dataset
#' @param m number of subsamples in total
#'
#' @export
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' compute the estimates
#'
#' @param formula formula
#' @param data dataset
#' @param n total sample size
#' @param B number of bootstrap performed
#' @export
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}

#' compute the regression estimates for a blb dataset
#'
#' @param formula formula
#' @param data dataset
#' @param n total sample size
#'
#' @export
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}

#' estimate the regression estimates based on given the number of repetitions
#'
#' @param formula formula
#' @param data dataset
#' @param freqs frequents how many times observations are repeated
#'
#' @export
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' compute the coefficients from fit
#'
#' @param fit fitted model
#'
#' @export
blbcoef <- function(fit) {
  coef(fit)
}

#' compute sigma from fit
#'
#' @param fit fitted model
#'
#' @export
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' print method on blblm
#'
#' @param x blblm object
#' @param ... other parameters
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

##################

#' print method on blblm
#'
#' @param object blblm object
#' @param confidence need confidence interval or not, default is FALSE
#' @param level level of the confidence intervals, default 0.95
#' @param ... other parameters
#'
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' get coefficients from blblm object
#'
#' @param object blblm output
#' @param ... other parameters
#'
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' get confidence interval of blblm object
#'
#' @param object blblm output
#' @param parm parameters
#' @param level confidence level, default 0.95
#' @param ... other parameters
#'
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' prediciton on blblm object
#'
#' @param object blblm object
#' @param new_data dataset we want to predict
#' @param confidence need confidence interval or not
#' @param level confidence level
#' @param ... other parameters
#'
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}

mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
