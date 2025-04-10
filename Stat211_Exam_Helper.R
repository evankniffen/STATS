################################################################################
# STAT 211 Exam Helper - Comprehensive R Script with Full Diagnostics
# This file implements functions that fully process exam problems dynamically.
# Every function prints a full diagnostic report (if verbose = TRUE) and returns
# a detailed list of computed values.
#
# To use during the exam:
#    source("Stat211_Exam_Helper.R")
# Then, call functions such as:
#    oneSampleTTest(x, mu0, alpha, alternative)
#    processDiscreteJointPMF(joint_pmf, outcome_values_row, condition_on, condition_index)
#    processJointPDF(f_xy, c(x_lower, x_upper), c(y_lower, y_upper))
#    etc.
#
################################################################################

########################################
# SECTION 1: Hypothesis Testing for the Mean
########################################

oneSampleTTest <- function(x, mu0, alpha = 0.05, alternative = "two.sided", verbose = TRUE) {
  n <- length(x)
  xbar <- mean(x)
  s <- sd(x)
  SE <- s / sqrt(n)
  t_obs <- (xbar - mu0) / SE
  
  if (alternative == "two.sided") {
    p_value <- 2 * (1 - pt(abs(t_obs), df = n - 1))
    conf_int <- t.test(x, conf.level = 1 - alpha)$conf.int
  } else if (alternative == "less") {
    p_value <- pt(t_obs, df = n - 1)
    conf_int <- t.test(x, alternative = "less", conf.level = 1 - alpha)$conf.int
  } else if (alternative == "greater") {
    p_value <- 1 - pt(t_obs, df = n - 1)
    conf_int <- t.test(x, alternative = "greater", conf.level = 1 - alpha)$conf.int
  } else {
    stop("Invalid alternative. Choose 'two.sided', 'less', or 'greater'.")
  }
  
  if (verbose) {
    cat("One-Sample t-Test Diagnostics:\n")
    cat("  n =", n, "\n")
    cat("  Sample mean =", round(xbar, 4), "\n")
    cat("  Sample SD =", round(s, 4), "\n")
    cat("  Standard Error =", round(SE, 4), "\n")
    cat("  Test statistic (t) =", round(t_obs, 4), "\n")
    cat("  p-value =", round(p_value, 4), "\n")
    cat("  ", 100*(1 - alpha), "% Confidence Interval =", paste0("[", 
        round(conf_int[1], 4), ", ", round(conf_int[2], 4), "]"), "\n\n")
  }
  
  return(list(n = n, sample_mean = xbar, sample_sd = s,
              SE = SE, t_statistic = t_obs, p_value = p_value, conf_int = conf_int))
}

zTestMean <- function(xbar, mu0, sigma, n, alpha = 0.05, alternative = "two.sided", verbose = TRUE) {
  SE <- sigma / sqrt(n)
  z_obs <- (xbar - mu0) / SE
  if (alternative == "two.sided") {
    p_value <- 2 * (1 - pnorm(abs(z_obs)))
    quantile_used <- qnorm(1 - alpha/2)
  } else if (alternative == "less") {
    p_value <- pnorm(z_obs)
    quantile_used <- qnorm(alpha)
  } else if (alternative == "greater") {
    p_value <- 1 - pnorm(z_obs)
    quantile_used <- qnorm(1 - alpha)
  } else {
    stop("Invalid alternative. Choose 'two.sided', 'less', or 'greater'.")
  }
  
  if (verbose) {
    cat("One-Sample z-Test Diagnostics (sigma known):\n")
    cat("  n =", n, "\n")
    cat("  Sample mean =", round(xbar, 4), "\n")
    cat("  Hypothesized mean =", mu0, "\n")
    cat("  Sigma =", sigma, "\n")
    cat("  Standard Error =", round(SE, 4), "\n")
    cat("  z-statistic =", round(z_obs, 4), "\n")
    cat("  p-value =", round(p_value, 4), "\n")
    cat("  Critical quantile used =", round(quantile_used, 4), "\n\n")
  }
  
  return(list(n = n, sample_mean = xbar, sigma = sigma, SE = SE, z_statistic = z_obs,
              p_value = p_value, critical_quantile = quantile_used))
}

pairedTTest <- function(x1, x2, alpha = 0.05, alternative = "two.sided", verbose = TRUE) {
  if (length(x1) != length(x2)) stop("x1 and x2 must be the same length.")
  diff <- x2 - x1
  n <- length(diff)
  d_bar <- mean(diff)
  s_diff <- sd(diff)
  SE_diff <- s_diff / sqrt(n)
  t_obs <- d_bar / SE_diff
  
  if (alternative == "two.sided") {
    p_value <- 2 * min(pt(t_obs, df = n - 1), 1 - pt(t_obs, df = n - 1))
    conf_int <- t.test(x2, x1, paired = TRUE, conf.level = 1 - alpha)$conf.int
  } else if (alternative == "less") {
    p_value <- pt(t_obs, df = n - 1)
    conf_int <- t.test(x2, x1, paired = TRUE, alternative = "less", conf.level = 1 - alpha)$conf.int
  } else if (alternative == "greater") {
    p_value <- 1 - pt(t_obs, df = n - 1)
    conf_int <- t.test(x2, x1, paired = TRUE, alternative = "greater", conf.level = 1 - alpha)$conf.int
  } else {
    stop("Invalid alternative. Choose 'two.sided', 'less', or 'greater'.")
  }
  
  if (verbose) {
    cat("Paired t-Test Diagnostics:\n")
    cat("  n =", n, "\n")
    cat("  Mean difference =", round(d_bar, 4), "\n")
    cat("  SD of differences =", round(s_diff, 4), "\n")
    cat("  Standard Error =", round(SE_diff, 4), "\n")
    cat("  t-statistic =", round(t_obs, 4), "\n")
    cat("  p-value =", round(p_value, 4), "\n")
    cat("  ", 100*(1 - alpha), "% Confidence Interval =", paste0("[", 
        round(conf_int[1], 4), ", ", round(conf_int[2], 4), "]"), "\n\n")
  }
  
  return(list(n = n, mean_difference = d_bar, sd_difference = s_diff,
              SE_difference = SE_diff, t_statistic = t_obs, p_value = p_value, conf_int = conf_int))
}

########################################
# SECTION 2: Hypothesis Testing for Proportions
########################################

oneSamplePropTest <- function(n, x, p0, alpha = 0.05, alternative = "two.sided", verbose = TRUE) {
  p_hat <- x / n
  SE <- sqrt(p0 * (1 - p0) / n)
  z_obs <- (p_hat - p0) / SE
  
  if (alternative == "two.sided") {
    p_value <- 2 * (1 - pnorm(abs(z_obs)))
    margin <- qnorm(1 - alpha/2) * sqrt(p_hat * (1 - p_hat) / n)
    conf_int <- c(p_hat - margin, p_hat + margin)
  } else if (alternative == "less") {
    p_value <- pnorm(z_obs)
    margin <- qnorm(1 - alpha) * sqrt(p_hat * (1 - p_hat) / n)
    conf_int <- c(NA, p_hat + margin)
  } else if (alternative == "greater") {
    p_value <- 1 - pnorm(z_obs)
    margin <- qnorm(1 - alpha) * sqrt(p_hat * (1 - p_hat) / n)
    conf_int <- c(p_hat - margin, NA)
  } else {
    stop("Invalid alternative. Choose 'two.sided', 'less', or 'greater'.")
  }
  
  if (verbose) {
    cat("One-Sample Proportion Test Diagnostics:\n")
    cat("  Sample size (n) =", n, "\n")
    cat("  Number of successes =", x, "\n")
    cat("  Sample proportion (p̂) =", round(p_hat, 4), "\n")
    cat("  Hypothesized proportion =", p0, "\n")
    cat("  Standard Error =", round(SE, 4), "\n")
    cat("  z-statistic =", round(z_obs, 4), "\n")
    cat("  p-value =", round(p_value, 4), "\n")
    cat("  ", 100*(1 - alpha), "% Confidence Interval =", paste0("[", 
        round(conf_int[1], 4), ", ", round(conf_int[2], 4), "]"), "\n\n")
  }
  
  return(list(n = n, successes = x, p_hat = p_hat, SE = SE, z_statistic = z_obs,
              p_value = p_value, conf_int = conf_int))
}

########################################
# SECTION 3: Confidence Intervals & Sample Size Calculations
########################################

ciMeanKnownSigma <- function(xbar, sigma, n, conf_level = 0.90, verbose = TRUE) {
  alpha <- 1 - conf_level
  z_alpha2 <- qnorm(1 - alpha/2)
  margin <- z_alpha2 * sigma / sqrt(n)
  ci <- c(xbar - margin, xbar + margin)
  
  if (verbose) {
    cat("Confidence Interval Diagnostics (σ known):\n")
    cat("  Sample mean =", round(xbar, 4), "\n")
    cat("  Sigma =", sigma, "\n")
    cat("  n =", n, "\n")
    cat("  Confidence level =", conf_level, "\n")
    cat("  z-critical value =", round(z_alpha2, 4), "\n")
    cat("  Margin of error =", round(margin, 4), "\n")
    cat("  Confidence Interval =", paste0("[", round(ci[1],4), ", ", round(ci[2],4), "]"), "\n\n")
  }
  
  return(ci)
}

sampleSizeMean <- function(sigma, margin_error, conf_level = 0.95, verbose = TRUE) {
  alpha <- 1 - conf_level
  z_alpha2 <- qnorm(1 - alpha/2)
  n_req <- ceiling((z_alpha2 * sigma / margin_error)^2)
  
  if (verbose) {
    cat("Sample Size Calculation for Mean (σ known):\n")
    cat("  Sigma =", sigma, "\n")
    cat("  Desired margin of error =", margin_error, "\n")
    cat("  Confidence level =", conf_level, "\n")
    cat("  z-critical value =", round(z_alpha2, 4), "\n")
    cat("  Required sample size (n) =", n_req, "\n\n")
  }
  
  return(n_req)
}

sampleSizeProportion <- function(margin_error, conf_level = 0.95, verbose = TRUE) {
  alpha <- 1 - conf_level
  z_alpha2 <- qnorm(1 - alpha/2)
  n_req <- ceiling((z_alpha2 * 0.5 / margin_error)^2)
  
  if (verbose) {
    cat("Sample Size Calculation for Proportion (worst-case, p = 0.5):\n")
    cat("  Desired margin of error =", margin_error, "\n")
    cat("  Confidence level =", conf_level, "\n")
    cat("  z-critical value =", round(z_alpha2, 4), "\n")
    cat("  Required sample size (n) =", n_req, "\n\n")
  }
  
  return(n_req)
}

########################################
# SECTION 4: Power Calculations
########################################

powerOneSampleZ <- function(mu0, mu_alt, sigma, n, alpha = 0.05, alternative = "greater", verbose = TRUE) {
  SE <- sigma / sqrt(n)
  if (alternative == "greater") {
    z_alpha <- qnorm(1 - alpha)
    rejection_boundary <- mu0 + z_alpha * SE
    power <- 1 - pnorm(rejection_boundary, mean = mu_alt, sd = SE)
  } else if (alternative == "less") {
    z_alpha <- qnorm(alpha)
    rejection_boundary <- mu0 + z_alpha * SE
    power <- pnorm(rejection_boundary, mean = mu_alt, sd = SE)
  } else if (alternative == "two.sided") {
    z_alpha2 <- qnorm(1 - alpha/2)
    power <- 1 - (pnorm(mu0 + z_alpha2 * SE, mean = mu_alt, sd = SE) -
                  pnorm(mu0 - z_alpha2 * SE, mean = mu_alt, sd = SE))
  } else {
    stop("Invalid alternative. Choose 'greater', 'less', or 'two.sided'.")
  }
  
  if (verbose) {
    cat("Power Calculation Diagnostics:\n")
    cat("  H0: mu =", mu0, "\n")
    cat("  Alternative mean =", mu_alt, "\n")
    cat("  Sigma =", sigma, ", n =", n, "\n")
    cat("  Standard Error =", round(SE, 4), "\n")
    if (alternative %in% c("greater", "less")) {
      cat("  Critical (z) value =", round(qnorm(ifelse(alternative == "greater", 1 - alpha, alpha)), 4), "\n")
      cat("  Rejection boundary =", round(rejection_boundary, 4), "\n")
    } else {
      cat("  Critical t-quantile (two-sided) =", round(qnorm(1 - alpha/2), 4), "\n")
    }
    cat("  Estimated power =", round(power, 4), "\n\n")
  }
  
  return(power)
}

########################################
# SECTION 5: Sampling Distribution Simulation (CLT)
########################################

simulateSamplingDistribution <- function(population_fun, n, num_samples = 1000, verbose = TRUE, ...) {
  set.seed(123)
  sample_means <- replicate(num_samples, mean(do.call(population_fun, c(list(n), list(...)))))
  if (verbose) {
    cat("Sampling Distribution Simulation (CLT):\n")
    cat("  Function: ", deparse(substitute(population_fun)), "\n")
    cat("  Sample size (n) =", n, "\n")
    cat("  Number of samples =", num_samples, "\n")
    cat("  Mean of sample means =", round(mean(sample_means), 4), "\n")
    cat("  SD of sample means =", round(sd(sample_means), 4), "\n")
    hist(sample_means, breaks = 30, probability = TRUE,
         main = "Sampling Distribution of the Mean",
         xlab = "Sample Mean")
    curve(dnorm(x, mean = mean(sample_means), sd = sd(sample_means)),
          add = TRUE, col = "red", lwd = 2)
    cat("\n")
  }
  return(sample_means)
}

########################################
# SECTION 6: Joint, Marginal, and Conditional Distributions (Continuous)
########################################

processJointPDF <- function(f_xy, x_bounds, y_bounds, compute_marginals = TRUE, 
                            compute_expectations = TRUE, region = NULL, verbose = TRUE) {
  # Total probability mass over the region (should be 1 for a proper pdf)
  total_mass <- integrate(function(y) {
    sapply(y, function(z) integrate(function(x) f_xy(x, z), lower = x_bounds[1], upper = x_bounds[2])$value)
  }, lower = y_bounds[1], upper = y_bounds[2])$value
  
  result <- list(total_mass = total_mass)
  
  if (compute_expectations) {
    # Compute E[X]:
    E_x <- integrate(function(y) {
      sapply(y, function(z) integrate(function(x) x * f_xy(x, z), lower = x_bounds[1], upper = x_bounds[2])$value)
    }, lower = y_bounds[1], upper = y_bounds[2])$value
    # Compute E[Y]:
    E_y <- integrate(function(x) {
      sapply(x, function(z) integrate(function(y) y * f_xy(z, y), lower = y_bounds[1], upper = y_bounds[2])$value)
    }, lower = x_bounds[1], upper = x_bounds[2])$value
    result$E_x <- E_x
    result$E_y <- E_y
  }
  
  if (compute_marginals) {
    # Marginal fX as a function:
    fX <- function(x) {
      if (x < x_bounds[1] || x > x_bounds[2]) return(0)
      integrate(function(y) f_xy(x, y), lower = y_bounds[1], upper = y_bounds[2])$value
    }
    # Marginal fY as a function:
    fY <- function(y) {
      if (y < y_bounds[1] || y > y_bounds[2]) return(0)
      integrate(function(x) f_xy(x, y), lower = x_bounds[1], upper = x_bounds[2])$value
    }
    result$fX <- fX
    result$fY <- fY
  }
  
  if (!is.null(region)) {
    # region is a list: region$x = c(ax, bx), region$y = c(cy, dy)
    region_prob <- integrate(function(y) {
      sapply(y, function(z) integrate(function(x) f_xy(x, z), lower = region$x[1], upper = region$x[2])$value)
    }, lower = region$y[1], upper = region$y[2])$value
    result$region_prob <- region_prob
  }
  
  if (verbose) {
    cat("Joint PDF Diagnostics:\n")
    cat("  x bounds: [", x_bounds[1], ",", x_bounds[2], "]\n")
    cat("  y bounds: [", y_bounds[1], ",", y_bounds[2], "]\n")
    cat("  Total mass =", round(total_mass, 4), "\n")
    if (compute_expectations) {
      cat("  E[X] =", round(result$E_x, 4), "\n")
      cat("  E[Y] =", round(result$E_y, 4), "\n")
    }
    if (compute_marginals) {
      cat("  Marginal functions for X and Y are available in the output list.\n")
    }
    if (!is.null(region)) {
      cat("  Probability over specified region =", round(result$region_prob, 4), "\n")
    }
    cat("\n")
  }
  
  return(result)
}

########################################
# SECTION 7: Covariance and Correlation
########################################

computeCovCor <- function(X, Y, verbose = TRUE) {
  cov_value <- cov(X, Y)
  cor_value <- cor(X, Y)
  if (verbose) {
    cat("Covariance and Correlation Diagnostics:\n")
    cat("  Covariance =", round(cov_value, 4), "\n")
    cat("  Correlation =", round(cor_value, 4), "\n\n")
  }
  return(list(covariance = cov_value, correlation = cor_value))
}

covLinear <- function(a, c, cov_XY, verbose = TRUE) {
  result <- a * c * cov_XY
  if (verbose) {
    cat("Covariance of Linear Combinations:\n")
    cat("  Cov(", a, "*X, ", c, "*Y) = ", result, "\n\n", sep = "")
  }
  return(result)
}

########################################
# SECTION 8: Discrete Joint and Conditional Distributions
########################################

processDiscreteJointPMF <- function(joint_pmf, outcome_values_row = NULL, outcome_values_col = NULL,
                                      condition_on = NULL, condition_index = NULL, verbose = TRUE) {
  # Compute marginal distributions
  marg_row <- rowSums(joint_pmf)
  marg_col <- colSums(joint_pmf)
  
  result <- list(marginal_rows = marg_row, marginal_columns = marg_col)
  
  if (!is.null(condition_on) && !is.null(condition_index)) {
    if (condition_on == "col") {
      cond_pmf <- joint_pmf[, condition_index] / sum(joint_pmf[, condition_index])
      result$conditional_PMF <- cond_pmf
      if (!is.null(outcome_values_row)) {
        exp_val <- sum(outcome_values_row * cond_pmf)
        result$conditional_expected_value <- exp_val
      }
    } else if (condition_on == "row") {
      cond_pmf <- joint_pmf[condition_index, ] / sum(joint_pmf[condition_index, ])
      result$conditional_PMF <- cond_pmf
      if (!is.null(outcome_values_col)) {
        exp_val <- sum(outcome_values_col * cond_pmf)
        result$conditional_expected_value <- exp_val
      }
    } else {
      stop("condition_on must be 'row' or 'col'.")
    }
  }
  
  if (verbose) {
    cat("Discrete Joint PMF Diagnostics:\n")
    cat("  Joint PMF matrix:\n")
    print(joint_pmf)
    cat("  Marginal PMF (rows):\n")
    print(marg_row)
    cat("  Marginal PMF (columns):\n")
    print(marg_col)
    if (!is.null(condition_on) && !is.null(condition_index)) {
      cat("  Conditional PMF (given", condition_on, "=", condition_index, "):\n")
      print(result$conditional_PMF)
      if (!is.null(result$conditional_expected_value))
        cat("  Conditional Expected Value =", round(result$conditional_expected_value, 4), "\n")
    }
    cat("\n")
  }
  
  return(result)
}

########################################
# SECTION 9: Helper for Normal Sample Mean Probabilities
########################################

sampleMeanProbability <- function(mu, sigma, n, target, type = c("equal", "less"), verbose = TRUE) {
  type <- match.arg(type)
  if (type == "equal") {
    prob <- 0  # For continuous distributions, P(exactly equals target) = 0.
  } else if (type == "less") {
    prob <- pnorm(target, mean = mu, sd = sigma / sqrt(n))
  }
  if (verbose) {
    cat("Sample Mean Probability Diagnostics:\n")
    cat("  Population mean =", mu, "; Sigma =", sigma, "; Sample size =", n, "\n")
    cat("  Target value =", target, "; Type =", type, "\n")
    cat("  Computed probability =", round(prob, 4), "\n\n")
  }
  return(prob)
}

################################################################################
# END OF STAT211 EXAM HELPER FUNCTIONS
#
# All calculations are performed inside each function and a full diagnostic
# output is printed. You only need to call these functions with dynamic inputs.
#
# Example usage:
# 
# --- For Problem 1 (Discrete Joint PMF):
# joint_pmf <- matrix(c(
#   0.02, 0.08, 0.09,
#   0.10, 0.07, 0.04,
#   0.10, 0.09, 0.01,
#   0.07, 0.01, 0.11,
#   0.10, 0.06, 0.05
# ), nrow = 5, byrow = TRUE)
# rownames(joint_pmf) <- 1:5
# colnames(joint_pmf) <- 1:3
# result1 <- processDiscreteJointPMF(joint_pmf, outcome_values_row = 1:5,
#                                      condition_on = "col", condition_index = "1")
#
# --- For Problem 2 (Continuous Joint PDF):
# f_xy <- function(x, y) { (1/25)*(21*x^2 + 27*y^2 + 36*x*y) }
# result2 <- processJointPDF(f_xy, x_bounds = c(0,1), y_bounds = c(0,1))
#
# --- For a hypothesis test, call oneSampleTTest, pairedTTest, or oneSamplePropTest.
#
# --- For CIs and sample sizes, call ciMeanKnownSigma, sampleSizeMean, sampleSizeProportion.
#
# --- For power, call powerOneSampleZ.
#
# --- For simulating sampling distributions, call simulateSamplingDistribution.
#
################################################################################
