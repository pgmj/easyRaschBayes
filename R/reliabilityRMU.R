### The code below is borrowed from https://github.com/giac01/gbtoolbox/blob/main/R/reliability.R
### and slightly modified to be used with `RIreliabilityRMU()`.
### version used: https://github.com/giac01/gbtoolbox/commit/0e5c5d1547abc61c44f7c72c74c36d829241a33a

#' Estimate reliability (Relative Measurement Uncertainty) from Bayesian measurement models
#'
#' This function measures reliability using posterior draws from a fitted Bayesian model.
#'
#' To use this function, you will need to provide a matrix (input_draws) that contains the posterior draws for the parameter you wish to calculate reliability. The function assumes that rows of input_draws represent subjects and columns represent posterior draws.
#'
#' For an example of how to apply this function to calculate mean score reliability using brms, see \href{https://www.bignardi.co.uk/8_bayes_reliability/tutorial_rmu_sum_score_reliability.html}{this tutorial}.
#'
#' For an example of how to apply this function to go/go-no task data using brms, see \href{https://www.bignardi.co.uk/8_bayes_reliability/tutorial_calculating_rmu_gonogo.html}{this tutorial}.
#'
#' @param input_draws A matrix or data frame of posterior draws. Rows represent subjects and columns represent draws.
#' @param verbose Logical. Print detailed information about the input data. Default is TRUE.
#' @param level Numeric. Credibility level for the highest density continuous interval. Default is 0.95.
#'
#' @return A list containing:
#' \itemize{
#'   \item hdci: A data frame with a point-estimate (posterior mean) and highest density continuous interval for reliability, calculated using the ggdist::mean_hdci function
#'   \item reliability_posterior_draws: A numeric vector of posterior draws for reliability, of length K/2 (K = number of columns/draws in your input_draws matrix)
#' }
#'
#' @references
#' Bignardi, G., Kievit, R., & Bürkner, P. C. (2025). A general method for estimating reliability using Bayesian Measurement Uncertainty. PsyArXiv. \href{https://osf.io/preprints/psyarxiv/h54k8}{doi:10.31234/osf.io/h54k8}
#'
#' @examples
#' \dontrun{
#' # See https://www.bignardi.co.uk/8_bayes_reliability/tutorial_rmu_sum_score_reliability.html for more details on this example
#'
#' # Simulate data
#'
#' set.seed(1)
#' N                   = 5000 # number of subjects (mice)
#' J                   = 3    # number of measurements per subject
#' true_score_variance = 1
#' error_variance      = 10
#'
#' df = expand.grid(j = 1:J, mouse = 1:N)
#'
#' true_scores       = rnorm(N, mean = 10, sd = sqrt(true_score_variance))
#' measurement_error = rnorm(N*J, mean = 0, sd = sqrt(error_variance))
#'
#' df$measurement = true_scores[df$mouse] + measurement_error
#'
#' df_average_lengths = df %>%
#'   group_by(mouse) %>%
#'   summarise(average_measurement = mean(measurement))
#'
#' # Reliability should equal this:
#'
#' true_score_variance/(true_score_variance+error_variance/J)
#'
#' # Approximately the same as:
#'
#' cor(df_average_lengths$average_measurement, true_scores)^2
#'
#' # Fit model and calculate RMU
#'
#' brms_model = brm(
#'   measurement ~ 1 + (1 | mouse),
#'   data    = df
#' )
#'
#' # Extract posterior draws from brms model
#'
#' posterior_draws = brms_model %>%
#'   as_draws_df() %>%
#'   select(starts_with("r_mouse")) %>%
#'   t()
#'
#' # Calculate RMU
#'
#' reliability(posterior_draws)$hdci
#' }
#'
#'
#' @export
RMUreliability <- function(
    input_draws,
    verbose = FALSE,
    level   = .95
) {
  if (!is.numeric(level) || level <= 0 || level >= 1) stop("level must be a numeric value between 0 and 1")

  input_draws <- as.matrix(input_draws)

  # Check for columns with zero SD
  sds <- apply(input_draws, 2, stats::sd, na.rm = TRUE)
  zero_sd_cols <- which(sds == 0)
  if (length(zero_sd_cols) > 0) {
    warning(sprintf("Found %d column(s) with zero standard deviation (columns: %s)",
                    length(zero_sd_cols),
                    paste(zero_sd_cols, collapse = ", ")))
  }

  # Check for NAs in input
  na_count <- base::sum(base::is.na(input_draws))
  if (na_count > 0) {
    warning(sprintf("Found %d NA value(s) in input_draws", na_count))
  }

  if (verbose) {
    base::print(paste0("Number of subjects: ", nrow(input_draws)))
    base::print(paste0("Number of posterior draws: ", ncol(input_draws)))
  }

  col_select <- base::sample(1:ncol(input_draws), replace = FALSE)
  input_draws_1 <- input_draws[, col_select[1:floor(length(col_select) / 2)]]
  input_draws_2 <- input_draws[, col_select[(floor(length(col_select) / 2) + 1):length(col_select)]]

  # Calculate correlations and handle NAs
  reliability_posterior_draws <- sapply(1:ncol(input_draws_1), function(i) {
    x <- input_draws_1[, i]
    y <- input_draws_2[, i]

    # Return 0 if either column has zero variance
    if (stats::var(x, na.rm = TRUE) == 0 || stats::var(y, na.rm = TRUE) == 0) {
      return(0)
    }

    # Calculate correlation and handle NA
    cor_val <- stats::cor(x, y, method = "pearson")
    return(cor_val)
  })

  hdci <- ggdist::mean_hdci(reliability_posterior_draws, .width = level)

  colnames(hdci)[1] = "rmu_estimate"
  colnames(hdci)[2] = "hdci_lowerbound"
  colnames(hdci)[3] = "hdci_upperbound"

  return(hdci)
}
