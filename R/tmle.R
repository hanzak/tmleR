# ---- TMLE ----

#' TMLE
#'
#' TMLE class that computes ATE and confidence interval of the ATE assuming
#'
#' @field W_A Dataframe or matrix of predictors (includes covariates W and treatments A).
#' @field Y Vector of outcomes.
#' @field outcome_superLearner SuperLearner used for expected outcome predictions.
#' @field propensity_superLearner SuperLearner used for propensity score predictions.
#' @field treatment_name Name of the treatment column in W_A.
#' @field full_tibble Tibble containing all relevant measures.
#' @field fit glm fit of the estimating equation.
#' @field eps Coefficient of clever covariate of the fit.
#' @field tmle_ate ATE value estimation.
#' @export
TMLE <- R6::R6Class("TMLE",
  public = list(
    W_A = NULL,
    Y = NULL,
    outcome_superLearner = NULL,
    propensity_superLearner = NULL,
    treatment_name = NULL,
    full_tibble = NULL,
    fit = NULL,
    eps = NULL,
    tmle_ate = NULL,

    #' @description Initilizes a TMLE object.
    #'
    #' @param W_A Dataframe or matrix of predictors (includes covariates W and treatments A).
    #' @param Y Vector of outcomes.
    #' @param outcome_superLearner SuperLearner used for expected outcome predictions.
    #' @param propensity_superLearner SuperLearner used for propensity score predictions.
    #' @param treatment_name Name of the treatment column in W_A.
    initialize = function(W_A, Y, outcome_superLearner, propensity_superLearner, treatment_name = "A") {
      self$W_A <- W_A
      self$Y <- Y
      self$outcome_superLearner <- outcome_superLearner
      self$propensity_superLearner <- propensity_superLearner
      self$treatment_name <- treatment_name
    },

    #' @description Computes the ATE and confidence interval of the ATE.
    #'
    #' @return Named list containing estimated_ATE as the estimation of the ATE,
    #' ci_low as the lower bound of the confidence interval and ci_high as the
    #' upper bound of the confidence interval.
    compute_ATE = function() {
      ## Initialize only treated vs only control dataframes
      W_A1 <- self$W_A |> dplyr::mutate(!!self$treatment_name := 1)
      W_A0 <- self$W_A |> dplyr::mutate(!!self$treatment_name := 0)

      ## Predict expected outcome given covariates for observed treatment (Q_A)
      ## Predict expected outcome given covariates for all treated (Q_1)
      ## Predict expected outcome given covariates for all control (Q_0)
      Q_A <- as.vector(self$outcome_superLearner$predict(self$W_A))
      Q_1 <- as.vector(self$outcome_superLearner$predict(W_A1))
      Q_0 <- as.vector(self$outcome_superLearner$predict(W_A0))

      Q_A <- as.vector(self$outcome_superLearner$predict(self$W_A))
      Q_1 <- as.vector(self$outcome_superLearner$predict(W_A1))
      Q_0 <- as.vector(self$outcome_superLearner$predict(W_A0))

      ## Bound predictions between 0 and 1
      Q_A <- pmax(1e-4, pmin(1 - 1e-4, Q_A))
      Q_1 <- pmax(1e-4, pmin(1 - 1e-4, Q_1))
      Q_0 <- pmax(1e-4, pmin(1 - 1e-4, Q_0))

      ## Predict propensity score and get clever covariate terms
      W <- self$W_A |> dplyr::select(!A)
      A <- self$W_A |> dplyr::pull(A)
      g_W <- as.vector(self$propensity_superLearner$predict(W))
      H_1 <- 1 / g_W
      H_0 <- -1 / (1 - g_W)

      ## Get clever covariate based on observed treatment
      binded_cols <- cbind(A, H_1, H_0)
      colnames(binded_cols) <- c("A", "H_1", "H_0")
      H_A <- dplyr::case_when(A == 1 ~ H_1, A == 0 ~ H_0)

      ## Construct full_tibble
      self$full_tibble <- tibble(Y = self$Y) |>
        dplyr::mutate(
          W = W,
          A = A,
          Q_A = Q_A,
          Q_1 = Q_1,
          Q_0 = Q_0,
          g_W = g_W,
          H_A = H_A,
          H_1 = H_1,
          H_0 = H_0,
        )

      ## Estimate the fluctuation parameter by solving the estimating
      ## function for the EIF of the ATE estimand.
      self$fit <- glm(self$Y ~ -1 + offset(qlogis(Q_A)) + H_A, family = binomial())

      ## Save the fluctuation parameter (eps) coefficient from the logistic regression
      self$eps <- coef(self$fit)

      ## Update the initial expected outcomes
      Q_A_update <- plogis(qlogis(Q_A) + self$eps * H_A)
      Q_1_update <- plogis(qlogis(Q_1) + self$eps * H_1)
      Q_0_update <- plogis(qlogis(Q_0) + self$eps * H_0)

      ## Compute estimand of interest
      self$tmle_ate <- mean(Q_1_update - Q_0_update)

      empirical_influence_function <- (self$Y - Q_A_update) * H_A + Q_1_update - Q_0_update - self$tmle_ate
      tmle_se <- sqrt(var(empirical_influence_function) / nrow(self$W_A))
      ci_low <- self$tmle_ate - qnorm(0.975) * tmle_se
      ci_high <- self$tmle_ate + qnorm(0.975) * tmle_se

      return(list(estimated_ATE = self$tmle_ate, ci_low = ci_low, ci_high = ci_high))
    }
  )
)
