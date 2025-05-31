# ---- SuperLearner R6 Class ----

#' SuperLearner
#'
#' Ensemble learning class that supports fitting multiple base learners and
#' a meta-learner using cross-validation.
#'
#' @field learner_library Library of base learners.
#' @field X Dataframe or matrix of predictors.
#' @field Y Vector of outcomes.
#' @field learners Base learners used for the SuperLearner.
#' @field fits Fitted base learners.
#' @field meta_learner The meta_learner learned using cross-validation.
#' @export
SuperLearner <- R6::R6Class("Superlearner",
  public = list(
    learner_library = list(
      glmnet = list(
        name = "glmnet",
        fit = function(X, Y) glmnet::cv.glmnet(x = X, y = Y, family = binomial()),
        predict = function(model, X) predict(model, newx = X, s = "lambda.min", type = "response"),
        x_type = "matrix"
      ),
      randomForest = list(
        name = "randomForest",
        fit = function(X, Y) randomForest::randomForest(as.factor(Y) ~ ., data = cbind(X, Y), ntree = 100),
        predict = function(model, X) predict(model, X, type = "prob")[, 2],
        x_type = "data.frame"
      ),
      earth = list(
        name = "earth",
        fit = function(X, Y) earth::earth(Y ~ ., data = cbind(X, Y), glm = list(family = binomial)),
        predict = function(model, X) predict(model, X, type = "response"),
        x_type = "data.frame"
      ),
      glm = list(
        name = "glm",
        fit = function(X, Y) glm(Y ~ ., data = cbind(X, Y), family = binomial()),
        predict = function(model, X) predict(model, X, type = "response"),
        x_type = "data.frame"
      )
    ),
    X = NULL,
    Y = NULL,
    learners = list(),
    fits = list(),
    meta_learner = NULL,


    #' @description Initializes the SuperLearner.
    #'
    #' @param X Dataframe or matrix of predictors.
    #' @param Y Vector of outcomes.
    #' @param learner_library library of base learners to use.
    initialize = function(X, Y, learner_library = c("glmnet", "randomForest", "earth")) {
      self$X <- X
      self$Y <- Y

      for (name in learner_library) {
        self$learners[[name]] <- list(
          name = name,
          fit = self$learner_library[[name]]$fit,
          predict = self$learner_library[[name]]$predict,
          x_type = self$learner_library[[name]]$x_type
        )
      }
    },

    #' @description Manually adds base learners to SuperLearner object.
    #'
    #' @param name Name of the base learner.
    #' @param fit_fun The fit function of the base learner.
    #' @param predict_fun The predict function of the base learner.
    #' @param x_type Data type for the predictors, default is data.frame.
    add_learner = function(name, fit_fun, predict_fun, x_type = "data.frame") {
      self$learners[[name]] <- list(
        name = name,
        fit = fit_fun,
        predict = predict_fun,
        x_type = x_type
      )
    },

    #' @description Fits all base learners using their corresponding fit function.
    fit_learners = function() {
      for (name in names(self$learners)) {
        learner <- self$learners[[name]]

        ## Check X type and convert accordingly
        if (learner$x_type == "matrix") {
          X_to_use <- as.matrix(self$X)
        } else {
          X_to_use <- as.data.frame(self$X)
        }

        ## Fit learner and save the fit
        model <- learner$fit(X_to_use, self$Y)
        self$fits[[name]] <- model
      }
    },

    #' @description Perform cross-validation using all base learners.
    #'
    #' @param V Number of folds for the cross-validation, default is 10.
    #' @return A matrix with n rows and one column per base learner. Each entry
    #' corresponds to the prediction made by a base learner on the held-out sample.
    cross_validate = function(V = 10) {
      n <- nrow(self$X)
      folds <- sample(rep(1:V, length.out = n))
      results_matrix <- matrix(NA, n, length(self$learners))
      colnames(results_matrix) <- names(self$learners)

      ## Progress bar
      pb <- progress::progress_bar$new(
        format = "  Cross-validating [:bar] :percent eta: :eta",
        total = V * length(self$learners), clear = FALSE, width = 100
      )

      ## Cross-validation loop
      for (v in 1:V) {
        idx_train <- which(folds != v)
        idx_valid <- which(folds == v)

        X_train <- self$X[idx_train, ]
        Y_train <- self$Y[idx_train]

        X_valid <- self$X[idx_valid, ]

        for (name in names(self$learners)) {
          learner <- self$learners[[name]]

          ## Check X type and convert accordingly
          if (learner$x_type == "matrix") {
            X_train <- as.matrix(X_train)
            X_valid <- as.matrix(X_valid)
          } else {
            X_train <- as.data.frame(X_train)
            X_valid <- as.data.frame(X_valid)
          }

          model <- learner$fit(X_train, Y_train)
          preds <- learner$predict(model, X_valid)
          results_matrix[idx_valid, name] <- preds
          pb$tick()
        }
      }
      return(results_matrix)
    },

    #' @description Fits the meta learner using glm (minimizes CV-MSE).
    #'
    #' @param results_matrix The resulting matrix from the cross_validate function.
    fit_meta_learner = function(results_matrix) {
      meta_data <- as.data.frame(results_matrix)
      meta_data$Y <- self$Y
      self$meta_learner <- glm(Y ~ ., data = meta_data)
    },


    #' @description Computes prediction on predictors using the meta-learner.
    #'
    #' @param X A dataframe or matrix of predictors.
    #' @return A vector of predictions
    predict = function(X) {
      preds <- list()
      for (name in names(self$learners)) {
        learner <- self$learners[[name]]
        model <- self$fits[[name]]


        ## Check X type and convert accordingly
        if (learner$x_type == "matrix") {
          X_to_use <- as.matrix(X)
        } else {
          X_to_use <- as.data.frame(X)
        }

        preds[[name]] <- learner$predict(model, X_to_use)
      }

      meta_features <- as.data.frame(preds)
      colnames(meta_features) <- names(self$learners)

      final_pred <- predict(self$meta_learner, meta_features, "response")

      return(final_pred)
    },

    #' @description Computes the full stack of functions to learn the meta-learner and prints
    #' the glm output for the meta-learner.
    #'
    #' @param V Number of folds for cross-validation, default is 10
    learn = function(V = 10) {
      self$fit_learners()
      results <- self$cross_validate(V)
      self$fit_meta_learner(results)
      print(self$meta_learner)
    }
  )
)
