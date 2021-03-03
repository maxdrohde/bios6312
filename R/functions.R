##################################################
#' Gets the robust variance from a linear model
#'
#' Inputs a linear model from the lm function and prints / returns the
#' robust variance and robust standard errors
#'
#'
#' @param model An lm object
#'
#' @return A list containing the robust variance and robust Std Err respectively
#' @export
#' @examples
#' # get_robust_variance(my_linear_model)
get_robust_variance <- function(model){
  robust.var <- sandwich::vcovHC(model, type = "HC1")
  stderr <- sqrt(diag(robust.var))

  cat("Robust Variance \n")
  print(robust.var)

  cat("\n\n")

  cat("Robust Standard Error \n")
  print(stderr)

  return(invisible(list(`Robust Variance` = robust.var,
                        `Robust Standard Error` = stderr)))
}

##################################################
#' Gets the studentized residuals from a linear model
#'
#' Inputs a linear model from the lm function and returns the
#' studentized residuals as a vector
#'
#'
#' @param model An lm object
#'
#'
#' @return Numeric vector of studentized residuals
#' @export
#' @examples
#' # get_studentized_residuals(my_linear_model)
get_studentized_residuals <- function(model){

  ## Residuals from regression model
  resid <- model$residuals
  ## Estimate error variance
  sigma.hat <- stats::sd(model$residuals)
  ## Create hat matrix
  dsn.X <- cbind(1, model$model$age)
  H <- dsn.X %*% solve(t(dsn.X) %*% dsn.X) %*% t(dsn.X)
  ## Diagonal entries (leverage)
  lvg <- diag(H)
  ## Create studentized residuals
  st.resid <- resid/(sigma.hat * sqrt(1 - lvg))
  return(st.resid)
}

##################################################
#' Creates a normal quantile-quantile plot for the studentized residuals
#'
#' Inputs a linear model from the lm function and returns a
#' normal quantile-quantile plot for the studentized residuals
#'
#' @param model An lm object
#'
#'
#' @return A normal quantile-quantile plot for the studentized residuals
#' @export
#' @examples
#' # qqnorm_studentized_residuals(my_linear_model)
qqnorm_studentized_residuals <- function(model){

  st.resid <- get_studentized_residuals(model)

  stats::qqnorm(st.resid, frame = FALSE,
         cex = 0.8, pch = 20, col = "gray40",
         xlim = c(-4,4), ylim = c(-4, 4))
  stats::qqline(st.resid, lwd = 1.5)
}

##################################################
#' Creates a plot of studentized residuals vs fitted values
#'
#' Inputs a linear model from the lm function and creates a plot of
#' studentized residuals vs fitted values
#'
#'
#' @param model An lm object
#'
#'
#' @return A plot of studentized residuals vs fitted values
#' @export
#' @examples
#' # studentized_residuals_vs_fitted_plot(my_linear_model)
studentized_residuals_vs_fitted_plot <- function(model){

  st.resid <- get_studentized_residuals(model)

  stats::scatter.smooth(model$fitted.values, st.resid,
                 xlab = "Predicted values",
                 ylab = "Studentized Residuals",
                 cex = 0.8, pch = 20, col = "gray40",
                 lpars = list(lwd = 3, col = "blue"),
                 frame.plot = FALSE)
  graphics::abline(0,0, lty = 3, lwd = 2, col = "red")

}

##################################################
#' Creates a plot of studentized residuals vs predictors
#'
#' Inputs a linear model from the lm function and creates a plot of
#' studentized residuals vs predictors
#'
#'
#' @param model An lm object
#'
#'
#' @return A plot of studentized residuals vs predictors
#' @export
#' @examples
#' # studentized_residuals_vs_predictors_plot(my_linear_model)
studentized_residuals_vs_predictors_plot <- function(model){

  st.resid <- get_studentized_residuals(model)

  stats::scatter.smooth(model$model[[2]], st.resid,
                 xlab = names(model$model[2]),
                 ylab = "Studentized Residuals",
                 cex = 0.8, pch = 20, col = "gray40",
                 lpars = list(lwd = 3, col = "blue"),
                 frame.plot = FALSE)
  graphics::abline(0,0, lty = 3, lwd = 2, col = "red")
}

##################################################
#' Test the hypothesis that a set of coefficients are equal to zero
#'
#' Inputs a linear model from the lm function and set of indices for
#' coefficients. Tests the hypothesis the specified set of coefficients
#' are equal to zero. NOTE: beta0 = 1, beta1 = 2, and so on.
#'
#' @param par A vector specifying the indices of the coefficients to be tested
#' @param model An lm object
#' @param type Indicates whether an F-test or Wald test should be used
#'
#' @return A vector containing the test statistic and p-value
#' @export
#' @examples
#' # Tests if beta1, beta4, beta5, and beta7 == 0 using an F-test
#' # testparm(par = c(2,5,6,8),
#' #          model=my_linear_model,
#' #          type="F)
testparm <- function(par, model, type = "F") {

  coefs <- model$coefficients
  vcov <- sandwich::vcovHC(model, type = "HC1")
  N <- dim(model$model)[1]

  R <- matrix(0, nrow = length(par), ncol = length(coefs))
  for (q in 1:length(par)) {R[q,unlist(par[q])] <- 1}
  if (type == "F") {
    if (is.null(N)) {stop("Please provide a value for N")}
    f <- as.numeric(t(R %*% coefs) %*%
                      solve(R %*% vcov %*% t(R)) %*%
                      (R %*% coefs)/(length(par)))
    p <- 1 - stats::pf(f, df1 = length(par),
                df2 = N - (length(coefs)))
    return(c(F = f, P = p)) }
  if (type == "W") {
    w <- as.numeric(t(R %*% coefs) %*%
                      solve(R %*% vcov %*% t(R)) %*%
                      (R %*% coefs))
    p <- 1 - stats::pchisq(w, df = length(par))
    return(c(W = w, P = p))}
}

##################################################
#' Analyzing subgroup effects with a continuous modifier
#'
#' Inputs a linear model from the `lm` function and set of indices and
#' multipliers for a set of coefficients. Return the estimate, CI and p-value.
#' NOTE: For the indices in `par`, beta0 = 1, beta1 = 2, and so on.
#'
#' @param par  A vector specifying the indices of coefficients
#' @param mults A vector of multipliers
#' @param model An lm object
#' @param alpha The alpha level for the CIs and p-values
#'
#' @return A vector containing the estimate, CI, and p-value
#' @export
#' @examples
#' # Testing beta1 + 7.5*beta3
#' # lincom(par = c(2,4),
#' #        mults = c(1,7.5),
#' #        model=my_linear_model)

lincom <- function(par, mults, model, alpha = 0.05) {

  coefs <- model$coefficients
  vcov <- sandwich::vcovHC(model, type = "HC1")
  N <- dim(model$model)[1]


  R <- matrix(0, nrow = 1, ncol = length(coefs))
  for (q in 1:length(par)) {R[1,par[q]] <- mults[q]}
  w <- sqrt(as.numeric(t(R %*% coefs) %*%
                         solve(R %*% vcov %*% t(R)) %*%
                         (R %*% coefs)))
  p <- 2*(1 - stats::pt(w, df = N - length(coefs)))
  Est <- R %*% coefs
  tol <- stats::qt(1 - alpha/2, df = N - length(coefs))
  CI.Lo <- R %*% coefs - tol*sqrt(R %*% vcov %*% t(R))
  CI.Hi <- R %*% coefs + tol*sqrt(R %*% vcov %*% t(R))
  return(c(EST = Est, CI.LO = CI.Lo, CI.HI = CI.Hi, P = p))
}
