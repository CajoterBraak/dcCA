
calculate_b_se_tval <- function(X_or_qr_decomp_of_X, y, w = NULL, scale2 = 0, name = "SNC", fitted_only = FALSE) {
  # specify  X_or_qr_decomp_of_X is an (yet unweigthed to-be weigthed) matrix or
  # the qr_decomp of the weigthed X matrix
  # y is the (yet unweigthed to-be weigthed) response vector or matrix

  if (is.null(w)){
    if(is.matrix(X_or_qr_decomp_of_X))w=rep(1,nrow(X_or_qr_decomp_of_X))else w = rep(1, nrow(X_or_qr_decomp_of_X$qr))
  }
  w <- w/sum(w); sqrtw <- sqrt(w)
  y_weightedw <-  as.matrix(y) * sqrtw

  TSS <- colSums(y_weightedw^2)


  # Perform QR decomposition on the weighted design matrix
  if (!is.qr(X_or_qr_decomp_of_X)) {
    X_weighted <-  X_or_qr_decomp_of_X *sqrtw
    X_or_qr_decomp_of_X <- qr(X_weighted)
  }

  # Compute estimated regression coefficients
  beta_hat <- qr.coef(X_or_qr_decomp_of_X, y_weightedw)

  # Calculate fitted values
  fitted_valuesw <- qr.fitted(X_or_qr_decomp_of_X, y_weightedw)
  #species_traits_correlations<- diag(cor(fitted_valuesw, y_weightedw))# sqrt of R2 below

  # Calculate residuals
  residuals <- y_weightedw - fitted_valuesw

  fitted_values <- fitted_valuesw / sqrtw

  if (fitted_only) {out1 <- fitted_values} else {

  # Compute residual sum of squares (RSS)
  RSS <- colSums(residuals^2)

  # Estimate variance of the errors
  n <- length(w)
  p <- X_or_qr_decomp_of_X$rank
  sigma_hat_sq <- RSS / (n - p - 1 )

  # Calculate variance-covariance matrix of the estimated regression coefficients
  #R_inv <- solve(X_or_qr_decomp_of_X)

  R <- qr.R(X_or_qr_decomp_of_X)
  R_inv <- solve(R)
  XtX_inv <- R_inv %*% t(R_inv)


  se <- matrix(NA, nrow = nrow(beta_hat), ncol = length(sigma_hat_sq))
  for (i in seq_along(sigma_hat_sq)){
    var_covar_matrix <- XtX_inv * sigma_hat_sq[i]
    # Calculate standard errors
    se[,i] <- sqrt(diag(var_covar_matrix))
  }


  TSSfit <- colSums(fitted_valuesw^2)


  if (scale2){
    sqrtTSS <- sqrt(TSSfit/scale2)
    if (length(TSS)>1) sTSS <- diag(1/sqrtTSS) else sTSS =  diag(1) * (1/sqrtTSS)
    fitted_values <- fitted_values %*% sTSS
    y <- y %*% sTSS
    beta_hat <- beta_hat %*% sTSS
    se <- se %*% sTSS
  }


  if (length(TSS)>1) colnames(y)<- paste(name, seq_len(ncol(y)), sep ="")
  if (length(TSS)>1) colnames(fitted_values)<- paste(name,"_lc", seq_len(ncol(y)), sep ="")

  #SNC <- cbind(y, fitted_values)
  tval <- beta_hat/se

  sds <- sqrt(colSums(qr.X(X_or_qr_decomp_of_X)^2))
  avg <- attr(X_or_qr_decomp_of_X$qr, which= "scaled:center")
  VIF <- diag(XtX_inv)*sds^2
  beta_stan <- beta_hat * sds
  colnames(beta_stan) <- paste("Regr", seq_len(ncol(beta_stan)),sep ="")

  colnames(tval) <- paste("tval", seq_len(ncol(tval)), sep="")
  b_se <- data.frame(beta= beta_hat,  se = se)
  coef_normed <- cbind(Avg = avg , SDS = sds, VIF = VIF, beta_stan, tval)
  attr(coef_normed, "meaning") <- "mean, sd, VIF, standardized regression coefficients and their t-ratio"
  out1 <- list(fitted = fitted_values, y = y, coef_normed=coef_normed, b_se = b_se,  R2 = 1-RSS/TSS)
  }
  return(out1)
}




