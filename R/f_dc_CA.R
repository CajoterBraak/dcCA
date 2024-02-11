scores_dcCA <- function(out, display= NULL, scaling = "sites"){



  if (scaling == "sites")  myconst <- sqrt(nobs(out$RDAonEnv)*out$RDAonEnv$tot.chi) else
    if (scaling == "species") myconst <- sqrt(nobs(out$RDAonEnv))

    lc_env_scores  <- vegan::scores(out$RDAonEnv, display = c("lc"), scaling = scaling,
                                    choices = seq_len(out$RDAonEnv$CCA$rank), const = myconst)
    site_scores  <- vegan::scores(out$RDAonEnv, display = c("sites"), scaling = scaling,
                                  choices = seq_len(out$RDAonEnv$CCA$rank), const = myconst)

    out1 <- list(site_scores_unconstrained = site_scores,lc_env_scores=lc_env_scores)
    return(out1)
}

calculate_b_se_tval <- function(qr_decomp_or_X, y, w = NULL, scale2 = 0, name = "SNC", fitted_only = FALSE) {
  # specify  qr_decomp_or_X is an (yet unweigthed to-be weigthed) matrix or
  # the qr_decomp of the weigthed X matrix
  # y is the (yet unweigthed to-be weigthed) response vector or matrix

  if (is.null(w)){
    if(is.matrix(qr_decomp_or_X))w=rep(1,nrow(qr_decomp_or_X))else w = rep(1, nrow(qr_decomp_or_X$qr))
  }
  w <- w/sum(w); sqrtw <- sqrt(w)
  y_weightedw <-  as.matrix(y) * sqrtw

  TSS <- colSums(y_weightedw^2)


  # Perform QR decomposition on the weighted design matrix
  if (!is.qr(qr_decomp_or_X)) {
    X_weighted <-  X *sqrtw
    qr_decomp_or_X <- qr(X_weighted)
  }

  # Compute estimated regression coefficients
  beta_hat <- qr.coef(qr_decomp_or_X, y_weightedw)

  # Calculate fitted values
  fitted_valuesw <- qr.fitted(qr_decomp_or_X, y_weightedw)
  #species_traits_correlations<- diag(cor(fitted_valuesw, y_weightedw))# sqrt of R2 below

  # Calculate residuals
  residuals <- y_weightedw - fitted_valuesw

  fitted_values <- fitted_valuesw / sqrtw

  if (fitted_only) {out1 <- fitted_values} else {

  # Compute residual sum of squares (RSS)
  RSS <- colSums(residuals^2)

  # Estimate variance of the errors
  n <- length(w)
  p <- qr_decomp_or_X$rank
  sigma_hat_sq <- RSS / (n - p - 1 )

  # Calculate variance-covariance matrix of the estimated regression coefficients
  #R_inv <- solve(qr_decomp_or_X)

  R <- qr.R(qr_decomp_or_X)
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

  sds <- sqrt(colSums(qr.X(qr_decomp_or_X)^2))
  avg <- attr(qr_decomp_or_X$qr, which= "scaled:center")
  beta_stan <- beta_hat * sds
  colnames(beta_stan) <- paste("Regr", seq_len(ncol(beta_stan)),sep ="")

  colnames(tval) <- paste("tval", seq_len(ncol(tval)), sep="")
  b_se <- data.frame(beta= beta_hat,  se = se)
  coef_normed <- cbind(Avg = avg , SDS = sds, beta_stan, tval)
  attr(coef_normed, "meaning") <- "mean, sd, standardized regression coefficients and their t-ratio"
  out1 <- list(fitted = fitted_values, y = y, coef_normed=coef_normed, b_se = b_se,  R2 = 1-RSS/TSS)
  }
  return(out1)
}

f_trait_axes <- function(out){
# SNC lc_traits and trait regr, tval, cor
  lc_scores  <- vegan::scores(out$RDAonEnv, display = c("lc"), scaling = "species",
                              choices = seq_len(out$RDAonEnv$CCA$rank), const = sqrt(nobs(out$RDAonEnv)))

  SNC <-  (t(as.matrix(out$data$Y)) %*% lc_scores) / (out$CCAonTraits$rowsum * nobs(out$RDAonEnv))

  if (!is.null(out$CCAonTraits$pCCA)){  # orthogalize with respect to any covariate
    SNC <- SNC - calculate_b_se_tval(out$CCAonTraits$pCCA$QR, y=SNC,
                               w = out$CCAonTraits$rowsum,  scale2 = 0, name = "SNC", fitted_only = TRUE)
  }
  #print(names(out))
  res <- calculate_b_se_tval(out$CCAonTraits$CCA$QR, y=SNC,
                             w = out$CCAonTraits$rowsum,  scale2 = 1, name = "SNC")

  print(names(res))

  c_traits_normed <- res$coef_normed
  attr(c_traits_normed, which = "warning") <-"The t-values are optimistic, i.e. an underestimate of their true absolute value"


  # correlations of the dataTraits with the SNC wrt the first axis
  traits0 <-  model.matrix(~. -1, data = out$data$dataTraits)
  Cormat <- cov2cor(ade4::covwt(cbind( traits0, SNC), w= out$CCAonTraits$rowsum))
  Cor_Trait_SNC <- Cormat[seq_len(ncol(traits0)),ncol(traits0) + seq_len(out$RDAonEnv$CCA$rank) , drop = FALSE]
  colnames(Cor_Trait_SNC) <- paste("SNC-ax", seq_len(ncol(Cor_Trait_SNC)), sep= "")
  attr(Cor_Trait_SNC, which = "meaning")<- "inter set correlation, correlation between traits and SNC of axes"
  out2 <- list(species_scores = list(species_scores_unconstrained = res$y, lc_traits_scores = res$fitted), correlation = Cor_Trait_SNC, c_traits_normed= c_traits_normed, b_se= res$b_se, R2_traits = res$R2)
  return(out2)
}



