#' Obtain Column (Sites) Related Scores, notably regression coefs and their optimistic t-ratios
#' @param out object from \code{\link{dc_CA_vegan}}
#' @param which_cor character or names of environmental variables
#' for which inter-set correlations must calculated.
#' Default: "in_model" for all environmental variables in the model specified by \code{formulaEnv}.
#' @noRd
#' @export
f_env_axes <- function(out, which_cor = "in model"){
  # which_cor character or names of environmental variables
  # for which inter-set correlations must calculated.
  #  Default "in_model" for all environmental variables in the model

  #print(names(out))
  myconst <- sqrt(stats::nobs(out$RDAonEnv)*out$RDAonEnv$tot.chi)
  CWM <- as.matrix(vegan::scores(out$RDAonEnv, display = c("sites"), scaling =  "sites",
                choices = seq_len(Rank_mod(out$RDAonEnv)), const = myconst))
  res <- calculate_b_se_tval(get_QR(out$RDAonEnv), y=CWM,
                             w = stats::weights(out$RDAonEnv, "sites"),  scale2 = 0, name = "CWM")
  c_env_normed <- res$coef_normed
  attr(c_env_normed, which = "warning") <-"The t-values are optimistic, i.e. an underestimate of their true absolute value"


  # correlations of the dataEnv with the CWMs wrt the  axes
  if (which_cor == "in model") in_model <- colnames(out$data$dataEnv)%in% colnames(attr(stats::terms(out$RDAonEnv), which = "factors")) else
    in_model = which_cor
  env0 <-  stats::model.matrix(~.-1, constrasts = FALSE, data = out$data$dataEnv[, in_model, drop= FALSE])
  #env0 <-  model.matrix(~. -1, constrasts = FALSE, data = out$data$dataEnv)
  Cormat <- stats::cov2cor(ade4::covwt(cbind( env0, CWM), w= stats::weights(out$RDAonEnv, "sites")))
  Cor_Env_CWM <- Cormat[seq_len(ncol(env0)),ncol(env0) + seq_len(Rank_mod(out$RDAonEnv)) , drop = FALSE]
  colnames(Cor_Env_CWM) <- paste("CWM-ax", seq_len(ncol(Cor_Env_CWM)), sep= "")
  attr(Cor_Env_CWM, which = "meaning")<- "inter set correlation, correlation between environmental variables and CWM of axes"

  out2 <- list( c_env_normed= c_env_normed, b_se= res$b_se, R2_env = res$R2, correlation = Cor_Env_CWM)
  return(out2)
}
