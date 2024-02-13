#' @title Performs a vegan-based double constrained correspondence Analysis (dc-CA)
#'
#' @description
#' \code{dcCA-vegan} performs double constained correspondence analysis (dc-CA) by
#' the two-step algorithm of ter Braak et al. (2018). This algorithm
#' combines and extends community- (sample-) and species-level analyses.
#' The first step uses \code{\link[vegan]{cca}} to regress the transposed abundance data on to the traits
#' and \code{\link[vegan]{rda}} to regress the community-weighted means (CWMs)
#' of the ortho-normalized traits on to the environmental predictors.
#' In this \code{vegan}-based version, the abundance data are divided by the sample total,
#' i.e. 'closed' in the terminology of compositional data analysis (CoDa). This
#' has the advantage that the results of a dc-CA correspond with unweighted (multi-trait)
#' community-level analyses, instead of corresponding to a weighted analysis (dc-CA weighs the
#' rows by the sample total, which is 1 after closure).
#' The current vegan-based analysis is efficient for sample-based permutation tests but slow or
#' -not yet available- for species-based permutation tests.
#' The technical reason for the closure is that \code{vegan} \code{\link[vegan]{rda}} cannot do a weighted analysis.
#'
#' @param formulaEnv formula or one-sided formula for the rows (samples) with row predictors in \code{dataEnv}.
#' When two-sided, the left hand side of the formula is not used. Specify row covariates (if any ) by adding \code{+ Condition(covariate-formula)}
#' to \code{formulaEnv} as in \code{\link[vegan]{rda}}. Default: \code{~.}, i.e. all variables in \code{dataEnv} are predictor variables.
#' @param formulaTraits  formula or one-sided formula for the columns (species) with colum predictors in \code{dataTraits}.
#' When two-sided, the left hand side of the formula is not used. Specify column covariates (if any ) by  adding \code{+ Condition(covariate-formula)}
#' to \code{formulaTraits} as in \code{\link[vegan]{cca}}.Default: \code{~.}, i.e. all variables in \code{dataTraits} are predictor traits.
#' @param response matrix or data frame of the abundance data (dimension \emph{n} x \emph{m}).
#' Rownames of \code{response}, if any, are carried through.
#' @param dataEnv matrix or data frame of the row predictors, with rows corresponding to those in \code{response}.
#' ((dimension \emph{n} x \emph{p}).
#' @param dataTraits matrix or data frame of the column predictors,
#'  with rows corresponding to the columns in \code{response}.((dimension \emph{m} x \emph{q}).
#' @param dc_CA_vegan_object  optional object from an earlier run of this function. Useful if the
#' same formula for the columns (\code{formulaTraits}), \code{dataTraits} and \code{response} are used
#' with a new formula for the rows. If set, the data of the previous run is used and the result of its first step
#' is taken for the new analysis.
#' @param verbose logical for printing a simple summary (default: FALSE)
#
#' @details
#' Empty (all zero) rows and columns in \code{response} are removed from the \code{response} and the corresponding
#' rows from \code{dataEnv} and \code{dataTraits}.
#' After 'closure' (division of the values in \code{response}
#' by their row sums), the subsequent algorithm follows the two-step algorithm of ter Braak et al. (2018).
#' and consits of two steps. First, the transpose of the \code{response} is regressed on to the traits
#' (the column predictors) using \code{\link[vegan]{cca}}
#' with \code{formulaTraits}.
#' The column scores of this analysis (in scaling 2) are community weigthed means (CWM) of the
#' orthonormalized traits.
#' These are then regressed on the environmental (row) predictors using \code{\link[vegan]{rda}} with
#' with \code{formulaEnv}.
#'
#'  All row based (sample-based) subsequent analyses, sample scores and permutation tests
#' can be obtained by the appropriate functions with argument \code{object$RDAonEnv}.
#'
#' The reason for the closure are
#' two-fold: dc-CA is the efficient summary of CWM-based regression analyses, which are rarely weighted, and
#' \code{vegan} \code{\link[vegan]{rda}} cannot do a weighted analysis, whereas \code{\link[vegan]{cca}} uses
#' the weights implied by the \code{response} after closure.
#'
#' The scale-free statistics in the example \code{dune_dcCA.r},
#'  have been checked against the results in Canoco 5.15 (ter Braak & Smilauer, 1918).
#'  But, the sites and species scores and, hence, the regression coefficients, differ
#'  due to a difference in scaling of scores. But, note that the scaling is
#'   immaterial for the interpretation of the results.
#'
#' @returns
#' A list of \code{class} \code{dccaV}; that is a list with elements
#' \describe{
#' \item{CCAonTraits}{a \code{\link[vegan]{cca.object}} from the \code{\link[vegan]{cca}} analysis
#' of the transpose of the closed \code{response} using formula \code{formulaTraits}. }
#' \item{formalaTraits}{the argument \code{formulaTraits}. }
#' \item{data}{a list of \code{Y} (response data after removing empty rows and columns and after closure)
#' and \code{dataEnv} and \code{dataTraits}.}
#' \item{RDAonEnv}{a \code{\link[vegan]{cca.object}} from the \code{\link[vegan]{rda}} analysis
#' of the column scores of the \code{cca}, which are the CWMs of orthonormalized traits,
#' using formula \code{formulaEnv}. }
#' \item{formalaEnv}{the argument \code{formulaEnv}. }
#' \item{eigenvalues}{the dc-CA eigenvalues (same as those of the \code{\link[vegan]{rda}} analysis)}
#' \item{inertia}{a vector with four inertias (weighted variances):
#' \itemize{
#' \item total: the total inertia.
#' \item conditionT: the inertia explained by the condition in \code{formulaTraits} (neglecting row constraints).
#' \item traits_explain: the inertia explained by the traits (neglecting the row predictors and any
#' condition in \code{formulaTraits}).
#' This is the maximum that the row predictors could explain in dc-CA
#' (the sum of the following two items is thus less than this value).
#' \item conditionE: the trait-constrained inertia explained by the condition in \code{formulaEnv}.
#' \item constraintsE: the trait-constrained inertia explained by the predictors (without the row covariates).
#' }
#'  }
#' }
#' if \code{verbose} is \code{TRUE} there are three more items (in this version).
#' \itemize{
#' \item \code{c_traits_normed}: mean, sd, VIF and (regression) coefficients of
#'  the traits that define the dc-CA trait axes (composite traits), and their optimistic t-ratio.
#' \item \code{c_env_normed}:  mean, sd, VIF and (regression) coefficients of the environmental variables that define the dc-CA axes
#'  in terms of the environmental variables (composite gradients), and their optimistic t-ratio.
#' \item \code{species_axes}: a list with three items
#'  \itemize{
#'  \item \code{species_scores}: a list with names \code{c("species_scores_unconstrained", "lc_traits_scores")} with the
#'  matrix with species niche centroids along the dc-CA axes (composite gradients) and
#'  the matrix with linear combinations of traits.
#'  \item \code{correlation}: a matrix with inter-set correlations of the traits with their SNCs.
#'  \item \code{b_se}: a matrix with (unstandardized) regression coefficients for traits and their optimistic standard errors.
#'  \item \code{R2_traits}: a vector with coefficient of determination (R2) of the SNCs on to the traits.
#'  The square-root thereof could be called the species-trait correlation in analogy with
#'  the species-environment correlation in CCA.
#'  }
#'
#' }
#'
#' @references
#' ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
#' double constrained correspondence analysis.
#' Environmental and Ecological Statistics, 25(2), 171-197.
#' https://doi.org/10.1007/s10651-017-0395-x or
#' http://rdcu.be/ETPh
#'
#' ter Braak C.J.F. and  P. Šmilauer  (2018). Canoco reference manual
#' and user's guide: software for ordination (version 5.1x).
#' Microcomputer Power, Ithaca, USA, 536 pp.
#' @seealso \code{\link{scores.dccav}} and \code{\link{print.dccav}}
#' @example demo/dune_dcCA.R
#' @export

dc_CA_vegan <- function(formulaEnv = ~., formulaTraits = ~., response =NULL, dataEnv, dataTraits= NULL, dc_CA_vegan_object  = NULL, verbose = FALSE) {
  # response matrix or data frame, dataEnv and dataTraits data frames in which formualaE and formulaT are evaluated
  #dc_CA_vegan_object = result (value) of a previous run, can be used to save computing time for
  # runs that modify the formula for samples (step2: RDAonEnv) only
  # The step1 (CCAonTraits and the data and formulaTraits) are taken from dc_CA_vegan_object into the new result.
  # If set, formulaTraits, response, dataEnv, dataTraits are not used at all and have no efffect on the result
  change_reponse  <- function(f, response){ # used in dc_CA_vegan
    # response : character
    ft <- as.character(f)
    if (!is.character(response)) stop("response must be character") else if (length(response)>1) stop("response must be of length one")
    # ft2 <- paste(response, ft[1], ft[-c(1,2)], collapse = " ")
    if (ft[1] == "~")  ft2 <- paste0(response, ft[1], ft[-1], collapse = " ") else
      ft2 <- paste0(response, ft[1], ft[-c(1,2)], collapse = " ")

    f2 <- as.formula(ft2)
    return(f2)
  }

  if (is.null(dc_CA_vegan_object)){
    #  check and amend: make sure there are no empty rows or columns -----------------------------------------------------------------------
    if (is.null(dataTraits)) stop("dataTraits must be specified in dc_CA_vegan")
    if (!is.matrix(response)) response <- as.matrix(response) else stop("response (matrix or df) must specified")
    id0 <-1
    while(length(id0)){
      TotR <- rowSums(response)
      id0 <- which(TotR == 0)
      if (length(id0)){
        response <- response[-id0,]
        dataEnv  <- dataEnv[-id0,]
      }
      TotC <- colSums(response)
      id0 <- which(TotC == 0)
      if (length(id0)){
        response <- response[,-id0]
        dataTraits <- dataTraits[-id0, ]
      }
    }
    # end of check -----------------------------------------------------------------------
    # close the data (divide by the row total, to get strictly compositional data) -----------------------------------------------------------------------


    Abun_frac <- sweep(response, 1, STATS = TotR, FUN = '/')
    TotRfrac <- rowSums(Abun_frac)
    tY <- t(Abun_frac)
    formulaTraits <- change_reponse(formulaTraits, "tY")
    environment(formulaTraits)<- environment()
    step1 <-vegan::cca(formulaTraits, data = dataTraits)
    data= list(Y = Abun_frac, dataEnv = dataEnv, dataTraits = dataTraits)
    out1 <- list(CCAonTraits = step1,
                 formulaTraits= formulaTraits,
                 data = list(Y = Abun_frac, dataEnv = dataEnv, dataTraits = dataTraits)
    )
  } else {
    step1 <- dc_CA_vegan_object$CCAonTraits
    out1 <- dc_CA_vegan_object[c("CCAonTraits", "formulaTraits","data")]
  }
  n <- nrow(out1$data$Y)
  CWMs_orthonormal_traits <- vegan::scores(step1, display= "species",
                scaling = "species", choices = 1:step1$CCA$rank) * sqrt((n-1)/n)
  if (rownames(CWMs_orthonormal_traits)[1]=="col1") rownames(CWMs_orthonormal_traits) <- paste("Sam", seq_len((nrow(out1$data$dataEnv))),sep="")
  formulaEnv <- change_reponse(formulaEnv, "CWMs_orthonormal_traits")
  environment(formulaEnv)<- environment()
  step2 <- vegan::rda(formulaEnv, data = out1$data$dataEnv)
  inertia <- cbind(c(total= step1$tot.chi,conditionT = step1$pCCA$tot.chi, traits_explain = step1$CCA$tot.chi, conditionE = step2$pCCA$tot.chi, constraintsE = step2$CCA$tot.chi ))
  colnames(inertia)<- "weigthed variance"
  expla <- c("total inertia","inertia of the trait condition", "trait-constrained inertia",
       "trait-constrained inertia explained by the condition in formulaEnv",
       "trait-constrained inertia explained by the predictors in formulaEnv")
  names(expla) <- c("total","conditionT","traits_explain","conditionE","constraintsE")
  attr(inertia, which = "meaning") <-  matrix( expla[rownames(inertia)], ncol=1,
                                                dimnames = list(rownames(inertia),"meaning"))


  out <- c(out1, list(RDAonEnv = step2,
                      formulaEnv = formulaEnv,
                      eigenvalues =  step2$CCA$eig,
                      inertia = inertia))
  class(out) <- c("dccav", "list")
  if (verbose) out<-print.dccav(out)
  return(out)
}

