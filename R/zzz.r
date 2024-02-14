#' @title The package dcCA performs double constrained correspondence analysis for trait-environment analysis in ecology
#'
#' @description
#' The \code{dcCA} package analyzes multi-trait multi-environment ecological data by
#' double constrained correspondence analysis (ter Braak et al. 2018) using \code{vegan}, \code{ade4} and native R code.
#' Throughout the two step algorithm of ter Braak et al. (2018) is used. This algorithm
#' combines and extends community- (sample-) and species-level analyses, i.e.
#' the usual community weigthed means (CWM)-based regression analysis and the
#' species-level analysis of species-niche centroids (SNC)-based regression analysis.
#' The two steps use \code{\link[vegan]{cca}} to regress the abundance data on to the traits
#' and \code{\link[vegan]{rda}} to regress the CWM of the orthonormalized traits on to the environemtal predictors.
#' The abundance data are divided by the sample total
#' (i.e. 'closed') in the vegan-based version. This
#' has the advantage that this multivariate analysis corresponds with an unweighted (multi-trait)
#' community-level analysis, instead of being weighted.
#' The current vegan-based analysis is efficient for sample-based permutation tests but slow or
#' -not yet available- for species-based permutation tests.
#' The technical reason for the closure, is that \code{vegan} \code{\link[vegan]{rda}} cannot do a weighted analysis.
#'
#'
#' @references
#'
#' ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
#' double constrained correspondence analysis.
#' Environmental and Ecological Statistics, 25(2), 171-197.
#' https://doi.org/10.1007/s10651-017-0395-x or
#' http://rdcu.be/ETPh
#'
#' @seealso \code{\link[vegan]{cca}}
#' @name dcCA
NULL