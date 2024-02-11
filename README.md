# dcCA

<!-- badges: start -->
<!-- badges: end -->

 dcCA analyzes multi-trait multi-environment ecological data by
double constrained correspondence analysis (ter Braak et al. 2018) 
using \code{vegan}, \code{ade4} and native R code.
Throughout the two step algorithm of ter Braak et al. (2018) is used. This algorithm
combines and extends community- (sample-) and species-level analyses, i.e.
the usual community weigthed means (CWM)-based regression analysis and the
species-level analysis of species-niche centroids (SNC)-based regression analysis.
The two steps use \code{\link[vegan]{cca}} to regress the abundance data on to the traits
and \code{\link[vegan]{rda}} to regress the CWM of the orthonormalized traits on to the environemtal predictors.
The abundance data are divided by the sample total
(i.e. 'closed') in the vegan-based version. This
has the advantage that this multivariate analysis corresponds with an unweighted (multi-trait)
community-level analysis, instead of being weighted, which may give a puzzling differences
between common univarite and this multivariate analysis.
Reference: ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots fordouble constrained correspondence analysis. Environmental and Ecological Statistics, 25(2), 171-197. https://doi.org/10.1007/s10651-017-0395-x

## Installation

You can install the released version of dcCA from github by
invoking the R-code within an R-console:

``` r
install.packages("remotes")
remotes::install_github("CajoterBraak/dcCA")
```

