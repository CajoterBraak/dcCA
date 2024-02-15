#' Takes the scores function to the dcCA library
#'
#' @param x in package \code{dcCA}, an object of \code{\link{dc_CA_vegan}}, in package \code{vegan}
#' an  ordination object, see also \code{\link[vegan]{cca.object}}.
#' @seealso \code{\link{scores.dccav}},\code{\link[vegan]{scores}}
#' @export
"scores" <-
function(x, ...) UseMethod("scores")

