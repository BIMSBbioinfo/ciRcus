# ---------------------------------------------------------------------------- #
#' circLinRatio
#'
#' @param se a SummarizedExperiment object
#' @param ... other arguments
#'
#' @return returns a data.table
#'
#'
#'
#' @docType methods
#' @rdname circLinRatio-methods
#'
#' @export
setGeneric("circLinRatio",
           function(se,
                    ...)
             standardGeneric("circLinRatio"))

#' @aliases circLinRatio,RangedSummarizedExperiment-method
#' @rdname circLinRatio-methods
setMethod("circLinRatio",
          signature("RangedSummarizedExperiment"),
          definition = function(se, ...) {

            assays(se)$ratio <- round(assays(se)$circ /
                                        pmax(assays(se)$linear.start,
                                             assays(se)$linear.end),
                                      2)

            return(se)
          })
