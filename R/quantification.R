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

            assays(se)$ratio <- round(assay(se, 1) / pmax(assay(se, 2),
                                                          assay(se, 3)), 2)

            return(se)
          })
