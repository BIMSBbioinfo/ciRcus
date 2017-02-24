# ---------------------------------------------------------------------------- #
#' circLinRatio
#'
#' @param bla
#'
#' @return returns a data.table
#'
#'
#'
#' @docType methods
#' @rdname summarizeCircs-methods
#'
#' @export
setGeneric("circLinRatio",
           function(se,
                    ...)
             standardGeneric("circLinRatio"))
setMethod("circLinRatio",
          signature("RangedSummarizedExperiment"),
          definition=function(se, ...) {

            assays(se)$ratio <- round(assay(se, 1) / pmax(assay(se, 2), assay(se, 3)), 2)

            return(se)
          })
