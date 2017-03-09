# ---------------------------------------------------------------------------- #
#' resTable
#'
#' Function that reads multiple find_circ output files and returns a
#' SummarizedExperiment object
#'
#'
#' @param se a SummarizedExperiment object
#' @param ... other arguments
#'
#' @return returns a data.table
#'
#'
#' @docType methods
#' @rdname resTable-methods
#'
#' @export
setGeneric("resTable",
           function(se,
                    ...)
           standardGeneric("resTable"))

#' @aliases resTable,RangedSummarizedExperiment-method
#' @rdname resTable-methods
setMethod("resTable",
          signature("RangedSummarizedExperiment"),
          definition=function(se, ...) {
            DT <- data.table(as.data.frame(rowRanges(se)))
            setnames(DT, "seqnames", "chr")
            for (i in 1:length(colData(se)$sample)) {

              DT <- cbind(DT, sapply(names(assays(se)), function(x) assays(se)[[x]][,i]))
              setnames(DT, tail(names(DT), length(assays(se))), paste0(colData(se)$sample[i], "_", tail(names(DT), length(assays(se)))))

            }

            return(DT)
          })
