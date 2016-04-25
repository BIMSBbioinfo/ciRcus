# ---------------------------------------------------------------------------- #
#' resTable
#'
#' Function that reads multiple find_circ output files and returns a
#' SummarizedExperiment object
#'
#'
#'
#' @return returns a data.table
#'
#'
#' @docType methods
#' @rdname summarizeCircs-methods
#'
#' @export
setGeneric("resTable",
           function(se,
                    ...)
           standardGeneric("resTable"))
setMethod("resTable",
          signature("RangedSummarizedExperiment"),
          definition=function(se, ...) {
            DT <- data.table(as.data.frame(rowRanges(se)))
            setnames(DT, "seqnames", "chr")
            for (i in 1:length(colData(se)$sample)) {
              DT <- cbind(DT, assays(se)$circ[,i], assays(se)$linear.start[,i], assays(se)$linear.end[,i])
              setnames(DT, c("V2", "V3", "V4"), paste0(colData(se)$sample[i], c("_circ", "_lin.start", "_lin.end")))
            }

            return(DT)
          })
