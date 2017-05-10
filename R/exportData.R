# ---------------------------------------------------------------------------- #
#' bedTracks
#'
#' Export circRNA information from given \code{SummarizedExperiment} object to
#' \code{BED6}-style \code{GRangesList}.
#' This is useful to generate genome browser tracks showing circRNA candidates
#' per sample.
#'
#'
#' @param se \code{SummarizedExperiment} object with circRNA information
#' @param score \code{character} vector (only first element will be used) naming
#'              the assay to use as \code{BED} score (\code{NULL} to omit BED
#'              score and use `.` as a placeholder instead)
#' @param min.score \code{numeric} vector (only first element will be used) with
#'                  minimal score a circRNA must have to be included in the
#'                  \code{BED} output for a given sample (\code{NULL} to include
#'                  all circRNAs for all samples, even if they were not detected
#'                  in that sample)
#' @param max.score \code{numeric} vector (only first element will be used) with
#'                  maximal score to be used in \code{BED} output (higher score
#'                  will be truncated) (\code{NULL} to keep scores unlimited,
#"                  ignoring the \code{BED} format definition)
#'
#' @return returns a GRangesList object with BED6 circRNA data per sample
#'
#'
#' @docType methods
#' @rdname bedTracks-methods
#'
#' @export
setGeneric("bedTracks",
           function(se,
                    score = "circ.uniq",
                    min.score = unlist(ifelse(is.null(score), list(NULL),
                                              list(1))),
                    max.score = unlist(ifelse(is.null(score), list(NULL),
                                              list(1000))))
           standardGeneric("bedTracks"))

#' @aliases bedTracks,RangedSummarizedExperiment-method
#' @rdname bedTracks-methods
setMethod("bedTracks",
          signature("RangedSummarizedExperiment"),
          definition = function(se, score, min.score, max.score) {

            # check input
            if (length(score) > 1) {
              warning("length(score) > 1; only first entry will be used")
              score <- score[1]
            }
            if (length(min.score) > 1) {
              warning("length(min.score) > 1; only first entry will be used")
              min.score <- min.score[1]
            }
            if (length(max.score) > 1) {
              warning("length(max.score) > 1; only first entry will be used")
              max.score <- max.score[1]
            }
            if (!is.null(score)) {
              if (!(score %in% names(assays(se)))) {
                warning(paste0("no assay named '", score, "'",
                               "; BED output will be generated without scores"))
                score <- NULL
              }
            }
            if (is.null(score) & !is.null(min.score)) {
              warning("no BED score defined; circRNAs will not be filtered")
              min.score <- NULL
            }
            if (is.null(score) & !is.null(max.score)) {
              warning("no BED score defined; BED scores will not be truncated")
              min.score <- NULL
            }

            # get circRNA coordinates
            ranges <- rowRanges(se)

            # drop metadata columns
            mcols(ranges) <- NULL

            # get score matrix (circRNA x sample)
            if (!is.null(score)) {
              scores <- assays(se)[[score]]

            # use `.` as score of no score assay was specified (1 x sample)
            } else {
              message("no BED score defined; using `.` as placeholder")
              scores <- matrix(rep(".", ncol(se)), nrow = 1)
              colnames(scores) <- colnames(se)
            }

            # setup GRangesList with circRNA coordinates scored per sample
            sample.names <- colnames(scores)
            ranges <- lapply(sample.names,
                             function(sample.name) {
                               mcols(ranges)$score <- scores[, sample.name]
                               ranges
                             })
            names(ranges) <- sample.names
            ranges <- GRangesList(ranges)

            # filter out circRNA not passing the minimal score (per sample)
            if (!is.null(min.score))
              ranges <- endoapply(ranges,
                                  function(gr) {
                                    subset(gr,
                                           rtracklayer::score(gr) >= min.score)
                                  })

            # truncate score at given maximum
            if (!is.null(max.score))
              ranges <- endoapply(ranges,
                                  function(gr) {
                                    mcols(gr)$score[rtracklayer::score(gr) >
                                                    max.score] <- max.score
                                    gr
                                  })

            # return GRangesList
            return(ranges)
          })


# ---------------------------------------------------------------------------- #
#' writeBedTracks
#'
#' Write circRNA information from given \code{BED6}-style \code{GRangesList} (or
#' \code{SummarizedExperiment} object) to \code{BED} files (one per sample).
#' These can be loaded as genome browser tracks showing circRNA candidates.
#'
#'
#' @param circs \code{GRangesList} or \code{SummarizedExperiment} object with
#'              circRNA information
#' @param out.prefix \code{character} vector (only first element will be used)
#'                   specifying the prefix to use for output files (before the
#'                   sample name)
#' @param out.suffix \code{character} vector (only first element will be used)
#'                   specifying the suffix to use for output files (after the
#'                   sample name)
#' @param seqlevels.style \code{character} vector (only first element will be
#'                        used) specifying the seqlevels style to use for the
#'                        output BED files
#' @param ... named arguments defined above to be passed on from
#'            \{code{SummarizedExperiment}-method to \code{GRangesList}-method
#'
#' @return None
#'
#'
#' @docType methods
#' @rdname writeBedTracks-methods
#'
#' @export
setGeneric("writeBedTracks",
           function(circs, ...)
           standardGeneric("writeBedTracks"))

#' @aliases writeBedTracks,GRangesList-method
#' @rdname writeBedTracks-methods
setMethod("writeBedTracks",
          signature("GRangesList"),
          definition = function(circs,
                                out.prefix = "ciRcus_",
                                out.suffix = ".bed",
                                seqlevels.style = "UCSC") {

            # check input
            if (length(seqlevels.style) > 1){
              warning(paste("length(seqlevels.style) > 1;",
                            "only first entry will be used"))
              seqlevels.style <- seqlevels.style[1]
            }

            # adjust seqlevels style
            if (!is.null(seqlevels.style)) {
              if (!(seqlevels.style) %in% seqlevelsStyle(circs)){
                message(paste0("Changing seqlevels style to '", seqlevels.style,
                               "' for BED track output file",
                               ifelse(length(circs) > 1, "s", ""), "."))
                seqlevelsStyle(circs) <- seqlevels.style
              }
            }

            # export GRanges to BED files one sample at a time
            for (sample.name in names(circs)) {
              export(circs[[sample.name]],
                     con = paste0(out.prefix, sample.name, out.suffix),
                     format = "BED"
                     )
            }
          })

#' @aliases writeBedTracks,RangedSummarizedExperiment-method
#' @rdname writeBedTracks-methods
setMethod("writeBedTracks",
          signature("RangedSummarizedExperiment"),
          definition = function(circs, ...) {

            # extract BED6-style GRangesList from SummarizedExperiment and
            # write BED tracks based on that GRangesList
            writeBedTracks(bedTracks(circs), ...)
          })
