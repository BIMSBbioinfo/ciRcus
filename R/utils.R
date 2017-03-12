#' test if input and annotation have the same coordinate indexing
#'
#' @param g1 a GenomicRanges object
#' @param g2 a GenomicRanges object
#'
#' @return a list
#' @export
#'
testCoordinateIndexing <- function(g1, g2) {
  start.hits <- c(equal    = sum( start(g1)      %in% start(g2)),
                  plusone  = sum((start(g1) + 1) %in% start(g2)),
                  minusone = sum((start(g1) - 1) %in% start(g2)))
  end.hits   <- c(equal    = sum( end(g1)        %in%   end(g2)),
                  plusone  = sum((end(g1) + 1)   %in%   end(g2)),
                  minusone = sum((end(g1) - 1)   %in%   end(g2)))
  if        (start.hits["plusone"] > sum(start.hits) * 0.9) {
    warning(round(start.hits["plusone"] / sum(start.hits) * 100, 2),  "% start coordinates are off by one when compared to annotation. Looks like input is zero- and annotation is one-based.")
  } else if (start.hits["minusone"] > sum(start.hits) * 0.9) {
    warning(round(start.hits["minusone"] / sum(start.hits) * 100, 2), "% start coordinates are off by one. Looks like input is one- and annotation is zero-based.")
  } else if (start.hits["minusone"] > sum(start.hits) * 0.9) {
    warning(round(end.hits["plusone"] / sum(end.hits) * 100, 2),  "% end coordinates are off by one when compared to annotation. Looks like input is zero- and annotation is one-based.")
  } else if (  end.hits["plusone"]  > sum(end.hits) * 0.9) {
    warning(round(end.hits["minusone"] / sum(end.hits) * 100, 2),  "% start coordinates are off by one when compared to annotation. Looks like input is one- and annotation is zero-based.")
  }

  return(list(start.hits = start.hits,
              end.hits = end.hits))
}
