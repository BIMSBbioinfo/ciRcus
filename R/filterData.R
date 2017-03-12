# ---------------------------------------------------------------------------- #
#' filter find_circ.py-generated splice candidates
#'
#' description
#'
#' details
#'
#' @param sites splice site candidates, as read from find_circ.py output by readCircs
#' @param n_uniq_thr minimal number of unique reads required for calling a splice site
#'
qualFilter <- function(sites, n_uniq_thr = 2) {

  sites <- sites[!(grepl("circ", name) & (end - start > 100000))]
  sites <- sites[breakpoints <= 1]
  sites <- sites[anchor_overlap <= 2]
  sites <- sites[edits <= 2]
  sites <- sites[chrom != "chrM"]
  sites <- sites[n_uniq >= n_uniq_thr]

  # pipeline-specific filters
  if (ncol(sites) == 19) {
    qual.ind <- rowSums(sites[, grepl('best_qual',colnames(sites)), with = FALSE] >= 35) > 0
    sites <- sites[qual.ind]
  } else if (ncol(sites) == 21) {
    sites <- sites[!grepl("HUGE|SHORT", category)]
  }

  return(sites)
}
