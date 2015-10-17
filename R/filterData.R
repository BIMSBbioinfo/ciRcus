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
qualFilter <- function(sites, n_uniq_thr=2) {

  sites <- sites[!(grepl("circ", sites$name) & (sites$end - sites$start > 100000)),]
  sites <- sites[sites$breakpoints <= 1,]
  sites <- sites[sites$anchor_overlap <= 2,]
  sites <- sites[sites$edits <= 2,]
  sites <- sites[sites$chrom != "chrM",]
  sites <- sites[sites$n_uniq >= n_uniq_thr,]
  sites <- sites[sites$best_qual_A >= 35 | sites$best_qual_B >= 35,]

  return(sites)
}
