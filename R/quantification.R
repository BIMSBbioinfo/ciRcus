# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#'
#' details
#'
#' @param sites splice site candidates, as read from find_circ.py output by readCircs
#' @param label.circ a string, label used in the name column to distinguish circular from linear splice sites
#' @param label.norm a string, label used in the name column to distinguish linear from circular splice sites
#'
#' @export
circLinRatio <- function(sites, label.circ="circ", label.norm="norm") {

  # split sites.bed into circular and linear subset
  sites.circ <- sites[grep(label.circ, sites$name)]
  setkeyv(sites.circ, cols=c("chrom", "start", "end"))
  sites.norm <- sites[grep(label.norm, sites$name)]
  setkeyv(sites.norm, cols=c("chrom", "start", "end"))

  # intron.start is the start coordinate of norms,
  # and the end coordinate of circs
  sites.circ$intron.start <- paste(sites.circ[,chrom], sites.circ[,end], sep=":")
  sites.circ$intron.end   <- paste(sites.circ[,chrom], sites.circ[,start], sep=":")

  sites.norm$intron.start <- paste(sites.norm[,chrom], sites.norm[,start], sep=":")
  sites.norm$intron.end   <- paste(sites.norm[,chrom], sites.norm[,end], sep=":")

  sites.norm <- merge(sites.norm, sites.norm[, sum(n_reads), by=intron.start], by="intron.start")
  setnames(sites.norm, "V1", "n_right")
  sites.norm <- merge(sites.norm, sites.norm[, sum(n_reads), by=intron.end], by="intron.end")
  setnames(sites.norm, "V1", "n_left")
  sites.starts  <- sites.norm[,.(intron.start, n_right)]
  sites.ends <- sites.norm[,.(intron.end, n_left)]

  sites.circ <- merge(sites.circ, sites.starts, by="intron.start", all.x=T)
  setkeyv(sites.circ, c("chrom", "start", "end"))
  sites.circ <- merge(sites.circ, sites.ends, by="intron.end", all.x=T)
  setkeyv(sites.circ, c("chrom", "start", "end"))
  sites.circ <- sites.circ[!duplicated(sites.circ)]

  sites.circ[["n_left"]][is.na(sites.circ[["n_left"]])] <- 0
  sites.circ[["n_right"]][is.na(sites.circ[["n_right"]])] <- 0
  sites.circ$max <- apply(sites.circ[, c("n_right", "n_left"), with=F], 1, max)

  sites.circ$ratio <- round(sites.circ$n_reads / sites.circ$max, digits = 2)

  return(sites.circ[, !(names(sites.circ) %in% c("intron.start", "intron.end")), with=F])
}
