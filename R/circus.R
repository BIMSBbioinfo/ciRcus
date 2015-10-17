#' ciRcus: An R package for circRNAs munging
#'
#' The foo package provides three categories of important functions:
#' foo, bar and baz.
#'
#' @section Foo functions:
#' The foo functions ...
#'
#' @docType package
#' @name ciRcus
#' @import AnnotationHub
#' @import biomaRt
#' @import data.table
#' @import DBI
#' @import hash
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import ggplot2
#' @import IRanges
#' @import RMySQL
NULL


.onLoad <- function(libname, pkgname) {
  op <- options()
  op.circus <- list(
    circbase.host = "localhost",
    circbase.user = "webuser",
    circbase.pass = "w3b_u5er",
    # circbase.host = "141.80.186.72",
    # circbase.user = "circbaseuser",
    # circbase.pass = "01circbaseuser",
    circbase.db   = "circbase",

    ensembl.release = list("54"      = "may2009.archive.ensembl.org",
                           "67"      = "may2012.archive.ensembl.org",
                           "68"      = "jul2012.archive.ensembl.org",
                           "69"      = "oct2012.archive.ensembl.org",
                           "70"      = "jan2013.archive.ensembl.org",
                           "71"      = "apr2013.archive.ensembl.org",
                           "72"      = "jun2013.archive.ensembl.org",
                           "73"      = "sep2013.archive.ensembl.org",
                           "74"      = "dec2013.archive.ensembl.org",
                           "75"      = "feb2014.archive.ensembl.org",
                           "76"      = "aug2014.archive.ensembl.org",
                           "77"      = "oct2014.archive.ensembl.org",
                           "78"      = "dec2014.archive.ensembl.org",
                           "79"      = "mar2015.archive.ensembl.org",
                           "80"      = "may2015.archive.ensembl.org",
                           "current" = "ensembl.org"
                           ),

    ensembl.organism = list("hsa" = "hsapiens",
                            "mmu" = "mmusculus",
                            "cel" = "celegans",
                            "rno" = "rnorvegicus",
                            "dme" = "dmelanogaster"
                            ),

    assembly2annhub = list("hg19" = "AH10684",
                           "hg38" = "AH47963",
                           "mm10" = "AH47973",
                           "dm6"  = "AH47953",
                           "rn5"  = "AH28841"
                           ),
    assembly2release = list("hg19" = "75",
                            "hg38" = "current",
                            "mm10" = "current",
                            "dm6"  = "current",
                            "rn5"  = "79"
                            ),
    assembly2organism = list("hg19" = "hsa",
                             "hg38" = "hsa",
                             "mm10" = "mmu",
                             "dm6"  = "dme"
                             )

  )
  toset <- !(names(op.circus) %in% names(op))
  if(any(toset)) options(op.circus[toset])

  invisible()
}
