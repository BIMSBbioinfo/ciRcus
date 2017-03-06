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
#' @importFrom AnnotationDbi saveDb loadDb
#' @import AnnotationHub
#' @importFrom biomaRt useMart useDataset getBM
#' @importFrom data.table data.table rbindlist dcast.data.table set setnames fread
#' @import DBI
#' @importFrom GenomicRanges makeGRangesFromDataFrame resize reduce
#' @importFrom GenomicFeatures makeTxDbFromGRanges
#' @import ggplot2
#' @importFrom hash hash keys
#' @importFrom IRanges findOverlaps
#' @import methods
#' @importFrom RColorBrewer brewer.pal
#' @import RMySQL
#' @import S4Vectors
#' @import stringr
#' @import SummarizedExperiment
NULL


.onLoad <- function(libname, pkgname) {
  op <- options()
  op.circus <- list(
    # circbase.host = "localhost",
    # circbase.user = "webuser",
    # circbase.pass = "w3b_u5er",
    circbase.host = "141.80.181.74",
    circbase.user = "circbase",
    circbase.pass = "circbase",
    circbase.db   = "circbase",
    circbase.port = 3306,

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
                           "81"      = "jul2015.archive.ensembl.org",
                           "current" = "ensembl.org"
                           ),

    ensembl.organism = list("hsa" = "hsapiens",
                            "mmu" = "mmusculus",
                            "cel" = "celegans",
                            "rno" = "rnorvegicus",
                            "dme" = "dmelanogaster"
                            ),

    assembly2annhub = list("hg19"     = "AH10684",
                           "hg38"     = "AH47963",
                           "mm10"     = "AH47973",
                           "dm6"      = "AH47953",
                           "rn5"      = "AH28841",
                           "WBcel235" = "AH47942"
                           ),

    assembly2release = list("hg19"     = "75",
                            "hg38"     = "current",
                            "mm10"     = "current",
                            "dm6"      = "current",
                            "rn5"      = "79",
                            "WBcel235" = "81"
                            ),

    assembly2organism = list("hg19"     = "hsa",
                             "hg38"     = "hsa",
                             "mm10"     = "mmu",
                             "dm6"      = "dme",
                             "WBcel235" = "cel"
                             )

  )
  toset <- !(names(op.circus) %in% names(op))
  if(any(toset)) options(op.circus[toset])

  invisible()
}
