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
#' @import data.table
#' @import DBI
#' @import RMySQL
#' @import ensembldb
#' @import EnsDb.Hsapiens.v75
#' @import hash
#' @import GenomicRanges
#' @import IRanges
# @import ensembldb
# @import EnsDb.Hsapiens.v75
NULL


.onLoad <- function(libname, pkgname) {
  op <- options()
  op.circus <- list(
    circbase.host = "localhost",
    circbase.user = "webuser",
    circbase.pass = "w3b_u5er",
    circbase.db   = "circbase"
  )
  toset <- !(names(op.circus) %in% names(op))
  if(any(toset)) options(op.circus[toset])

  invisible()
}
