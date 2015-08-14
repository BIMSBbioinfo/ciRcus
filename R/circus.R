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
#' @import hash
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import IRanges
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
