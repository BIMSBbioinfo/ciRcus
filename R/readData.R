# ---------------------------------------------------------------------------- #
#' read a tabular circRNA candidate list
#'
#' description
#'
#' details
#' @param file location of the input file, a character string such as:
#'             "/home/user/find_circ_sites.bed"
#' @param subs a character string, keep only lines containing it in the name column
#' @param qualFilter a boolean, tells whether the quality filtering should be performed
#' @param keepCols a vector of column numbers return
#'
readCircs <- function(file, subs="all", qualfilter=TRUE, keepCols=1:6, ...) {

  suppressWarnings(
    DT <- fread(file, sep="\t", header = T) # maybe add colClasses later
  )
  setnames(DT, "# chrom", "chrom")
  DT <- DT[!grepl("#", DT$chrom)]

  # change column classes where needed
  # ***due to find_circ.py logic of putting a header line
  #    in the middle of the output file, all columns are
  #    character after fread()
  char.class = c('chrom','name','strand','tissues','signal','strandmatch','category')
  for (col in setdiff(colnames(DT),char.class)){
    set(DT, j=col, value=as.integer(DT[[col]]))
  }

  if (subs != "all") {
    DT <- DT[grep(subs, DT$name)]
  }

  if (qualfilter == TRUE) {
    DT <- qualFilter(DT)
  }

  DT <- DT[, keepCols, with=F]
  DT$id <- paste(DT$chrom, ":", DT$start, "-", DT$end, sep="")

  return(DT[, !"id", with=F])
}
