# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#'
#' details
#'
#' @param chrom
#' @param start
#' @param end
#'
#' @export
#'
getIDs <- function(chrom, start, end, strand, organism, assembly) {

  options("ciRcus-circBase-host" = "localhost")
  options("ciRcus-circBase-user" = "webuser")
  options("ciRcus-circBase-pass" = "w3b_u5er")
  options("ciRcus-circBase-db"   = "circbase")

  db <- dbConnect(drv    = dbDriver("MySQL"),
                  host   = getOption("ciRcus-circBase-host"),
                  user   = getOption("ciRcus-circBase-user"),
                  pass   = getOption("ciRcus-circBase-pass"),
                  dbname = getOption("ciRcus-circBase-db"))

  query <- paste( "SELECT circID, chrom, pos_start, pos_end, strand FROM ",
                  organism, "_", assembly, "_circles ",
                  "WHERE chrom=\"", chrom, "\"",
                  " limit 10", sep="")
  print(query)
  rs <- dbSendQuery(db, query)
  chunk <- data.table(fetch(rs, n=-1))
  chunk$id <- paste(chunk$chrom, ":", chunk$pos_start, "-", chunk$pos_end, sep="")

  return(chunk)
}

#getIDs(chrom="chr3", organism="mmu", assembly="mm9")


#all_cons <- dbListConnections(MySQL())
#for(con in all_cons) dbDisconnect(con)
