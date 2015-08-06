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
getIDs <- function(circs, organism, assembly, chrom="chrom", start="start", end="end", strand="strand") {

  con <- dbConnect(drv    = dbDriver("MySQL"),
                   host   = getOption("circbase.host"),
                   user   = getOption("circbase.user"),
                   pass   = getOption("circbase.pass"),
                   dbname = getOption("circbase.db"))

  query <- paste( "SELECT circID, chrom, pos_start, pos_end, strand FROM ",
                  organism, "_", assembly, "_circles", sep="")

  rs <- dbSendQuery(con, query)
  chunk <- data.table(fetch(rs, n=-1))
  chunk$id <- paste(chunk$chrom, ":", chunk$pos_start, "-", chunk$pos_end, sep="")

  circs$id <- paste(circs[[chrom]], ":", circs[[start]], "-", circs[[end]], sep="")
  out <- merge(circs, chunk[,.(id, circID)], by="id", all.x=T)

  dbHasCompleted(rs)
  dbClearResult(rs)
  dbListTables(con)
  dbDisconnect(con)

  return(out[, !"id", with=F])
}

# getIDs(chrom="chr3", organism="mmu", assembly="mm9")
# getIDs(chrom=circs.f$chrom, start=circs.f$start, end=circs.f$end, strand=circs.f$strand, organism="mmu", assembly="mm9")


#all_cons <- dbListConnections(MySQL())
#for(con in all_cons) dbDisconnect(con)
