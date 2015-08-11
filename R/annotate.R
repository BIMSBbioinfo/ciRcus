# library(ensembldb)
# library(EnsDb.Hsapiens.v75)
# library(data.table)

# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#'
#' details
#'
#' @param circs
#'
#' @export
annotateCircs <- function(circs) {

  require(GenomicRanges)

  genes.gr <- genes(EnsDb.Hsapiens.v75)
  genes.gr <- keepStandardChromosomes(genes.gr)
  seqlevels(genes.gr) <- paste("chr", seqlevels(genes.gr), sep="")
  seqlevels(genes.gr)[which(seqlevels(genes.gr) == "chrMT")] <- "chrM"

  # circs to GR
  circs.gr <- GRanges(seqnames=circs$chrom,
                      ranges=IRanges(start=circs$start,
                                     end=circs$end),
                      strand=circs$strand,
                      id=paste(circs$chrom, ":", circs$start, "-", circs$end, sep=""))
  circs.gr <- sort(circs.gr)


  # GR with circ starts and ends only (left and right flanks, actually)
  circ.starts.gr <- circs.gr
  end(circ.starts.gr) <- start(circ.starts.gr)
  circ.ends.gr <- circs.gr
  start(circ.ends.gr) <- end(circ.ends.gr)

  olap.start <- findOverlaps(circ.starts.gr, genes.gr, type="within")
  olap.end   <- findOverlaps(circ.ends.gr, genes.gr, type="within")

  circs$id <- paste(circs$chrom, ":", circs$start, "-", circs$end, sep="")
  circs$start.hit <- circs$id %in% circs.gr$id[queryHits(olap.start)]
  circs$end.hit   <- circs$id %in% circs.gr$id[queryHits(olap.end)]

  matches.start <- data.table(id=circs.gr$id[queryHits(olap.start)], gene=names(genes.gr)[subjectHits(olap.start)] )
  matches.end   <- data.table(id=circs.gr$id[queryHits(olap.end)],   gene=names(genes.gr)[subjectHits(olap.end)] )

  start.list <- lapply(split(matches.start, matches.start$id), function(x) x$gene)
  end.list   <- lapply(split(matches.end,   matches.end$id),   function(x) x$gene)

  circs <- merge(circs, data.table(id=names(start.list), starts=start.list), by="id", all.x=T)
  circs <- merge(circs, data.table(id=names(end.list), ends=end.list), by="id", all.x=T)

  hs <- hash(start.list)
  he <- hash(end.list)

  circs$hit.ctrl <- circs$id %in% unique(c(keys(hs), keys(he)))

  #ptm <- proc.time()
  tmphits <- integer()
  tmpgenes <- integer()
  host.candidates <- integer()
  for (circ in circs$id) {
    tmphits <- append(tmphits, sum(hs[[circ]] %in% he[[circ]]))
    tmpgenes <- append(tmpgenes, paste(hs[[circ]][hs[[circ]] %in% he[[circ]]], collapse=","))
    host.candidates <- append(host.candidates, length(unique(c(hs[[circ]], he[[circ]]))))
  }
  #proc.time() - ptm
  circs$hitcnt <- tmphits
  circs$hitgenes <- tmpgenes
  circs$host.candidates <- host.candidates

  circs$host[circs$hitcnt == 1] <- circs$hitgenes[circs$hitcnt == 1]
  circs$host[circs$hitcnt > 1]  <- "ambiguous"
  circs$host[circs$hitcnt == 0 & circs$start.hit == FALSE & circs$end.hit == FALSE] <- "intergenic" # TODO: actually, some of them may have a putative host gene within, I was only testing starts/ends
  circs$host[circs$hitcnt == 0 & circs$start.hit == TRUE  & circs$end.hit == TRUE]  <- "no_single_host"
  circs$host[circs$hitcnt == 0 & xor(circs$start.hit, circs$end.hit) & circs$host.candidates > 1] <- "ambiguous"
  circs$host[circs$hitcnt == 0 & circs$start.hit == TRUE  & circs$end.hit == FALSE & circs$host.candidates == 1] <- circs$starts[circs$hitcnt == 0 & circs$start.hit == TRUE   & circs$end.hit == FALSE & circs$host.candidates == 1]
  circs$host[circs$hitcnt == 0 & circs$start.hit == FALSE & circs$end.hit == TRUE & circs$host.candidates == 1]  <- circs$ends[circs$hitcnt == 0   & circs$start.hit == FALSE  & circs$end.hit == TRUE  & circs$host.candidates == 1]

  return(circs)
}
