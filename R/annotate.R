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
loadAnnotation <- function(txdb.file) {

  require(GenomicFeatures)

  txdb <- loadDb(txdb.file)

  return(txdb)
}

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
annotateHostGenes <- function(circs, txdb) {

  require(GenomicFeatures)

  genes.gr <- genes(txdb)

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

  circs <- merge(circs, data.table(id=names(start.list), starts=sapply(start.list, function(x) paste(x, collapse=","))), by="id", all.x=T)
  circs <- merge(circs, data.table(id=names(end.list), ends=sapply(end.list, function(x) paste(x, collapse=","))), by="id", all.x=T)

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

  #return(circs)
  return(circs[, !c("id", "start.hit", "end.hit", "starts", "ends", "hit.ctrl", "hitcnt", "hitgenes", "host.candidates"), with=F])
}

# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#'
#'
#' details
#'
#' @param circs
#'
#' @export
annotateFlanks <- function(circs, txdb) {

  require(data.table)
  require(GenomicRanges)

  cat('Munging input data...\n')
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

  cat('Loading annotation...\n')
  annot.list <- GRangesList(utr5   = reduce(unlist(fiveUTRsByTranscript(txdb))),
                            utr3   = reduce(unlist(threeUTRsByTranscript(txdb))),
                            cds    = reduce(cds(txdb)),
                            intron = reduce(unlist(intronsByTranscript(txdb))))

  cat('Annotating circRNAs...\n')
  circ.starts.gr$feat_start     <- AnnotateRanges(r1 = circ.starts.gr, l = annot.list, type="precedence")
  circ.starts.gr$feat_start_all <- AnnotateRanges(r1 = circ.starts.gr, l = annot.list, type="all")
  circ.ends.gr$feat_end         <- AnnotateRanges(r1 = circ.ends.gr,   l = annot.list, type="precedence")
  circ.ends.gr$feat_end_all     <- AnnotateRanges(r1 = circ.ends.gr,   l = annot.list, type="all")

  cat('Merging data')
  circs$id <- paste(circs$chrom, ":", circs$start, "-", circs$end, sep="")
  circs <- merge(circs, data.table(as.data.frame(values(circ.starts.gr))), by="id")
  circs <- merge(circs, data.table(as.data.frame(values(circ.ends.gr))), by="id")

  return(circs)
}

# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#' annotates the ranges with the corresponding list
#'
#' details
#'
#' @param circs
#'
#' @export
AnnotateRanges = function(r1, l, ignore.strand=FALSE, type = 'precedence', null.fact = 'None', collapse.char=':') {

  if(! class(r1) == 'GRanges')
    stop('Ranges to be annotated need to be GRanges')

  if(! all(sapply(l, class) == 'GRanges'))
    stop('Annotating ranges need to be GRanges')

  if(!type %in% c('precedence','all'))
    stop('type may only be precedence and all')

  require(data.table)
  require(GenomicRanges)
  cat('Overlapping...\n')
  if(class(l) != 'GRangesList')
    l = GRangesList(lapply(l, function(x){values(x)=NULL;x}))
  a = suppressWarnings(data.table(as.matrix(findOverlaps(r1, l, ignore.strand=ignore.strand))))
  a$id = names(l)[a$subjectHits]
  a$precedence = match(a$id,names(l))[a$subjectHits]
  a = a[order(a$precedence)]

  if(type == 'precedence'){
    cat('precedence...\n')
    a = a[!duplicated(a$queryHits)]
    annot = rep(null.fact, length(r1))
    annot[a$queryHits] = a$id
  }
  if(type == 'all'){
    cat('all...\n')
    a = a[,list(id=paste(unique(id),collapse=collapse.char)),by='queryHits']
    annot = rep(null.fact, length(r1))
    annot[a$queryHits] = a$id

  }

  return(annot)
}

