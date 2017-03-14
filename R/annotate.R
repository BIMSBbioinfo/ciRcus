# ---------------------------------------------------------------------------- #
#' Build and save a TxDb object containing gene annotation
#'
#' Loads the GTF file for the selected assembly, drops non-standard chromosomes,
#' adds "chr" prefix to Ensembl chromosome names and saves the SQLite database
#' to \code{db.file}
#'
#' @param assembly abbreviation for one of the supported assemblies
#' @param db.file a file to save SQLite database to
#'
#' @export
gtf2sqlite <-
  function(assembly = c("hg19", "hg38", "mm10", "rn5", "dm6", "WBcel235"),
           db.file) {

  ah <- AnnotationHub()
  gtf.gr <- ah[[getOption("assembly2annhub")[[assembly]]]]
  gtf.gr <- keepStandardChromosomes(gtf.gr)
  seqlevels(gtf.gr) <- paste("chr", seqlevels(gtf.gr), sep = "")
  seqlevels(gtf.gr)[which(seqlevels(gtf.gr) == "chrMT")] <- "chrM"
  txdb <- makeTxDbFromGRanges(gtf.gr, drop.stop.codons = FALSE,
                              metadata = data.frame(name = "Genome",
                                                    value = assembly))
  saveDb(txdb, file = db.file)

}


# ---------------------------------------------------------------------------- #
#' Load annotation and prepare a list of features for later
#'
#' Loads a local .sqlite TxDb annotation file (e.g. created using
#' \code{gtf2sqlite}), returns a list of features needed for circRNA annotation.
#'
#' @param txdb.file path to the TxDb gene annotation file saved as SQLite
#' database
#'
#' @export
loadAnnotation <- function(txdb.file) {

  txdb <- loadDb(txdb.file)
  exns <- unique(exons(txdb))
  junct.start <- exns
  end(junct.start) <- start(junct.start)
  junct.end <- exns
  start(junct.end) <- end(junct.end)

  gene.feats <-
    GRangesList(cds    = reduce(cds(txdb)),
                utr5   = reduce(unlist(fiveUTRsByTranscript(txdb))),
                utr3   = reduce(unlist(threeUTRsByTranscript(txdb))),
                intron = reduce(unlist(intronsByTranscript(txdb))),
                tx     = reduce(exns))
  junctions <- GRangesList(start = junct.start,
                           end   = junct.end)

  genes <- genes(txdb)

  return(list(genes = genes, gene.feats = gene.feats, junctions = junctions))
}

# ---------------------------------------------------------------------------- #
#' Load and annotate a list of circRNA candidates
#'
#' Loads a list of splice junctions detected using \code{find_circ.py}
#' (Memczak et al. 2013; www.circbase.org), applies quality filters,
#' calculates circular-to-linear ratios, and extends the input with
#' genomic features
#'
#' @param se a SummarizedExperiment object
#' @param annot.list list of relevant genomic features generated using
#' \code{loadAnnotation()}
#' @param assembly what genome assembly the input data are coming from
#' @param fixCoordIndexing check and fix genomic coordinate indexing?
#' @param ... other arguments
#'
#' @docType methods
#' @rdname annotateCircs-methods
#' @export
setGeneric("annotateCircs",
           function(se,
                    annot.list,
                    assembly = c("hg19", "hg38", "mm10", "rn5", "dm6",
                                 "WBcel235"),
                    fixCoordIndexing = TRUE,
                    ...)
             standardGeneric("annotateCircs"))

#' @aliases annotateCircs, RangedSummarizedExperiment-method
#' @rdname annotateCircs-methods
setMethod("annotateCircs", signature("RangedSummarizedExperiment"),
          function(se, annot.list,
                   assembly = c("hg19", "hg38", "mm10", "rn5", "dm6",
                                "WBcel235"),
                   fixCoordIndexing = TRUE, ...) {

            if (fixCoordIndexing == TRUE) {
              message("checking out coordinate indexing...")
              coordfix <- testCoordinateIndexing(rowRanges(se),
                                                 annot.list$gene.feats$cds)
              tophits <- sapply(coordfix, function(x) which(x / sum(x) > 0.9))
              if (!is.na(max(coordfix[[1]][2:3] / sum(coordfix[[1]])) > 0.9) &
                  max(coordfix[[1]][2:3] / sum(coordfix[[1]])) > 0.9) {
                if (tophits[1] == 2) {
                  start(se) <- start(se) + 1
                } else if (tophits[1] == 3) {
                  start(se) <- start(se) - 1
                }
                warning("start coordinates modified to match annotation.")
              }

              if (!is.na(max(coordfix[[2]][2:3] / sum(coordfix[[2]])) > 0.9) &
                  max(coordfix[[2]][2:3] / sum(coordfix[[2]])) > 0.9) {
                if (tophits[2] == 2) {
                  end(se) <- end(se) + 1
                } else if (tophits[2] == 3) {
                  end(se) <- end(se) - 1
                }
                warning("end coordinates modified to match annotation.")
              }
            } else {
              testCoordinateIndexing(rowRanges(se), annot.list$gene.feats$cds)
              warning("input coordinates were not modified.")
            }

            # check seqlevel style, match input
            if (!any(seqlevelsStyle(annot.list$genes) %in%
                     seqlevelsStyle(se))) {
              warning("nonmatching seqlevel styles, will fix automatically")
              for (i in 1:length(annot.list)) {
                seqlevelsStyle(annot.list[[i]]) <- seqlevelsStyle(se)
              }
            }

            message("annotating host genes...")
            se <- annotateHostGenes(se, annot.list$genes)
            message("annotating splice junctions...")
            se <- annotateFlanks(se, annot.list$gene.feats)
            se <- annotateJunctions(se, annot.list$junctions)

            if (grepl("linear", names(assays(se)))) {
              message("calculating circular/linear ratios...")
              se <- circLinRatio(se)
            } else {
              message(paste("no linear splicing info found, skipping",
                            "circular/linear ratios..."))
            }

            return(se)

          })
# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#'
#' details
#'
#' @param se a SummarizedExperiment object
#' @param genes.gr a GenomicRanges object with gene models
#'
#' @export
annotateHostGenes <- function(se, genes.gr) {

  # circs to GR
  circs.gr <- rowRanges(se)
  circs.gr$id <- paste0(seqnames(circs.gr), ":", start(circs.gr), "-",
                        end(circs.gr))


  # GR with circ starts and ends only (left and right flanks, actually)
  circ.starts.gr <- circs.gr
  end(circ.starts.gr) <- start(circ.starts.gr)
  circ.ends.gr <- circs.gr
  start(circ.ends.gr) <- end(circ.ends.gr)

  olap.start <- findOverlaps(circ.starts.gr, genes.gr, type = "within")
  olap.end   <- findOverlaps(circ.ends.gr, genes.gr, type = "within")

  circs <- data.table(start.hit = names(circs.gr) %in% queryHits(olap.start),
                      end.hit   = names(circs.gr) %in% queryHits(olap.end),
                      id        = circs.gr$id,
                      ord       = 1:length(circs.gr))

  matches.start <- data.table(id = circs.gr$id[queryHits(olap.start)],
                              gene = names(genes.gr)[subjectHits(olap.start)])
  matches.end   <- data.table(id = circs.gr$id[queryHits(olap.end)],
                              gene = names(genes.gr)[subjectHits(olap.end)])

  start.list <- lapply(split(matches.start, matches.start$id),
                       function(x) x$gene)
  end.list   <- lapply(split(matches.end,   matches.end$id),
                       function(x) x$gene)

  circs <- merge(circs,
                 data.table(id = names(start.list),
                            starts = sapply(start.list,
                                            function(x)
                                              paste(x, collapse = ","))),
                 by = "id", all.x = T)
  circs <- merge(circs,
                 data.table(id = names(end.list),
                            ends = sapply(end.list,
                                          function(x)
                                            paste(x, collapse = ","))),
                 by = "id", all.x = T)

  hs <- hash(start.list)
  he <- hash(end.list)

  circs$hit.ctrl <- circs$id %in% unique(c(keys(hs), keys(he)))

  tmphits <- integer()
  tmpgenes <- integer()
  host.candidates <- integer()
  for (circ in circs$id) {
    tmphits <- append(tmphits, sum(hs[[circ]] %in% he[[circ]]))
    tmpgenes <- append(tmpgenes, paste(hs[[circ]][hs[[circ]] %in% he[[circ]]],
                                       collapse = ","))
    host.candidates <- append(host.candidates, length(unique(c(hs[[circ]],
                                                               he[[circ]]))))
  }

  circs$hitcnt <- tmphits
  circs$hitgenes <- tmpgenes
  circs$host.candidates <- host.candidates

  circs$host[circs$hitcnt == 1] <- circs$hitgenes[circs$hitcnt == 1]
  circs$host[circs$hitcnt > 1]  <- "ambiguous"
  # TODO: actually, some of them may have a putative host gene within, I was
  # only testing starts / ends
  circs$host[circs$hitcnt == 0 & circs$start.hit == FALSE &
             circs$end.hit == FALSE] <- "intergenic"
  circs$host[circs$hitcnt == 0 & circs$start.hit == TRUE  &
             circs$end.hit == TRUE]  <- "no_single_host"
  circs$host[circs$hitcnt == 0 & xor(circs$start.hit, circs$end.hit) &
             circs$host.candidates > 1] <- "ambiguous"
  circs$host[circs$hitcnt == 0 & circs$start.hit == TRUE  &
             circs$end.hit == FALSE & circs$host.candidates == 1] <-
    circs$starts[circs$hitcnt == 0 & circs$start.hit == TRUE   &
                 circs$end.hit == FALSE & circs$host.candidates == 1]
  circs$host[circs$hitcnt == 0 & circs$start.hit == FALSE &
             circs$end.hit == TRUE & circs$host.candidates == 1]  <-
    circs$ends[circs$hitcnt == 0   & circs$start.hit == FALSE  &
               circs$end.hit == TRUE  & circs$host.candidates == 1]

  rowRanges(se)$gene_id <- circs$host[order(circs$ord)]

  return(se)
}


# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#'
#'
#' details
#'
#' @param se a SummarizedExperiment object
#' @param annot.list an annotation list, as returned by \code{loadAnnotation()}
#'
#' @export
annotateFlanks <- function(se, annot.list) {

  circs.gr <- rowRanges(se)

  # GR with circ starts and ends only (left and right flanks, actually)
  circ.starts.gr <- circs.gr
  end(circ.starts.gr) <- start(circ.starts.gr)
  circ.ends.gr <- circs.gr
  start(circ.ends.gr) <- end(circ.ends.gr)

  circ.starts.gr$feat_start <- AnnotateRanges(r1 = circ.starts.gr,
                                              l = annot.list,
                                              null.fact = "intergenic",
                                              type = "precedence")
  circ.ends.gr$feat_end     <- AnnotateRanges(r1 = circ.ends.gr,
                                              l = annot.list,
                                              null.fact = "intergenic",
                                              type = "precedence")


  rowRanges(se)$feature <- ifelse(as.character(strand(rowRanges(se))) == "+",
                                  paste(circ.starts.gr$feat_start,
                                        circ.ends.gr$feat_end,     sep = ":"),
                                  paste(circ.ends.gr$feat_end,
                                        circ.starts.gr$feat_start, sep = ":"))


  return(se)
}



# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#'
#'
#' details
#'
#' @param se a SummarizedExperiment object
#' @param annot.list an annotation list, as returned by \code{loadAnnotation()}
#'
#' @export
annotateJunctions <- function(se, annot.list) {

  circs.gr <- rowRanges(se)

  # GR with circ starts and ends only (left and right flanks, actually)
  circ.starts.gr <- circs.gr
  end(circ.starts.gr) <- start(circ.starts.gr)
  circ.ends.gr <- circs.gr
  start(circ.ends.gr) <- end(circ.ends.gr)

  circ.starts.gr$annotated_start_junction <- AnnotateRanges(r1 = circ.starts.gr,
                                                            l = annot.list,
                                                            type = "precedence")
  circ.ends.gr$annotated_end_junction     <- AnnotateRanges(r1 = circ.ends.gr,
                                                            l = annot.list,
                                                            type = "precedence")


  circ.starts.gr$annotated_start_junction <-
    ifelse(circ.starts.gr$annotated_start_junction == "None", FALSE, TRUE)
  circ.ends.gr$annotated_end_junction     <-
    ifelse(circ.ends.gr$annotated_end_junction == "None", FALSE, TRUE)

  junct.known <- rep("", length(circs.gr))
  junct.known[circ.starts.gr$annotated_start_junction == TRUE  &
              circ.ends.gr$annotated_end_junction == TRUE] <- "both"
  junct.known[circ.starts.gr$annotated_start_junction == FALSE &
              circ.ends.gr$annotated_end_junction == FALSE] <- "none"
  junct.known[circ.starts.gr$annotated_start_junction == TRUE  &
              circ.ends.gr$annotated_end_junction == FALSE &
              as.logical(strand(circs.gr) == "+")] <- "5pr"
  junct.known[circ.starts.gr$annotated_start_junction == TRUE  &
              circ.ends.gr$annotated_end_junction == FALSE &
              as.logical(strand(circs.gr) == "-")] <- "3pr"
  junct.known[circ.starts.gr$annotated_start_junction == FALSE &
              circ.ends.gr$annotated_end_junction == TRUE  &
              as.logical(strand(circs.gr) == "+")] <- "3pr"
  junct.known[circ.starts.gr$annotated_start_junction == FALSE &
              circ.ends.gr$annotated_end_junction == TRUE  &
              as.logical(strand(circs.gr) == "-")] <- "5pr"

  rowRanges(se)$junct.known <- junct.known

  return(se)
}

# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#' annotates the ranges with the corresponding list
#'
#' details
#'
#' @param r1 a GenomicRanges object with ranges to be annotated
#' @param l a GRangesList with genomic features
#' @param ignore.strand should strand info be ignored
#' @param type annotation mode, \code{precedence} or \code{all}
#' @param null.fact a value to be used for ranges that do not hit any annotation
#' @param collapse.char a collapse character for multiple hits in \code{all} mode
#'
#' @export
AnnotateRanges <- function(r1, l, ignore.strand = FALSE, type = "precedence",
                           null.fact = "None", collapse.char = ":") {

  if (! class(r1) == "GRanges")
    stop("Ranges to be annotated need to be GRanges")

  if (! all(sapply(l, class) == "GRanges"))
    stop("Annotating ranges need to be GRanges")

  if (!type %in% c("precedence", "all"))
    stop("type may only be precedence and all")

  if (class(l) != "GRangesList")
    l <- GRangesList(lapply(l, function(x){
                                 values(x) <- NULL;x
                               }))
  a <- suppressWarnings(data.table(as.matrix(findOverlaps(r1, l,
                                                          ignore.strand =
                                                            ignore.strand))))
  a$id <- names(l)[a$subjectHits]
  a$precedence <- match(a$id, names(l))
  a <- a[order(a$precedence)]

  if (type == "precedence"){
    a <- a[!duplicated(a$queryHits)]

  }

  if (type == "all"){
    a <- a[, list(id = paste(unique(id), collapse = collapse.char)),
           by = "queryHits"]
  }

  annot <- rep(null.fact, length(r1))
  annot[a$queryHits] <- a$id

  return(annot)
}

# ---------------------------------------------------------------------------- #
#' title
#'
#' description
#'
#'
#' details
#'
#' @param ensg character vector of ensembl gene ids
#' @param organism three-letter organism abbreviation, such as hsa, mmu, rno...
#' @param release which ensembl release to use?
#'
ensg2name <- function(ensg, organism, release = "current") {

  ensembl.host <- getOption("ensembl.release")[[release]]

  ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = ensembl.host)
  ensembl <-
    useDataset(dataset = paste(getOption("ensembl.organism")[[organism]],
                               "_gene_ensembl", sep = ""),
               mart = ensembl)

  xrefs <- getBM(attributes = c("external_gene_id", "ensembl_gene_id"),
                 filter     = "ensembl_gene_id",
                 values     = ensg,
                 mart       = ensembl)

  out.dt <- data.table(ensembl_gene_id = ensg)
  out.dt <- merge(out.dt, xrefs, by = "ensembl_gene_id", all.x = T)
  out.dt <- out.dt[match(ensg, out.dt$ensembl_gene_id)]

  return(out.dt$external_gene_id)
}
