
all.sites.f <- readCircs(file = "/data/circrna/Human/Sy5y_diff/D0/Sy5y_D0_sites.bed", qualfilter = TRUE, n_uniq_thr = 2, keepCols = 1:19)
circs.f <- circLinRatio(sites = all.sites.f[,.(chrom, start, end, name, n_reads, strand, n_uniq)])
circs.f <- getIDs(circs.f, "hsa", "hg19")

system.time(
  txdb.ens75 <- loadAnnotation(txdb.file = "data/hsa_GRCh37_ens75.sqlite")
)

circs.f$start <- circs.f$start + 1

system.time(
  circs.fa <- annotateHostGenes(circs = circs.f, txdb = txdb.ens75)
)

ptm <- proc.time()
  exns <- unique(exons(txdb.ens75))
  junct.start <- exns
  end(junct.start)   <- start(junct.start) + 2
  start(junct.start) <- start(junct.start) - 2
  junct.end <- exns
  start(junct.end) <- end(junct.end) - 2
  end(junct.end)  <- end(junct.end) + 2

  gene.feats <- GRangesList(utr5   = reduce(unlist(fiveUTRsByTranscript(txdb.ens75))),
                            utr3   = reduce(unlist(threeUTRsByTranscript(txdb.ens75))),
                            cds    = reduce(cds(txdb.ens75)),
                            intron = reduce(unlist(intronsByTranscript(txdb.ens75))))
  junctions <- GRangesList(start = junct.start,
                           end   = junct.end)
proc.time() - ptm

system.time(
  circs.faf <- annotateFlanks(circs = circs.fa, annot.list = gene.feats)
)

system.time(
  circs.fafj <- annotateJunctions(circs = circs.fa, annot.list = junctions)
)

circs.fafj$name <- ensg2name(circs.fafj$host, release = "75", organism = "hsa")
