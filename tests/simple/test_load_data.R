
txdb.ens75 <- loadAnnotation(txdb.file = "data/hsa_GRCh37_ens75.sqlite")
exns <- unique(exons(txdb.ens75))
junct.start <- exns
#   end(junct.start)   <- start(junct.start) + 2
#   start(junct.start) <- start(junct.start) - 2
junct.end <- exns
#   start(junct.end) <- end(junct.end) - 2
#   end(junct.end)  <- end(junct.end) + 2

gene.feats <- GRangesList(utr5   = reduce(unlist(fiveUTRsByTranscript(txdb.ens75))),
                          utr3   = reduce(unlist(threeUTRsByTranscript(txdb.ens75))),
                          cds    = reduce(cds(txdb.ens75)),
                          intron = reduce(unlist(intronsByTranscript(txdb.ens75))))
junctions <- GRangesList(start = junct.start,
                         end   = junct.end)


circs.f <- readCircs(file = "data/Sy5y_D0_sites.bed")
circs.f <- circLinRatio(sites = circs.f)
circs.f <- getIDs(circs.f, "hsa", "hg19")
circs.f$start <- circs.f$start + 1
circs.f <- annotateHostGenes(circs = circs.f, txdb = txdb.ens75)
circs.f <- annotateFlanks(circs = circs.f, annot.list = gene.feats)
circs.f <- annotateJunctions(circs = circs.f, annot.list = junctions)
circs.f$gene <- ensg2name(circs.f$host, release = "75", organism = "hsa")



circHist(circs.f, 0.5)
annotPie(circs.f, 0.02)
