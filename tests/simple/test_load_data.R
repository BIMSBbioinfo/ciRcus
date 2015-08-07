#all.sites <- readCircs(file = "/data/circrna/Mouse/P19_diff/undiff/undiff_sites.bed", qualfilter = FALSE)
#all.sites.f <- readCircs(file = "/data/circrna/Mouse/P19_diff/undiff/undiff_sites.bed", qualfilter = TRUE, n_uniq_thr = 2, keepCols = 1:19)
#circs.f <- circLinRatio(sites = all.sites.f[,.(chrom, start, end, name, n_reads, strand, n_uniq)])


all.sites.f <- readCircs(file = "/data/circrna/Human/Sy5y_diff/D0/Sy5y_D0_sites.bed", qualfilter = TRUE, n_uniq_thr = 2, keepCols = 1:19)
circs.f <- circLinRatio(sites = all.sites.f[,.(chrom, start, end, name, n_reads, strand, n_uniq)])

ptm <- proc.time()
circs.fa <- annotateCircs(circs.f)
proc.time() - ptm


#source("/clusterhome/pglazar/bin/lib/little_helpers.R")


# library(TxDb.Mmusculus.UCSC.mm10.ensGene)
# txdb  <- TxDb.Mmusculus.UCSC.mm10.ensGene
#
# txdb2 <- loadDb("/data/BIO2/pglazar/indices/txdb.ensembl67.mmu.sqlite")
#
# tx.1 <- transcriptsBy(txdb, by="gene")
# tx.2 <- transcriptsBy(txdb2, by="gene")
#
#
# circs <- circs.f
# read.counts <- readCountsToDT()
#
# function(circs, read.counts, biotype, genelengths, exonsByGene, return.subset="list", skip.genes=NA) {
#
#   # get gene expression data
#   tpms <- readCountsToDT(filename=read.counts, biotype=biotype, genelengths=genelengths, skip.genes=skip.genes)
#
#   # organize a GRanges object with all the exons
#   exons.gr <- unlist(exonsByGene, use.names=TRUE)
#   exons.gr <- exons.gr[names(exons.gr) %in% tpms$gene]
#   seqlevels(exons.gr) <- paste("chr", seqlevels(exons.gr), sep="")
#   exons.gr <- keepSeqlevels(exons.gr, paste("chr", c(1:22, "X", "Y"), sep=""))
#   start(exons.gr) <- start(exons.gr) - 1
#   elementMetadata(exons.gr) <- NULL
#   exons.gr <- unique(exons.gr)
#
#   # organize circs
#   circs.gr <- GRanges(seqnames=circs$chrom,
#                       ranges=IRanges(start=circs$start,
#                                      end=circs$end),
#                       strand=circs$strand,
#                       id=circs$id)
#   circs.gr <- sort(circs.gr)
#   circs.gr.exonic <- circs.gr[overlapsAny(circs.gr, exons.gr, type="start") | overlapsAny(circs.gr, exons.gr, type="end")]
#   circs.gr.other <- circs.gr[!(circs.gr$id %in% circs.gr.exonic$id)]
#   # kei, now i have a GR with all circs that did and did not match annotated exons
#
#   olaps.start <- findOverlaps(circs.gr.exonic, exons.gr, type="start")
#   olaps.end   <- findOverlaps(circs.gr.exonic, exons.gr, type="end")
#
#   matches.start <- data.table(id=circs.gr.exonic$id[queryHits(olaps.start)], gene=names(exons.gr)[subjectHits(olaps.start)] )
#   matches.end   <- data.table(id=circs.gr.exonic$id[queryHits(olaps.end)],   gene=names(exons.gr)[subjectHits(olaps.end)] )
#   matches <- rbind(matches.start, matches.end)
#   matches <- matches[!duplicated(matches)]
#
#   onetomany <- matches[matches$id %in% names(table(matches$id)[table(matches$id)>1])]
#   onetomany <- merge(onetomany, tpms[,c(1,2,5), with=FALSE], by="gene", all.x=TRUE)
#   matches <- matches[!(matches$id %in% onetomany$id),]
#   matches <- merge(matches, tpms[,c(1,2,5), with=FALSE], by="gene", all.x=TRUE)
#   # implement a logic for clearing up oneToMany relationships
#   # (picking the gene with higer readcount, longer overlap with a circ, or something)
#   circs.matches <- merge(circs, matches, by="id")
#   circs.matches[ratio == Inf, ratio := as.numeric(n_reads)]
#   circs.ambiguous <- merge(circs, onetomany, by="id")
#   circs.other <- circs[!(as.character(circs$name) %in% unique(c(as.character(circs.matches$name), as.character(circs.ambiguous$name)))),]
