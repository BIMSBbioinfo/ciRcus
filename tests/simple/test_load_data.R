all.sites <- readCircs(file = "/data/circrna/Mouse/P19_diff/undiff/undiff_sites.bed", qualfilter = FALSE)
all.sites.f <- readCircs(file = "/data/circrna/Mouse/P19_diff/undiff/undiff_sites.bed", qualfilter = TRUE, n_uniq_thr = 2, keepCols = 1:19)
circs.f <- circLinRatio(sites = all.sites.f[,.(chrom, start, end, name, n_reads, strand, n_uniq)])

source("/clusterhome/pglazar/bin/lib/little_helpers.R")


library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb  <- TxDb.Mmusculus.UCSC.mm10.ensGene
txdb2 <- loadDb("/data/BIO2/pglazar/indices/txdb.ensembl67.mmu.sqlite")

tx.1 <- transcriptsBy(txdb, by="gene")
tx.2 <- transcriptsBy(txdb2, by="gene")
