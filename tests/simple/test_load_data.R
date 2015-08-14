
all.sites.f <- readCircs(file = "/data/circrna/Human/Sy5y_diff/D0/Sy5y_D0_sites.bed", qualfilter = TRUE, n_uniq_thr = 2, keepCols = 1:19)
circs.f <- circLinRatio(sites = all.sites.f[,.(chrom, start, end, name, n_reads, strand, n_uniq)])

system.time(
  txdb.ens75 <- loadAnnotation(txdb.file = "data/hsa_GRCh37_ens75.sqlite")
)

system.time(
  circs.fa <- annotateCircs(circs = circs.f, txdb = txdb.ens75)
)



