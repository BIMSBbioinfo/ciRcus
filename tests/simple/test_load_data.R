
annot.list <- loadAnnotation("data/hsa_GRCh37_ens75.sqlite")

circs.f <- annotateCircs(circs.bed = "data/Sy5y_D0_sites.bed", annot.list = annot.list, assembly = "hg19")

circHist(circs.f, 0.5)
annotPie(circs.f, 0.02)
