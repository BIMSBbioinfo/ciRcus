
gtf2sqlite(assembly = "hg19", db.file="data/test.sqlite")

annot.list <- loadAnnotation("data/test.sqlite")

circs.f <- annotateCircs(circs.bed = "data/Sy5y_D0_sites.bed", annot.list = annot.list, assembly = "hg19")

circHist(circs.f, 0.5)
annotPie(circs.f, 0.02)
