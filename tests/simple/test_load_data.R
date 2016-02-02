
gtf2sqlite(assembly = "hg19", db.file="data/test.sqlite")

annot.list <- loadAnnotation("data/test.sqlite")

circs.f <- annotateCircs(circs.bed = "data/Sy5y_D0_sites.bed", annot.list = annot.list, assembly = "hg19")



circHist(circs.f, 0.5)
annotPie(circs.f, 0.02)

library(data.table)
se <- summarizeCircs(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))], wobble=1)

circs <- lapply(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))], readCircs)
circs <- lapply(circs, circLinRatio, return.readcounts=TRUE)


# circ.files <- dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))]
# keep.linear = TRUE
# wobble = 1
# subs = "all"
# qualfilter = TRUE
# keepCols = 1:6
# colData = NULL
