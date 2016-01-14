
gtf2sqlite(assembly = "hg19", db.file="data/test.sqlite")

annot.list <- loadAnnotation("data/test.sqlite")

circs.f <- annotateCircs(circs.bed = "data/Sy5y_D0_sites.bed", annot.list = annot.list, assembly = "hg19")

a <- summarizeCircs(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))])

a <- summarizeCircs(c("data/demo//FrontalCortex_rep1_sites.bed", "data/demo//FrontalCortex_rep2_sites.bed", "data/demo//HeartF_rep1_sites.bed"), wobble = 5)



circHist(circs.f, 0.5)
annotPie(circs.f, 0.02)
