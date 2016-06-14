# master
gtf2sqlite(assembly = "hg19", db.file="data/test.sqlite")

annot.list <- loadAnnotation("data/test.sqlite")

circs.f <- annotateCircs(circs.bed = "data/Sy5y_D0_sites.bed", annot.list = annot.list, assembly = "hg19")



circHist(circs.f, 0.5)
annotPie(circs.f, 0.02)

# development
#library(data.table)
#library(GenomicRanges)
#library(hash)
annot.list <- loadAnnotation("data/test.sqlite")
cdata <- data.frame(sample=c("FC1", "FC2", "H1", "H2", "L1", "L2"),
                    filename=basename(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))]))
se <- summarizeCircs(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))], wobble=1, colData = cdata)
se <- annotateHostGenes(se, annot.list$genes)
resTable(se)
se <- annotateFlanks(se, annot.list$gene.feats)
resTable(se)
se <- annotateJunctions(se, annot.list$junctions)
resTable(se)
se <- circLinRatio(se)
resTable(se)
histogram(se)

# SH-SY5Y test
cdata <- data.frame(sample=c("D0"),
                    filename="data/Sy5y_D0_sites.bed")
se <- summarizeCircs(as.character(cdata$filename), wobble=1, colData = cdata)
se <- annotateHostGenes(se, annot.list$genes)
se <- annotateFlanks(se, annot.list$gene.feats)
se <- annotateJunctions(se, annot.list$junctions)
se <- circLinRatio(se)
DT <- resTable(se)
histogram(se)


circ.ends.gr <- rowRanges(se)
start(circ.ends.gr) <- end(circ.ends.gr)

AnnotateRanges(r1 = circ.ends.gr,         l = annot.list$gene.feats,  null.fact = "intergenic", type="precedence")[1547]
AnnotateRanges(r1 = circ.ends.gr[1547],   l = annot.list$gene.feats,  null.fact = "intergenic", type="precedence")



circs <- lapply(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))], readCircs)
circs <- lapply(circs, circLinRatio, return.readcounts=TRUE)


circ.files <- dir("/clusterhome/pglazar/projects/circus/ciRcus/data/demo/", full.names=T)[grep("sites", dir("/clusterhome/pglazar/projects/circus/ciRcus/data/demo/"))]
keep.linear = TRUE
wobble = 1
subs = "all"
qualfilter = TRUE
keepCols = 1:6
colData = NULL
