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
#annot.list <- loadAnnotation("data/test.sqlite")
annot.list <- loadAnnotation("inst/extdata/hsa_ens75_minimal.sqlite")
cdata <- data.frame(sample=c("FC1", "FC2", "H1", "H2", "L1", "L2"),
                    filename=dir("inst/extdata/", full.names=T)[grep("rep", dir("inst/extdata/"))])
#se <- summarizeCircs(dir("inst/extdata/", full.names=T)[grep("rep", dir("inst/extdata/"))], wobble=1, colData = cdata)
se <- summarizeCircs(colData = cdata, wobble=1)
se <- annotateCircs(se, annot.list, "mm9", fixCoordIndexing = TRUE)
tab <- resTable(se)
tab

se <- annotateHostGenes(se, annot.list$genes)
resTable(se)
se <- annotateFlanks(se, annot.list$gene.feats)
resTable(se)
se <- annotateJunctions(se, annot.list$junctions)
resTable(se)
se <- circLinRatio(se)
resTable(se)
histogram(se)
annotPie(se)

# find_circ2
annot.list <- loadAnnotation("data/test.sqlite")
#annot.list <- loadAnnotation("inst/extdata/hsa_ens75_minimal.sqlite")
#cdata <- data.frame(sample=c("FC1"),
#                    filename="../data/fc2/FrontalCortex_rep1_circ_splice_sites.bed")
# cdata <- data.frame(sample=c("CB1", "CB2", "FC1", "FC2"),
#                     filename=c("../data/fc2_full/Cerebellum/rep1/circ_splice_sites.bed",
#                                "../data/fc2_full/Cerebellum/rep2/circ_splice_sites.bed",
#                                "../data/fc2_full/FrontalCortex/rep1/circ_splice_sites.bed",
#                                "../data/fc2_full/FrontalCortex/rep2/circ_splice_sites.bed"))
cdata <- data.frame(sample=c("CB1", "CB2", "FC1", "FC2"),
                    filename=c("/data/rajewsky/home/pglazar/projects/ENCODE_brain/human/Cerebellum/rep1/circ_splice_sites.bed",
                               "/data/rajewsky/home/pglazar/projects/ENCODE_brain/human/Cerebellum/rep2/circ_splice_sites.bed",
                               "/data/rajewsky/home/pglazar/projects/ENCODE_brain/human/FrontalCortex/rep1/circ_splice_sites.bed",
                               "/data/rajewsky/home/pglazar/projects/ENCODE_brain/human/FrontalCortex/rep2/circ_splice_sites.bed"))
#se <- summarizeCircs(dir("inst/extdata/", full.names=T)[grep("rep", dir("inst/extdata/"))], wobble=1, colData = cdata)
se <- summarizeCircs(colData = cdata, wobble=1)
se <- annotateHostGenes(se, annot.list$genes)
resTable(se)
se <- annotateFlanks(se, annot.list$gene.feats)
resTable(se)
se <- annotateJunctions(se, annot.list$junctions)
resTable(se)
se <- circLinRatio(se)
resTable(se)
histogram(se)
annotPie(se)

se2 <- summarizeCircs(colData = cdata, wobble=1)
se2 <- annotateCircs(se2, annot.list, "hg19")
histogram(se2)
annotPie(se2, other.threshold = 0.02)

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
annotPie(se, 0.1)


# wrapper
cdata <- data.frame(sample=c("FC1", "FC2", "H1", "H2", "L1", "L2"),
                    filename=dir("data/demo/", full.names=T))
se  <- summarizeCircs(circ.files = as.character(cdata$filename), wobble=1, colData = cdata)
se1 <- annotateCircs(se, annot.list=annot.list)
cdata <- data.frame(sample=c("D0"),
                    filename="data/Sy5y_D0_sites.bed")
se  <- summarizeCircs(circ.files = as.character(cdata$filename), wobble=1, colData = cdata)
se2 <- annotateCircs(se = se, annot.list=annot.list, fixCoordIndexing = TRUE)

# connect to the public instance of circBase
con <- dbConnect(drv    = dbDriver("MySQL"),
                 host   = "141.80.181.75",
                 user   = "circbase",
                 pass   = "circbase",
                 dbname = "circbase",
                 port   = 3306)


ah <- AnnotationHub()
gtf.gr <- ah[[getOption("assembly2annhub")[["hg19"]]]]
gtf.gr <- keepStandardChromosomes(gtf.gr)
seqlevels(gtf.gr) <- paste("chr", seqlevels(gtf.gr), sep="")
seqlevels(gtf.gr)[which(seqlevels(gtf.gr) == "chrMT")] <- "chrM"
gtf.gr <- subsetByOverlaps(gtf.gr, rowRanges(se))
txdb <- makeTxDbFromGRanges(gtf.gr, drop.stop.codons = FALSE, metadata = data.frame(name="Genome", value="GRCh37"))
saveDb(txdb, file="../data/hsa_ens75_minimal.gtf")


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
