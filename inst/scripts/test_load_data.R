#CIRI2
annot.list <- loadAnnotation("inst/extdata/db/test.sqlite")
cdata <- data.frame(sample = c("riboz", "RNaseR"),
                    filename = paste0("inst/extdata/ciri_demo_hek/",
                                      c("HEK_riborezo_CIRI_sites.txt",
                                        "HEK_RNaseR_CIRI_sites.txt")))

# colData <- cdata
# keep.linear <- TRUE
# wobble <- 1
# subs <- "all"
# qualfilter <- FALSE
# keepCols <- 1:12


circs.se <- summarizeCircs(colData = cdata, keep.linear = FALSE, wobble = 1,
                           subs = "all", qualfilter = FALSE, keepCols = 1:12)
circs.se <-  annotateCircs(se = circs.se, annot.list = annot.list,
                           assembly = "hg19")

# issue16
tmp <- resTable(circs.se)
tmp[feature == "intergenic:intergenic"  & gene_id != "intergenic"]
tmp.gr <- GRanges(seqnames = tmp[feature == "intergenic:intergenic"  & gene_id != "intergenic"][3:7,]$chr,
                  ranges   = IRanges(start=tmp[feature == "intergenic:intergenic"  & gene_id != "intergenic"][3:7]$start,
                                     end=tmp[feature == "intergenic:intergenic"  & gene_id != "intergenic"][3:7]$end),
                  strand = "-")

tmp.gr <- range(tmp.gr)


# find_circ2
annot.list <- loadAnnotation("inst/extdata/db/test.sqlite")

cdata <- data.frame(sample = c("D0", "D2", "D4"),
                    filename = c("inst/extdata/Sy5y_D0_sites.bed",
                                 paste0("/data/circrna/Human/Sy5y_diff/",
                                        c("D2/Sy5y_D2_sites.bed",
                                          "D4/Sy5y_D4_sites.bed"))))

circs.se <- summarizeCircs(colData = cdata, wobble = 1, keepCols = 1:19)

circs.se <- annotateCircs(se = circs.se, annot.list = annot.list,
                          assembly = "hg19", fixCoordIndexing = T)

# issue #16

histogram(circs.se, 0.5)
annotPie(circs.se, 0.02)
uniqReadsQC(circs.se, "all")

## issue #16 #
tmp <- resTable(circs.se)
tmp[feature=="intergenic:intergenic" & gene_id == "no_single_host"]
rowRanges(circs.se)[start(rowRanges(circs.se))==119107780]
rowRanges(circs.se)[start(rowRanges(circs.se))==180669210]

rowRanges(circs.se)[start(rowRanges(circs.se)) %in% c(119107780, 180669210)]

sub.se <- circs.se[c(which(start(rowRanges(circs.se)) %in% c(119107779, 180669209)), 200:300),]





# development
#library(data.table)
#library(GenomicRanges)
#library(hash)
#annot.list <- loadAnnotation("data/test.sqlite")
annot.list <- loadAnnotation("inst/extdata/db/hsa_ens75_minimal.sqlite")
cdata <- data.frame(sample = c("FC1", "FC2", "H1", "H2", "L1", "L2"),
                    filename = dir("inst/extdata/encode_demo_small/",
                                   full.names = T))
#se <- summarizeCircs(dir("inst/extdata/",
#                         full.names = T)[grep("rep", dir("inst/extdata/"))],
#                     wobble = 1, colData = cdata)
se <- summarizeCircs(colData = cdata, wobble = 1, keepCols = 1:7)
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
#cdata <-
#  data.frame(sample = c("FC1"),
#             filename = "../data/fc2/FrontalCortex_rep1_circ_splice_sites.bed")
cdata <- data.frame(sample = c("CB1", "CB2", "FC1", "FC2"),
                    filename = paste0("/data/rajewsky/home/pglazar/projects/",
                                      "ENCODE_brain/human/",
#                    filename = paste0("../data/fc2_full/",
                                      c(paste0(   "Cerebellum/rep", 1:2),
                                        paste0("FrontalCortex/rep", 1:2)),
                                      "/circ_splice_sites.bed"))
#se <- summarizeCircs(dir("inst/extdata/",
#                         full.names = T)[grep("rep", dir("inst/extdata/"))],
#                     wobble = 1, colData = cdata)
se <- summarizeCircs(colData = cdata, wobble = 1)
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

se2 <- summarizeCircs(colData = cdata, wobble = 1)
se2 <- annotateCircs(se2, annot.list, "hg19")
histogram(se2)
annotPie(se2, other.threshold = 0.02)

# SH-SY5Y test
cdata <- data.frame(sample = c("D0"),
                    filename = "data/Sy5y_D0_sites.bed")
se <- summarizeCircs(as.character(cdata$filename), wobble = 1, colData = cdata)
se <- annotateHostGenes(se, annot.list$genes)
se <- annotateFlanks(se, annot.list$gene.feats)
se <- annotateJunctions(se, annot.list$junctions)
se <- circLinRatio(se)
DT <- resTable(se)
histogram(se)
annotPie(se, 0.1)


# wrapper
cdata <- data.frame(sample = c("FC1", "FC2", "H1", "H2", "L1", "L2"),
                    filename = dir("data/demo/", full.names = T))
se  <- summarizeCircs(circ.files = as.character(cdata$filename), wobble = 1,
                      colData = cdata)
se1 <- annotateCircs(se, annot.list = annot.list)
cdata <- data.frame(sample = c("D0"),
                    filename = "data/Sy5y_D0_sites.bed")
se  <- summarizeCircs(circ.files = as.character(cdata$filename), wobble = 1,
                      colData = cdata)
se2 <- annotateCircs(se = se, annot.list = annot.list, fixCoordIndexing = TRUE)

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
seqlevels(gtf.gr) <- paste("chr", seqlevels(gtf.gr), sep = "")
seqlevels(gtf.gr)[which(seqlevels(gtf.gr) == "chrMT")] <- "chrM"
gtf.gr <- subsetByOverlaps(gtf.gr, rowRanges(se))
txdb <- makeTxDbFromGRanges(gtf.gr, drop.stop.codons = FALSE,
                            metadata = data.frame(name = "Genome",
                                                  value = "GRCh37"))
saveDb(txdb, file = "../data/hsa_ens75_minimal.gtf")


circ.ends.gr <- rowRanges(se)
start(circ.ends.gr) <- end(circ.ends.gr)

AnnotateRanges(r1 = circ.ends.gr,       l = annot.list$gene.feats,
               null.fact = "intergenic", type = "precedence")[1547]
AnnotateRanges(r1 = circ.ends.gr[1547], l = annot.list$gene.feats,
               null.fact = "intergenic", type = "precedence")



circs <- lapply(dir("data/demo/", full.names = T)[grep("sites",
                                                       dir("data/demo/"))],
                readCircs)
circs <- lapply(circs, circLinRatio, return.readcounts = TRUE)


circ.files <- dir("/clusterhome/pglazar/projects/circus/ciRcus/data/demo/",
                  full.names <- T)[grep("sites",
                                   dir(paste0("/clusterhome/pglazar/projects/",
                                              "circus/ciRcus/data/demo/")))]
keep.linear <- TRUE
wobble <- 1
subs <- "all"
qualfilter <- TRUE
keepCols <- 1:6
colData <- NULL


# bedran
annot.list <- loadAnnotation("inst/extdata/db/hsa_ens75_minimal.sqlite")
cdata <- data.frame(sample = c("s1", "s2"),
                    filename = dir("../data/vedran/",
                                   pattern = "circ_splice_sites.bed$",
                                   full.names = T,
                                   recursive = T))
#se <- summarizeCircs(dir("inst/extdata/",
#                         full.names = T)[grep("rep", dir("inst/extdata/"))],
#                     wobble = 1, colData = cdata)
se <- summarizeCircs(colData = cdata, wobble = 1, keepCols = 1:9)
se <- annotateCircs(se, annot.list, "mm9", fixCoordIndexing = TRUE)
tab <- resTable(se)
tab
