
gtf2sqlite(assembly = "hg19", db.file="data/test.sqlite")

annot.list <- loadAnnotation("data/test.sqlite")

circs.f <- annotateCircs(circs.bed = "data/Sy5y_D0_sites.bed", annot.list = annot.list, assembly = "hg19")



circHist(circs.f, 0.5)
annotPie(circs.f, 0.02)

library(data.table)
library(GenomicRanges)
cdata <- data.frame(sample=c("FC1", "FC2", "H1", "H2", "L1", "L2"),
                    filename=basename(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))]))
se <- summarizeCircs(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))], wobble=1, colData = cdata)
se <- annotateHostGenes(se, annot.list$genes)
resTable(se)
se <- circLinRatio(se)
histogram(se)

# SH-SY5Y test
cdata <- data.frame(sample=c("D0"),
                    filename="data/Sy5y_D0_sites.bed")
se <- summarizeCircs(as.character(cdata$filename), wobble=1, colData = cdata)
se <- annotateHostGenes(se, annot.list$genes)
resTable(se)
se <- circLinRatio(se)
histogram(se)


circs <- lapply(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))], readCircs)
circs <- lapply(circs, circLinRatio, return.readcounts=TRUE)


circ.files <- dir("/clusterhome/pglazar/projects/circus/ciRcus/data/demo/", full.names=T)[grep("sites", dir("/clusterhome/pglazar/projects/circus/ciRcus/data/demo/"))]
keep.linear = TRUE
wobble = 1
subs = "all"
qualfilter = TRUE
keepCols = 1:6
colData = NULL




  tmp <- GRanges(seqnames=rep("chr1", 2),
                 ranges = IRanges(start=c(1000000, 1000000),
                                  end = c(1000100, 1000100)),
                 strand = c("+", "-"))
  #tmp <- shift(tmp, -99)
  #tmp <- resize(tmp, fix="start", width = 100)
  tmp <- c(flank(tmp, 100, start=TRUE), flank(tmp, 100, start=FALSE))


  getFlanks <- function(gr, length=100, type=c("circ", "lin")) {
    if (is.null(names(gr))) {
      names(gr) <- 1:length(gr)
    }

    if (type == "circ") {
      gr2 <- c(resize(gr, width = length, fix = "start"), resize(gr, width = length, fix = "end"))
    }
  }


  mergeGRL <- function(x, y) {

#    shared_keys <- intersect(names(x), names(y))
#    merged <- c(GRangesList(mapply(c, x[shared_keys], y[shared_keys])), x[!(names(x) %in% shared_keys)], y[!(names(y) %in% shared_keys)])
    xg = unlist(x)
    xg$.name = rep(names(x), times=elementLengths(x))
    yg = unlist(y)
    yg$.name = rep(names(y), times=elementLengths(y))
    gg = c(xg, yg)
    gg = split(gg, names(gg))

    return(gg)
  }
