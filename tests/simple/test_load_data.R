
gtf2sqlite(assembly = "hg19", db.file="data/test.sqlite")

annot.list <- loadAnnotation("data/test.sqlite")

circs.f <- annotateCircs(circs.bed = "data/Sy5y_D0_sites.bed", annot.list = annot.list, assembly = "hg19")



circHist(circs.f, 0.5)
annotPie(circs.f, 0.02)

library(data.table)
library(GenomicRanges)
se <- summarizeCircs(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))], wobble=1)

circs <- lapply(dir("data/demo/", full.names=T)[grep("sites", dir("data/demo/"))], readCircs)
circs <- lapply(circs, circLinRatio, return.readcounts=TRUE)


circ.files <- dir("/clusterhome/pglazar/projects/circus/ciRcus/data/demo/", full.names=T)[grep("sites", dir("/clusterhome/pglazar/projects/circus/ciRcus/data/demo/"))]
keep.linear = TRUE
wobble = 1
subs = "all"
qualfilter = TRUE
keepCols = 1:6
colData = NULL


  circs <- lapply(circ.files, readCircs)
#linByCirc <- function(circs) {
  #circs = lapply(circ.files, readCircs, subs, qualfilter, keepCols)
  dcircs = rbindlist(circs)
  dcircs$type = ifelse(grepl('circ',dcircs$name),'circ','linear')
  dcircs = unique(dcircs[,.(chrom, start, end, strand, type)])
  dcircs = split(dcircs, dcircs$type)

  circ.gr =  makeGRangesFromDataFrame(as.data.frame(dcircs[['circ']]), keep.extra.columns=TRUE)
  lin.gr =  makeGRangesFromDataFrame(as.data.frame(dcircs[['linear']]), keep.extra.columns=TRUE)

  circ.gr.s = resize(resize(circ.gr, fix='start', width=1), fix='center', width=wobble)
  circ.gr.e = resize(resize(circ.gr, fix='end',   width=1), fix='center', width=wobble)

  cfos = data.table(as.matrix(findOverlaps(resize(lin.gr, fix='start', width=1), circ.gr.e, ignore.strand=FALSE)))
  cfoe = data.table(as.matrix(findOverlaps(resize(lin.gr, fix='end',   width=1), circ.gr.s, ignore.strand=FALSE)))

  cfos <- merge(cfos, data.table(rno=1:length(circ.gr), circID=paste0(seqnames(circ.gr), ":", start(circ.gr), "-", end(circ.gr))), by.x="subjectHits", by.y="rno")
  cfoe <- merge(cfoe, data.table(rno=1:length(circ.gr), circID=paste0(seqnames(circ.gr), ":", start(circ.gr), "-", end(circ.gr))), by.x="subjectHits", by.y="rno")

  s <- split(lin.gr[cfos$queryHits], factor(cfos$subjectHits, levels=1:length(circ.gr)))
  e <- split(lin.gr[cfoe$queryHits], factor(cfoe$subjectHits, levels=1:length(circ.gr)))

  names(s) <- cfos$circID[ match(names(s), cfos$subjectHits)]
  names(e) <- cfoe$circID[ match(names(e), cfoe$subjectHits)]

  se <- mergeGRL(s, e)

  circ.gr$circID <- paste0(seqnames(circ.gr), ":", start(circ.gr), "-", end(circ.gr))
  circ.gr$lin <- se



  tmp <- GRanges(seqnames=rep("chr1", 2),
                 ranges = IRanges(start=c(1000, 1000),
                                  end = c(2000, 2000)),
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
