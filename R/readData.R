# ---------------------------------------------------------------------------- #
#' read a tabular circRNA candidate list
#'
#' description
#'
#' details
#' @param file location of the input file, a character string such as:
#'             "/home/user/find_circ_sites.bed"
#' @param subs a character string, keep only lines containing it in the name column
#' @param qualFilter a boolean, tells whether the quality filtering should be performed
#' @param keepCols a vector of column numbers return
#' @export
readCircs <- function(file, subs="all", qualfilter=TRUE, keepCols=1:6, ...) {

  suppressWarnings(
    DT <- fread(file, sep="\t", header = T) # maybe add colClasses later
  )
  setnames(DT, "# chrom", "chrom")
  DT <- DT[!grepl("#", DT$chrom)]

  # change column classes where needed
  # ***due to find_circ.py logic of putting a header line
  #    in the middle of the output file, all columns are
  #    character after fread()
  char.class = c('chrom','name','strand','tissues','signal','strandmatch','category')
  for (col in setdiff(colnames(DT),char.class)){
    set(DT, j=col, value=as.integer(DT[[col]]))
  }

  if (subs != "all") {
    DT <- DT[grep(subs, DT$name)]
  }

  if (qualfilter == TRUE) {
    DT <- qualFilter(DT)
  }

  DT <- DT[, keepCols, with=F]
  DT$id <- paste(DT$chrom, ":", DT$start, "-", DT$end, sep="")

  return(DT[, !"id", with=F])
}

# ---------------------------------------------------------------------------- #
#' Function that reads multiple find_circ output files and returns a
#' SummarizedExperiment object
#'
#'
#' @param circ.files a character vector of paths to find_circ output files
#' @param keep.linear a boolean indicating whether to keep the counts of
#'        linearly spliced transcripts
#'
#' @param wobble number of nucleotides around the splicing border that should be
#'        considered when collapsing circular transcripts - helps with mapping
#'        imprecision
#' @param subs a character string, keep only lines containing it in the name column
#' @param qualfilter a boolean, tells whether the quality filtering should be performed
#' @param keepCols a vector of column numbers return
#' @param colData a \code{DataFrame} object that contains the experiment data. It has
#'        to have the same number of rows as the number of files
#'
#' @return returns a \code{SummarizedExperiment}
#' @examples
#'
#'
#' @docType methods
#' @rdname summarizeCircs-methods
#' @export
setGeneric("summarizeCircs",
           function(circ.files,
                    keep.linear=TRUE,
                    wobble = 5,
                    subs = 'all',
                    qualfilter=TRUE,
                    keepCols=1:6,
                    colData=NULL,
                    ...)
             standardGeneric("summarizeCircs"))


#' @rdname summarizeCircs-methods
#' @usage  \\S4method{summarizeCircs}{character}(files, keep.linear, wobble, subs, qualfilter, keepCols,colData)
setMethod("summarizeCircs",signature("character"),
          function(circ.files, keep.linear, wobble, subs, qualfilter,keepCols, colData){

            # -------------------------------------- #
            circs = lapply(circ.files, readCircs, subs, qualfilter, keepCols)
            dcircs = rbindlist(circs)
            dcircs$type = ifelse(grepl('circ',dcircs$name),'circ','linear')
            dcircs$set = factor(sub('norm_.+','',sub('circ_.+','',dcircs$name)))
            dcircs = split(dcircs, dcircs$type)

            # -------------------------------------- #
            message('Processing circular transcripts')
            circ.gr =  makeGRangesFromDataFrame(as.data.frame(dcircs[['circ']]), keep.extra.columns=TRUE)

            circ.gr.s = resize(resize(circ.gr, fix='start', width=1), fix='center', width=wobble)
            circ.gr.e = resize(resize(circ.gr, fix='end', width=1), fix='center',   width=wobble)

            circ.fos = data.table(as.matrix(findOverlaps(circ.gr.s, reduce(circ.gr.s, ignore.strand=FALSE))))
            circ.foe = data.table(as.matrix(findOverlaps(circ.gr.e, reduce(circ.gr.e, ignore.strand=FALSE))))

            merge.fos = merge(circ.fos, circ.foe, by='queryHits', all=TRUE)
            merge.fos$fac = with(merge.fos, as.numeric(factor(paste(subjectHits.x, subjectHits.y))))

            circ.gr.reduced = sort(unlist(range(split(circ.gr, merge.fos$fac), ignore.strand=FALSE)))

            # -------------------------------------- #
            message('Fetching circular expression')
            circ.ex = merge.fos[,c(1,4), with=FALSE]
            circ.ex$nreads = circ.gr$n_reads[circ.ex$queryHits]
            circ.ex$set = circ.gr$set[circ.ex$queryHits]
            circ.ex.matrix = dcast.data.table(formula=fac~set, fun.aggregate=sum, fill=0, value.var='nreads', data=circ.ex)
            circ.ex.matrix = circ.ex.matrix[match(names(circ.gr.reduced),circ.ex.matrix$fac)]
            assays = list()
            assays$circ = as.matrix(circ.ex.matrix[,-1,with=FALSE])

            if(keep.linear==TRUE){
                message('Processing linear transcripts')
                linear = ProcessLinear(dcircs, circ.gr.reduced, wobble)
                assays = c(assays, linear)
            }


            if(is.null(colData))
              colData = DataFrame(sample = sub('.candidates.bed','',basename(circ.files)))

            if(class(colData) == 'data.frame')
              colData = DataFrame(colData)

            sex = SummarizedExperiment(assays=assays,
                                       rowRanges=circ.gr.reduced,
                                       colData=colData)
            return(sex)

})

ProcessLinear = function(dcircs, circ.gr.reduced, wobble){

    lin.gr =  makeGRangesFromDataFrame(as.data.frame(dcircs[['linear']]), keep.extra.columns=TRUE)
    circ.gr.s = resize(resize(circ.gr.reduced, fix='start', width=1), fix='center', width=wobble)
    circ.gr.e = resize(resize(circ.gr.reduced, fix='end',   width=1), fix='center', width=wobble)

    cfos = data.table(as.matrix(findOverlaps(resize(lin.gr, fix='start', width=1), circ.gr.e, ignore.strand=FALSE)))
    cfoe = data.table(as.matrix(findOverlaps(resize(lin.gr, fix='end',   width=1), circ.gr.s, ignore.strand=FALSE)))

    cfos$nreads = lin.gr$n_reads[cfos$queryHits]
    cfos$set = lin.gr$set[cfos$queryHits]
    cfos$queryHits = factor(cfos$subjectHits, levels=1:length(circ.gr.reduced))
    cfos.cast = dcast.data.table(formula=queryHits~set, fun.aggregate=sum,fill=0, value.var='nreads', data=cfos, drop=FALSE)
    cfos.cast = cfos.cast[match(names(circ.gr.reduced),cfos.cast$queryHits)]

    cfoe$nreads = lin.gr$n_reads[cfoe$queryHits]
    cfoe$set = lin.gr$set[cfoe$queryHits]
    cfoe$queryHits = factor(cfoe$subjectHits, levels=1:length(circ.gr.reduced))
    cfoe.cast = dcast.data.table(formula=queryHits~set, fun.aggregate=sum,fill=0, value.var='nreads', data=cfoe, drop=FALSE)
    cfoe.cast = cfoe.cast[match(names(circ.gr.reduced),cfoe.cast$queryHits)]

    return(list(linear.start = as.matrix(cfos.cast[,-1,with=FALSE]),
                linear.end   = as.matrix(cfoe.cast[,-1,with=FALSE])))
}

