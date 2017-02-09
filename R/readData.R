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
#'
# readCircs <- function(file, subs="all", qualfilter=TRUE, keepCols=1:6, ...) {
#
#   suppressWarnings(
#     DT <- fread(file, sep="\t", header = T) # maybe add colClasses later
#   )
#   setnames(DT, "# chrom", "chrom")
#   DT <- DT[-grep("#", DT$chrom)]
#
#   # change column classes where needed
#   # ***due to find_circ.py logic of putting a header line
#   #    in the middle of the output file, all columns are
#   #    character after fread()
#   for (col in names(DT)[c(2,3,5,7,8,9,10,11,13,14,15,16)]){
#     set(DT, j=col, value=as.integer(DT[[col]]))
#   }
#
#   if (subs != "all") {
#     DT <- DT[grep(subs, DT$name)]
#   }
#
#   if (qualfilter == TRUE) {
#     DT <- qualFilter(DT, ...)
#   }
#
#   DT <- DT[, keepCols, with=F]
#   DT$id <- paste(DT$chrom, ":", DT$start, "-", DT$end, sep="")
#
#   return(DT[, !"id", with=F])
# }
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
    DT <- fread(file, sep="\t", header=T) # maybe add colClasses later
  )

  # try to figure out what tool the data are coming from
  if (ncol(DT) == 19 & names(DT)[1] == '# chrom') {
    # read find_circ
    setnames(DT, "# chrom", "chrom")
    DT <- DT[!grepl("#", DT$chrom)]
    # change column classes where needed
    # ***due to find_circ.py logic of putting a header line
    #    in the middle of the output file, all columns are
    #    character after fread()
    char.class = c('chrom', 'name', 'strand', 'tissues', 'signal', 'strandmatch', 'category')
    for (col in setdiff(colnames(DT), char.class)){
      set(DT, j=col, value=as.integer(DT[[col]]))
    }

  } else if (ncol(DT) >= 21 & names(DT)[1] == '#chrom'){

    DT.lin <- fread(sub("circ", "lin", file), sep="\t", header=T)
    DT.lin$name <- sub("lin", "norm", DT.lin$name)
    DT <- rbind(DT.lin, DT)

    # read find_circ2
    setnames(DT, "#chrom", "chrom")
    setnames(DT, names(DT)[5], "n_reads") # TODO: this is dirty af

    DT <- DT[!grepl("#", DT$chrom)]
  }


  if (subs != "all") {
    DT <- DT[grep(subs, DT$name)]
  }

  if (qualfilter == TRUE) {
    DT <- qualFilter(DT)
  }

  DT <- DT[, keepCols, with=F]

  return(DT)
}

# ---------------------------------------------------------------------------- #
#' Function that reads multiple find_circ output files and returns a
#' SummarizedExperiment object. The resulting object contains unified circ RNA
#' coordinates as GRanges and counts of back-spliced and linearly spliced reads
#' for each circ RNA in each sample
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
#'
#'
#' @examples
#' circ.files = list.files(system.files('extdata'), full.names=TRUE, pattern=bed)
#' circs = summarizeCircs(circ.files)
#'
#' @docType methods
#' @rdname summarizeCircs-methods
#' @export
setGeneric("summarizeCircs",
           function(colData=NULL,
                    keep.linear=TRUE,
                    wobble = 5,
                    subs = 'all',
                    qualfilter=TRUE,
                    keepCols=1:6,
                    ...)
             standardGeneric("summarizeCircs"))


#' @rdname summarizeCircs-methods
#' @usage  \\S4method{summarizeCircs}{character}(files, keep.linear, wobble, subs, qualfilter, keepCols,colData)
setMethod("summarizeCircs",signature("data.frame"),
          function(colData, keep.linear, wobble, subs, qualfilter,keepCols){


            # -------------------------------------- #
            coldata.cnams = c('sample','filename')
            if(!all(coldata.cnams %in% colnames(colData)))
               stop(paste(setdiff(coldata.cnams, colnames(colData),
                            'is missing from colData')))

            # -------------------------------------- #
            circ.files = as.character(colData$filename)

            # -------------------------------------- #
            if(!all(file.exists(circ.files)))
              stop('Supplied circ files do not exist')

            circs = lapply(circ.files, readCircs, subs, qualfilter, keepCols)
            names(circs) <- colData$sample

            dcircs = rbindlist(circs)
            if (grepl("_circ_norm_", dcircs$name[1]) | grepl("_circ_circ_", dcircs$name[1])) {
              message('funky naming scheme used, will convert _circ_norm_ to _norm_ and _circ_circ_ to _circ_ before everything crashes')
              dcircs$name <- sub("_circ_norm_", "_norm_", dcircs$name) # TODO: add exception somewhere, if there is _circ_norm_ in names, boom
              dcircs$name <- sub("_circ_circ_", "_circ_", dcircs$name)
            }
            dcircs$type = ifelse(grepl('circ', dcircs$name), 'circ', 'linear')
            #dcircs$set = factor(sub('norm_.+','',sub('circ_.+','',dcircs$name))) # old
            dcircs$set = factor(rep(names(circs), sapply(circs, nrow)))
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

            #circ.gr.reduced = sort(unlist(range(split(circ.gr, merge.fos$fac), ignore.strand=FALSE)))
            circ.gr.reduced = unlist(range(split(circ.gr, merge.fos$fac), ignore.strand=FALSE))
            # -------------------------------------- #
            message('Fetching circular expression')
            circ.ex = merge.fos[,c(1,4), with=FALSE]
            circ.ex$nreads = circ.gr$n_reads[circ.ex$queryHits]
            circ.ex$set = circ.gr$set[circ.ex$queryHits]
            circ.ex.matrix = dcast.data.table(formula=fac~set,
                                              fun.aggregate=sum,
                                              fill=0,
                                              value.var='nreads',
                                              data=circ.ex)
            circ.ex.matrix = circ.ex.matrix[match(names(circ.gr.reduced),circ.ex.matrix$fac)]
            assays = list()
            assays$circ = as.matrix(circ.ex.matrix[,-1,with=FALSE])

            if(keep.linear==TRUE){
              message('Processing linear transcripts')
              linear = ProcessLinear(dcircs, circ.gr.reduced, wobble)
              assays = c(assays, linear)
            }

            if(class(colData) == 'data.frame')
                colData = DataFrame(colData)


            sex = SummarizedExperiment(assays=assays,
                                       rowRanges=circ.gr.reduced,
                                       colData=colData)
            return(sex)


})

#' @rdname summarizeCircs-methods
#' @usage  \\S4method{summarizeCircs}{character}(files, keep.linear, wobble, subs, qualfilter, keepCols,colData)
setMethod("summarizeCircs",signature("character"),
          function(colData, keep.linear, wobble, subs, qualfilter,keepCols){



            message('Constructing colData...')
            colData = data.frame(sample = sub('.candidates.bed','',basename(colData)),
                                filename=colData,
                                stringsAsFactors=FALSE)

            summarizeCircs(colData=colData,
                           keep.linear=keep.linear,
                           wobble=wobble,
                           subs=subs,
                           qualfilter=qualfilter,
                           keepCols = keepCols)
})


# ---------------------------------------------------------------------------- #
# Function that extracts the linear splicing isoforms for each circ RNA
ProcessLinear = function(dcircs, circ.gr.reduced, wobble){

    lin.gr =  makeGRangesFromDataFrame(as.data.frame(dcircs[['linear']]),
                                       keep.extra.columns=TRUE)
    circ.gr.s = resize(resize(circ.gr.reduced, fix='start', width=1),
                       fix='center', width=wobble)
    circ.gr.e = resize(resize(circ.gr.reduced, fix='end',   width=1),
                       fix='center', width=wobble)

    cfos = data.table(as.matrix(findOverlaps(resize(lin.gr, fix='start', width=1),
                                             circ.gr.e, ignore.strand=FALSE)))
    cfoe = data.table(as.matrix(findOverlaps(resize(lin.gr, fix='end',   width=1),
                                             circ.gr.s, ignore.strand=FALSE)))

    cfos$nreads = lin.gr$n_reads[cfos$queryHits]
    cfos$set = lin.gr$set[cfos$queryHits]
    cfos$queryHits = factor(cfos$subjectHits, levels=1:length(circ.gr.reduced))
    cfos.cast = dcast.data.table(formula=queryHits~set, fun.aggregate=sum,fill=0,
                                 value.var='nreads', data=cfos, drop=FALSE)
    cfos.cast = cfos.cast[match(names(circ.gr.reduced),cfos.cast$queryHits)]

    cfoe$nreads = lin.gr$n_reads[cfoe$queryHits]
    cfoe$set = lin.gr$set[cfoe$queryHits]
    cfoe$queryHits = factor(cfoe$subjectHits, levels=1:length(circ.gr.reduced))
    cfoe.cast = dcast.data.table(formula=queryHits~set, fun.aggregate=sum,fill=0,
                                 value.var='nreads', data=cfoe, drop=FALSE)
    cfoe.cast = cfoe.cast[match(names(circ.gr.reduced),cfoe.cast$queryHits)]

    return(list(linear.start = as.matrix(cfos.cast[,-1,with=FALSE]),
                linear.end   = as.matrix(cfoe.cast[,-1,with=FALSE])))
}

