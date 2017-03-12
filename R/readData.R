# --------------------------------------------------------------------------- #
#' read a tabular circRNA candidate list
#'
#' \code{readCircs} loads a list of circRNA candidates into a \code{data.table}.
#' Currently supported input formats are find_circ.py and find_circ2.py.
#'
#' Not intended to be used directly, but will stay exported for the time being.
#' If find_circ2.py is used, fifth column will be renamed to \code{n_reads}
#' (from \code{n_frags} or \code{counts}, depending on the find_circ2 release) for
#' backwards compatibility and old time's sake. Will rename "lin" to "norm" in
#' find_circ2 input name column for the same reason.
#' Expects \code{<myproject>/fc2/circ_splice_sites.bed} as input, and assumes
#' the existence of \code{<myproject>/fc2/lin_splice_sites.bed}.
#'
#' @param file Input file location, a character string such as
#'             \code{/home/user/my_circRNA_project/circ_splice_sites.bed}
#' @param subs A character string, keep only lines containing it in the name column.
#' @param qualfilter should quality filtering be performed?
#' @param keepCols An integer vector. Which input columns should be returned?
#' @param ... other arguments
#'
#' @return A data table.
#'
#' @export
readCircs <- function(file, subs = "all", qualfilter = TRUE, keepCols = 1:6, ...) {

  suppressWarnings(
    DT <- fread(file, sep = "\t", header = T) # maybe add colClasses later
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
      set(DT, j = col, value = as.integer(DT[[col]]))
    }

  } else if (ncol(DT) >= 21 & names(DT)[1] == '#chrom'){
    # read find_circ2
    DT.lin <- fread(sub("circ_splice_sites.bed$", "lin_splice_sites.bed", file), sep = "\t", header = T)
    DT.lin$name <- sub("lin", "norm", DT.lin$name)
    DT <- rbind(DT.lin, DT)

    # read find_circ2
    setnames(DT, "#chrom", "chrom")
    setnames(DT, names(DT)[5], "n_reads") # TODO: this is dirty af

    DT <- DT[!grepl("#", DT$chrom)]
  } else if (ncol(DT) == 12 & names(DT)[1] == 'circRNA_ID') {
    # read CIRI2
    setnames(DT, c('name', 'chrom', 'start', 'end', 'n_reads', 'SM_MS_SMS', 'n_reads_nonjunction', 'junction_reads_ratio', 'circRNA_type', 'gene_id', 'strand', 'junction_reads_ID'))
  }


  if (subs != "all") {
    DT <- DT[grep(subs, DT$name)]
  }

  if (qualfilter == TRUE) {
    DT <- qualFilter(DT)
  }

  DT <- DT[, keepCols, with = F]

  return(DT)
}

# --------------------------------------------------------------------------- #
#' load circRNA detection output to a SummarizedExperiment object
#'
#' Function that reads an arbitrary number of circRNA lists and returns a
#' SummarizedExperiment object. The resulting object contains unified circ RNA
#' coordinates as GRanges and counts of back-spliced and linearly spliced reads
#' for each circRNA in each sample
#'
#'
#' @param keep.linear A boolean indicating whether to keep the counts of
#'        linearly spliced transcripts
#' @param wobble Number of nucleotides around the splicing border that should be
#'        considered when collapsing circular transcripts - helps with mapping
#'        imprecision
#' @param subs A character string, keep only lines containing it in the name column
#' @param qualfilter A boolean. Should quality filtering be performed?
#' @param keepCols An integer vector. Which input columns should be returned?
#' @param colData A \code{DataFrame} object that contains the input files
#'        and sample names to be used for further analysis
#' @param ... other arguments
#'
#' @return A \code{SummarizedExperiment} object.
#'
#'
#' @docType methods
#' @rdname summarizeCircs-methods
#'
#' @export
setGeneric("summarizeCircs",
           function(colData = NULL,
                    keep.linear = TRUE,
                    wobble = 5,
                    subs = 'all',
                    qualfilter = TRUE,
                    keepCols = 1:6,
                    ...)
             standardGeneric("summarizeCircs"))

#' @aliases summarizeCircs,data.frame-method
#' @rdname summarizeCircs-methods
setMethod("summarizeCircs", signature("data.frame"),
          function(colData, keep.linear, wobble, subs, qualfilter, keepCols){


            # munge input
            # -------------------------------------- #
            coldata.cnams = c('sample','filename')
            if(!all(coldata.cnams %in% colnames(colData)))
               stop(paste(setdiff(coldata.cnams, colnames(colData),
                            'is missing from colData')))

            circ.files = as.character(colData$filename)

            if(!all(file.exists(circ.files)))
              stop('Supplied circ files do not exist')

            # load circs
            # -------------------------------------- #
            circs = lapply(circ.files, readCircs, subs, qualfilter, keepCols)
            names(circs) <- colData$sample

            dcircs = rbindlist(circs)
            dcircs$set = factor(rep(names(circs), sapply(circs, nrow)), levels = names(circs))

            # process circular and linear if input is find_circ
            if (!("SM_MS_SMS" %in% names(circs[[1]]))) { #TODO: find a better way to recognize CIRI2
              # find_circ
              if (grepl("_circ_norm_", dcircs$name[1]) | grepl("_circ_circ_", dcircs$name[1])) {
                message('funky naming scheme used, will convert _circ_norm_ to _norm_ and _circ_circ_ to _circ_ before everything crashes')
                dcircs$name <- sub("_circ_norm_", "_norm_", dcircs$name)
                dcircs$name <- sub("_circ_circ_", "_circ_", dcircs$name)
              }
              dcircs$type = ifelse(grepl('circ', dcircs$name), 'circ', 'linear')

              dcircs = split(dcircs, dcircs$type)

              circ.gr =  makeGRangesFromDataFrame(as.data.frame(dcircs[['circ']]), keep.extra.columns = TRUE)
            } else {
              # CIRI2
              circ.gr =  makeGRangesFromDataFrame(as.data.frame(dcircs), keep.extra.columns = TRUE)
            }

            # prepare the wobble
            # -------------------------------------- #
            circ.gr.s = resize(resize(circ.gr, fix = 'start', width = 1), fix = 'center', width = wobble)
            circ.gr.e = resize(resize(circ.gr, fix = 'end',   width = 1), fix = 'center', width = wobble)

            circ.fos = data.table(as.matrix(findOverlaps(circ.gr.s, reduce(circ.gr.s, ignore.strand = FALSE))))
            circ.foe = data.table(as.matrix(findOverlaps(circ.gr.e, reduce(circ.gr.e, ignore.strand = FALSE))))

            merge.fos = merge(circ.fos, circ.foe, by = 'queryHits', all = TRUE)
            merge.fos$fac = with(merge.fos, as.numeric(factor(paste(subjectHits.x, subjectHits.y))))

            #circ.gr.reduced = sort(unlist(range(split(circ.gr, merge.fos$fac), ignore.strand = FALSE)))
            # TODO: wobble logic is naive, should be improved to select the best expressed circRNA as a referent one
            circ.gr.reduced = unlist(range(split(circ.gr, merge.fos$fac), ignore.strand = FALSE))

            # prepare the assays object
            # TODO: this screams refactoring
            # -------------------------------------- #
            message('Fetching circular expression')

            assays = list()
            if (!("SM_MS_SMS" %in% names(circs[[1]]))) { #TODO: find a better way to recognize CIRI2

              n_reads.dt <- MungeColumn(merge.fos, circ.gr, circ.gr.reduced, "n_reads")
              n_uniq.dt  <- MungeColumn(merge.fos, circ.gr, circ.gr.reduced, "n_uniq")

              assays$circ      = as.matrix(n_reads.dt[, -1, with = FALSE])
              assays$circ.uniq = as.matrix( n_uniq.dt[, -1, with = FALSE])

            } else {

              n_reads.dt <- MungeColumn(merge.fos, circ.gr, circ.gr.reduced, "n_reads")
              assays$circ      = as.matrix(n_reads.dt[, -1, with = FALSE])

            }

            if(keep.linear == TRUE){
              message('Processing linear transcripts')
              linear = ProcessLinear(dcircs, circ.gr.reduced, wobble)
              assays = c(assays, linear)
            }



            if(class(colData) == 'data.frame')
                colData = DataFrame(colData)


            sex = SummarizedExperiment(assays = assays,
                                       rowRanges = circ.gr.reduced,
                                       colData = colData)
            return(sex)


})

#' @aliases summarizeCircs,character-method
#' @rdname summarizeCircs-methods
setMethod("summarizeCircs", signature("character"),
          function(colData, keep.linear, wobble, subs, qualfilter,keepCols){

            message('Constructing colData...')
            colData = data.frame(sample   = sub('.candidates.bed', '', basename(colData)),
                                 filename = colData,
                                 stringsAsFactors = FALSE)

            summarizeCircs(colData = colData,
                           keep.linear = keep.linear,
                           wobble = wobble,
                           subs = subs,
                           qualfilter = qualfilter,
                           keepCols = keepCols)
})


#' Title
#'
#' Function that, based on circRNA candidate list and collapsed circRNA candidate list
#' summarizes a numeric input column into a matrix that can be hooked to SummarizedExperiment
#'
#' @param merge.fos merge fos
#' @param circ.gr circs
#' @param circ.gr.reduced reduced circs
#' @param column.name column to extract
#'
#' @return a matrix
#' @export
MungeColumn <- function(merge.fos, circ.gr, circ.gr.reduced, column.name) {

  if (!(column.name %in% colnames(elementMetadata(circ.gr)))) {
    stop('unknown column name: ', column.name)
  }

  circ.ex = merge.fos[,.(queryHits, fac)]
  circ.ex$nreads = values(circ.gr)[[column.name]][circ.ex$queryHits]
  circ.ex$set = circ.gr$set[circ.ex$queryHits]
  circ.ex.matrix = dcast.data.table(formula = fac~set,
                                    fun.aggregate = sum,
                                    fill = 0,
                                    value.var = 'nreads',
                                    data = circ.ex)
  circ.ex.matrix = circ.ex.matrix[match(names(circ.gr.reduced), circ.ex.matrix$fac)]

  return(circ.ex.matrix)
}
# ---------------------------------------------------------------------------- #
#' Title
#'
#' Function that extracts the linear splicing isoforms for each circ RNA
#'
#' @param dcircs dcircs
#' @param circ.gr.reduced reduced circs
#' @param wobble how many nucleotides of wobble to tolerate?
#'
#' @return a list
#' @export
ProcessLinear = function(dcircs, circ.gr.reduced, wobble){

    lin.gr =  makeGRangesFromDataFrame(as.data.frame(dcircs[['linear']]),
                                       keep.extra.columns = TRUE)
    circ.gr.s = resize(resize(circ.gr.reduced, fix = 'start', width = 1),
                       fix = 'center', width = wobble)
    circ.gr.e = resize(resize(circ.gr.reduced, fix = 'end',   width = 1),
                       fix = 'center', width = wobble)

    cfos = data.table(as.matrix(findOverlaps(resize(lin.gr, fix = 'start', width = 1),
                                             circ.gr.e, ignore.strand = FALSE)))
    cfoe = data.table(as.matrix(findOverlaps(resize(lin.gr, fix = 'end',   width = 1),
                                             circ.gr.s, ignore.strand = FALSE)))

    cfos$nreads = lin.gr$n_reads[cfos$queryHits]
    cfos$set = lin.gr$set[cfos$queryHits]
    cfos$queryHits = factor(cfos$subjectHits, levels = 1:length(circ.gr.reduced))
    cfos.cast = dcast.data.table(formula = queryHits~set, fun.aggregate = sum,fill = 0,
                                 value.var = 'nreads', data = cfos, drop = FALSE)
    cfos.cast = cfos.cast[match(names(circ.gr.reduced),cfos.cast$queryHits)]

    cfoe$nreads = lin.gr$n_reads[cfoe$queryHits]
    cfoe$set = lin.gr$set[cfoe$queryHits]
    cfoe$queryHits = factor(cfoe$subjectHits, levels = 1:length(circ.gr.reduced))
    cfoe.cast = dcast.data.table(formula = queryHits~set, fun.aggregate = sum,fill = 0,
                                 value.var = 'nreads', data = cfoe, drop = FALSE)
    cfoe.cast = cfoe.cast[match(names(circ.gr.reduced),cfoe.cast$queryHits)]

    return(list(linear.start = as.matrix(cfos.cast[,-1,with = FALSE]),
                linear.end   = as.matrix(cfoe.cast[,-1,with = FALSE])))
}

