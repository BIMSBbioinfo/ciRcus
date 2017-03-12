test_that("summarizeCircs", {


  circ.files <- list.files(system.file('extdata/encode_demo_small', package = 'ciRcus'),
                          pattern = 'sites.bed',
                          full.names = TRUE)
  circ.files <- circ.files[!grepl('Sy5y', circ.files)]
  se <- summarizeCircs(circ.files, wobble = 1, keepCols = 1:7)

  # score matrix should be 4 by 6 (4 circs, 6 samples)
  expect_equal(dim(assays(se)$circ),     c(4,6))
  # there should be 4 assays
  expect_equal(length(assays(se)), 4)
  # the total number of circRNA reads in FrontalCortex rep1 should be 1574
  expect_equal(sum(assays(se)$circ[,1]),  1574)
  # circ read counts
  expect_equal(assays(se)$circ[[3,3]],     14892)
  expect_equal(unname(colSums(assays(se)$circ)), c(1574, 757, 17266, 7221, 67, 93))
  # lin read counts
  expect_equal(assays(se)$linear.start[[3,3]], 2362)
  # lin read count for a non-existing circRNA
  expect_equal(assays(se)$linear.end[[2,1]], 90)
  # circular unique read counts
  expect_equal(unname(rowSums(assays(se)$circ.uniq)), c(1105, 40, 4866, 218))
  # resTable tests
  expect_equal(dim(resTable(se)), c(4, 29))

})


test_that('Annotation',{

  circ.files = list.files(system.file('extdata/encode_demo_small', package = 'ciRcus'),
                          pattern = 'sites.bed',
                          full.names = TRUE)
  circ.files = circ.files[!grepl('Sy5y', circ.files)]
  se <- summarizeCircs(circ.files, wobble = 1, keepCols = 1:7)


  # annotation tests
  annot.file = system.file('extdata/db/hsa_ens75_minimal.sqlite', package = 'ciRcus')
  annot.list <- suppressMessages(loadAnnotation(annot.file))
  se <- annotateHostGenes(se, annot.list$genes)
  expect_equal(resTable(se)$gene_id, c("ENSG00000183023",
                                       "ENSG00000183023",
                                       "ENSG00000183023",
                                       "ENSG00000180357"))
  se <- annotateFlanks(se, annot.list$gene.feats)
  expect_equal(resTable(se)$feature[2], "utr5:cds")
  expect_equal(resTable(se)$feature[4], "utr5:cds")
  se <- annotateJunctions(se, annot.list$junctions)
  expect_equal(resTable(se)$junct.known, c("5pr", "5pr", "5pr", "both"))
  se <- circLinRatio(se)
  expect_equal(unname(assays(se)$ratio[3,4]), 4.03)
})

test_that('sample labels are robust upon resorting colData', {

  circ.files = list.files(system.file('extdata/encode_demo_small', package = 'ciRcus'),
                          pattern = 'sites.bed',
                          full.names = TRUE)
  circ.files = circ.files[!grepl('Sy5y', circ.files)]

  se.sorted   <- summarizeCircs(circ.files,      wobble = 1, keepCols = 1:7)
  se.unsorted <- summarizeCircs(rev(circ.files), wobble = 1, keepCols = 1:7)

  # rownames should be equal to sample column
  expect_equal(rownames(colData(se.sorted)),   colData(se.sorted)$sample)
  expect_equal(rownames(colData(se.unsorted)), colData(se.unsorted)$sample)

  # sorted vs. sorted
  expect_equal(unname(assays(se.sorted)$circ[,'FrontalCortex_rep1_sites.bed']),   resTable(se.sorted)$FrontalCortex_rep1_sites.bed_circ)
  # unsorted vs. unsorted
  expect_equal(unname(assays(se.unsorted)$circ[,'FrontalCortex_rep1_sites.bed']), resTable(se.unsorted)$FrontalCortex_rep1_sites.bed_circ)

  # sorted vs. unsorted, assays()
  expect_equal(assays(se.sorted)$circ[,'FrontalCortex_rep1_sites.bed'],   assays(se.unsorted)$circ[,'FrontalCortex_rep1_sites.bed'])
  # sorted vs. unsorted, resTable()
  expect_equal(resTable(se.sorted)$FrontalCortex_rep1_sites.bed_circ, resTable(se.unsorted)$FrontalCortex_rep1_sites.bed_circ)
})

test_that('CIRI2 input can be digested by ciRcus', {

  circ.files = list.files(system.file('extdata/ciri_demo_hek', package = 'ciRcus'),
                          pattern = 'HEK',
                          full.names = TRUE)

  se <- summarizeCircs(circ.files, keep.linear = FALSE, wobble = 1, subs = "all", qualfilter = FALSE, keepCols = 1:12)

  # total circRNAs in both samples
  expect_equal(nrow(resTable(se)), 14936)
  # there should be 4127 ribozero circRNAs
  expect_equal(sum(resTable(se)[, 6, with = F] > 0), 4127)
  # there should be 14013 RNaseR circRNAs
  expect_equal(sum(resTable(se)[, 7, with = F] > 0), 14013)
})
