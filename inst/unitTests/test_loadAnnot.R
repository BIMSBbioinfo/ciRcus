# test_that("annotation list loads properly", {
#
#   annot.list <- loadAnnotation("../../data/test.sqlite")
#
#   expect_equal(length(annot.list$genes),                 57773)
#   expect_equal(length(annot.list$gene.feats$utr5),       58834)
#   expect_equal(length(annot.list$gene.feats$cds),       209804)
#   expect_equal(length(annot.list$junctions$start),      562242)
#   expect_equal(mean(start(annot.list$junctions$end)), 74316829)
#
#   circs.f <- annotateCircs(circs.bed = "../../data/Sy5y_D0_sites.bed", annot.list = annot.list, assembly = "hg19")
#
#   expect_equal(nrow(circs.f),             2166)
#   expect_equal(sum(circs.f$n_reads),      9182)
#   expect_equal(median(circs.f$start), 69383524)
#
# })
test_that("summarizeCircs", {


  circ.files = list.files(system.file('extdata', package='ciRcus'),
                          pattern='sites.bed',
                          full.names=TRUE)
  circ.files = circ.files[!grepl('Sy5y', circ.files)]
  se <- summarizeCircs(circ.files, wobble=1)

  # score matrix should be 4 by 6 (4 circs, 6 samples)
  expect_equal(dim(assays(se)$circ),     c(4,6))
  # the total number of circRNA reads in FrontalCortex rep1 should be 1574
  expect_equal(sum(assays(se)$circ[,1]),  1574)
  # circ read counts
  expect_equal(assays(se)$circ[[3,3]],     14892)
  expect_equal(apply(assays(se)$circ, 2, sum), c(1574, 757, 17266, 7221, 67, 93))
  # lin read counts
  expect_equal(assays(se)$linear.start[[3,3]], 2362)
  # lin read count for a non-existing circRNA
  expect_equal(assays(se)$linear.end[[2,1]], 90)

  # resTable tests
  #expect_equal(resTable(se[, se$sample=="FC1"])$FC1_lin.start[3], 449L)
  expect_equal(dim(resTable(se)), c(4, 23))

})


test_that('Annotation',{

  circ.files = list.files(system.file('extdata', package='ciRcus'),
                          pattern='sites.bed',
                          full.names=TRUE)
  circ.files = circ.files[!grepl('Sy5y', circ.files)]
  se <- summarizeCircs(circ.files, wobble=1)


  # annotation tests
  annot.file = system.file('extdata/hsa_ens75_minimal.sqlite', package='ciRcus')
  annot.list <- suppressMessages(loadAnnotation(annot.file))
  se <- annotateHostGenes(se, annot.list$genes)
  expect_equal(resTable(se)$gene_id, c("ENSG00000183023",
                                       "ENSG00000183023",
                                       "ENSG00000183023",
                                       "ENSG00000180357"))
  se <- annotateFlanks(se, annot.list$gene.feats)
  expect_equal(resTable(se)$feature[2], "cds:utr5")
  expect_equal(resTable(se)$feature[4], "utr5:cds")
  se <- annotateJunctions(se, annot.list$junctions)
  expect_equal(resTable(se)$junct.known, c("5pr", "5pr", "5pr", "both"))
  se <- circLinRatio(se)
  expect_equal(unname(assays(se)$ratio[3,4]), 4.03)
})
