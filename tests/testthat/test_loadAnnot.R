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
test_that("demo data are read properly", {

  se <- summarizeCircs(dir("../../data/demo/", full.names=T)[grep("sites", dir("../../data/demo/"))], wobble=1)

  # score matrix should be 4 by 6 (4 circs, 6 samples)
  expect_equal(dim(assays(se)$circ),     c(4,6))
  # the total number of circRNA reads in FrontalCortex rep1 should be 1574
  expect_equal(sum(assays(se)$circ[,1]),  1574)
  # circ read counts
  expect_equal(assays(se)$circ[[2,3]],     14892)
  expect_equal(apply(assays(se)$circ, 2, sum), c(1574, 757, 17266, 7221, 67, 93))

  }
)
