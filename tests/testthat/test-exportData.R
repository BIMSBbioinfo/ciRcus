test_that("bedTracks", {

  circ.files <- list.files(system.file("extdata/encode_demo_small",
                                       package = "ciRcus"),
                           pattern = "sites.bed",
                           full.names = TRUE)
  se <- summarizeCircs(circ.files, wobble = 1, keepCols = 1:7)

  bed.tracks <- bedTracks(se)

  # there should be one BED track per sample
  expect_equal(length(bed.tracks), ncol(se))

  # BED entry names should be identical to circRNA names
  expect_equal(sort(unique(names(unlist(unname(bed.tracks))))),
               sort(rownames(se)))

  # samples should only have BED entries for detected samples
  expect_equal(unname(elementNROWS(bed.tracks)), c(3, 3, 4, 4, 3, 3))

  # scores should be truncated at maximal score
  expect_equal(max(rtracklayer::score(unlist(bedTracks(se, max.score = 10)))),
               10)

  # score should be `.` if no assay was specified
  expect_equal(as.vector(as.matrix(unique(mcols(unlist(bedTracks(se, NULL)))))),
               ".")

  # not specifying a score assay should result in a message
  expect_message(bedTracks(se, NULL),
                 "no BED score defined; using `.` as placeholder",
                 fixed = TRUE)

  # specifying > 1 score assay should result in a warning
  expect_warning(bedTracks(se, c("circ", "circ.uniq")),
                 "length(score) > 1; only first entry will be used",
                 fixed = TRUE)

  # specifying > 1 minimal score should result in a warning
  expect_warning(bedTracks(se, min.score = c(50, 1)),
                 "length(min.score) > 1; only first entry will be used",
                 fixed = TRUE)

  # specifying > 1 maximal score should result in a warning
  expect_warning(bedTracks(se, max.score = c(100, 1000)),
                 "length(max.score) > 1; only first entry will be used",
                 fixed = TRUE)

  # specifying a non-existent score assay should result in a warning
  expect_warning(bedTracks(se, "dummy_assay"),
                 paste0("no assay named 'dummy_assay'",
                       "; BED output will be generated without scores"),
                 fixed = TRUE)

  # not specifying a score assay but a minimal score should result in warning
  expect_warning(bedTracks(se, NULL, 1),
                 "no BED score defined; circRNAs will not be filtered",
                 fixed = TRUE)

  # not specifying a score assay but a maximal score should result in warning
  expect_warning(bedTracks(se, NULL, max.score = 1000),
                 "no BED score defined; BED scores will not be truncated",
                 fixed = TRUE)
})
