test_that("region parser accepts standard VCF intervals", {
  x <- goatHapNet:::parse_region("chr5:12300000-12400000")
  expect_equal(x$chr, "chr5")
  expect_equal(x$start, 12300000L)
  expect_equal(x$end, 12400000L)
})

test_that("prefix sanitizer is file-system friendly", {
  expect_equal(goatHapNet:::sanitize_prefix("chr5:1-10"), "chr5_1-10")
})
