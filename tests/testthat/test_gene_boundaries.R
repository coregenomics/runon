library(runon)
context("Gene boundary discovery")

gr <- GenomicRanges::GRanges(
    seqnames = c("chr1:2500-2800",
                 "chr2:2900-3300",
                 "chrUn_gibberish:3500-3800",
                 "chr1_gibberish_random:3900-4300",
                 "chr2_gibberish_alt:4500-4800",
                 "chrX:4900-5300"))

test_that("limitChromosomes drops unmappable chromosome names", {
    expect_equal(
