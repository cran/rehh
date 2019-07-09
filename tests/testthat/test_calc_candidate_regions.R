#checks calculation of candidate regions of selection
context("calc_candidate_regions")

setup({
  library(rehh.data)
  data(wgscan.cgu)
  data(wgscan.eut)
})

test_that("checked_calc_candidate_regions", {
  sink("test_calc_candidate_regions.log")
  ihs.cgu <- ihh2ihs(wgscan.cgu)
  cr.ihs.cgu <-
    calc_candidate_regions(
      ihs.cgu,
      threshold = 4,
      pval = TRUE,
      overlap = 100000,
      min_n_extr_mrk = 2,
      min_perc_extr_mrk = 20
    )
  chr_numbers <- as.integer(levels(cr.ihs.cgu$CHR)[cr.ihs.cgu$CHR])
  expect_identical(chr_numbers, c(5L, 12L, 18L))
  
  cr.ihs.cgu <-
    calc_candidate_regions(
      ihs.cgu,
      threshold = 4,
      pval = TRUE,
      overlap = 100000,
      min_n_extr_mrk = 2
    )
  chr_numbers <- as.integer(levels(cr.ihs.cgu$CHR)[cr.ihs.cgu$CHR])
  expect_identical(chr_numbers, c(1L, 4L, 5L, 5L, 7L, 12L, 18L))
  
  rsb.cgu_eut <- ines2rsb(wgscan.cgu, wgscan.eut)
  cr.rsb.cgu_eut <-
    calc_candidate_regions(
      rsb.cgu_eut,
      threshold = 4,
      pval = TRUE,
      overlap = 100000,
      min_n_extr_mrk = 2
    )
  chr_numbers <-
    
    as.integer(levels(cr.rsb.cgu_eut$CHR)[cr.rsb.cgu_eut$CHR])
  expect_identical(chr_numbers,
                   c(2L,
                     3L,
                     5L,
                     5L,
                     5L,
                     7L,
                     7L,
                     10L,
                     11L,
                     12L,
                     13L,
                     14L,
                     14L,
                     16L,
                     16L,
                     18L))
  
  xpehh.cgu_eut <- ies2xpehh(wgscan.cgu, wgscan.eut)
  cr.xpehh.cgu_eut <-
    calc_candidate_regions(
      xpehh.cgu_eut,
      threshold = 4,
      pval = TRUE,
      overlap = 100000,
      min_n_extr_mrk = 2
    )
  chr_numbers <-
    as.integer(levels(cr.xpehh.cgu_eut$CHR)[cr.xpehh.cgu_eut$CHR])
  expect_identical(chr_numbers, c(2L, 3L, 5L, 5L, 10L, 12L, 13L, 13L, 14L, 14L, 16L))
  
  ## calc_region_stats should yield the same results if applied to the same score
  region_stats.cgu <- calc_region_stats(ihs.cgu,
                                        cr.ihs.cgu,
                                        threshold = 4,
                                        pval = TRUE)
  expect_identical(cr.ihs.cgu, region_stats.cgu)
  
  sink()
  file.remove("test_calc_candidate_regions.log")
})
