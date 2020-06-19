#compares iHH calculated from file with pre-calculated values in dataset
#compares iHH and iES from the functions calc_ehh(s) with scan_hh
context("scan_hh_full")

test_that("checked_scan_hh_full", {
  # compare calculation of scan_hh and scan_hh2
  sink(file = "test_scan_hh_full.log")
  hh <- data2haplohh(
    hap_file = "bta12_cgu.hap",
    map_file = "map.inp",
    allele_coding = "map",
    chr.name = 12,
    verbose = FALSE
  )
  
  #phased
  scan1 <- scan_hh(
    hh,
    interpolate = FALSE,
    limehh = 0,
    limehhs = 0,
    phased = TRUE,
    discard_integration_at_border = FALSE
  )
  
  scan2 <- scan_hh_full(hh,
                        phased = TRUE,
                        discard_integration_at_border = FALSE)
  
  expect_equal(scan1, scan2)
  
  #unphased
  scan1 <- scan_hh(
    hh,
    interpolate = FALSE,
    limehh = 0,
    limehhs = 0,
    phased = FALSE,
    discard_integration_at_border = FALSE
  )
  
  scan2 <- scan_hh_full(hh,
                        phased = FALSE,
                        discard_integration_at_border = FALSE)
  
  expect_equal(scan1, scan2)
  
  #test maxgap
  hh <- data2haplohh(
    hap_file = "example1.hap",
    map_file = "example1.map",
    allele_coding = "01",
    verbose = FALSE
  )
  hh@positions[5] <- 55000
  hh@positions[11] <- 115000
  for (maxgap in c(4000, 5000, 10000, 15000, NA)) {
    for (phased in c(TRUE, FALSE)) {
      scan1 <- scan_hh(
        hh,
        interpolate = FALSE,
        limehh = 0,
        limehhs = 0,
        discard_integration_at_border = FALSE,
        maxgap = maxgap,
        phased = phased
      )
      scan2 <- scan_hh_full(
        hh,
        discard_integration_at_border = FALSE,
        maxgap = maxgap,
        phased = phased,
        geometric.mean = FALSE
      )
      
      expect_equal(scan1, scan2)
    }
  }
  
  sink()
  
  file.remove("test_scan_hh_full.log")
})
