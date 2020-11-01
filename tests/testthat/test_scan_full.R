#checks equivalence between scan_hh_full and scan_hh
#and scan_hh_full and calc_pairwise_haplen
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
  
  #compare with element wise calculation
  IES <- vapply(1:nmrk(hh), function(x) {
    sum(calc_pairwise_haplen(hh, mrk = x)) / (nhap(hh) * (nhap(hh) - 1))
  }, FUN.VALUE = 0.0)
  
  expect_equal(scan2$IES, IES)
  
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
  
  #compare with element wise calculation
  IES <- vapply(1:nmrk(hh), function(x) {
    sum(calc_pairwise_haplen(hh, mrk = x, phased = FALSE)) / nhap(hh)
  }, FUN.VALUE = 0.0)
  #for unphased data IES is set to INES (normalized by focal marker)
  INES <- IES * nhap(hh) / (scan2$NHAPLO_A + scan2$NHAPLO_D)
  expect_equal(scan2$INES, INES)
  
  #test max_extend and maxgap
  hh <- data2haplohh(
    hap_file = "example1.hap",
    map_file = "example1.map",
    allele_coding = "01",
    verbose = FALSE
  )
  
  #test max_extend
  scan1 <-
    scan_hh_full(hh,
                 discard_integration_at_border = FALSE,
                 max_extend = 20000)
  
  IHH_A <-
    c(
      17142.8571428571,
      25238.0952380952,
      33333.3333333333,
      33333.3333333333,
      29523.8095238095,
      28333.3333333333,
      31333.3333333333,
      32380.9523809524,
      36666.6666666667,
      25238.0952380952,
      20000
    )
  IHH_D <- c(0, 0, 30, 30, 0, 40, 40, 0, 30, 0, 15) * 1000
  
  expect_equal(scan1$IHH_A, IHH_A)
  expect_equal(scan1$IHH_D, IHH_D)
  
  
  #test maxgap
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
