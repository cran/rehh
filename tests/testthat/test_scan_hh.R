#compares iHH calculated from file with pre-calculated values in dataset
#compares iHH and iES from the functions calc_ehh(s) with scan_hh
context("scan_hh")

test_that("checked_scan_hh", {
  skip_if_not_installed("rehh.data")
  # compare calculation with pre-calculated values in object wgscan.cgu
  sink(file = "test_scan_hh.log")
  hh <- data2haplohh(
    hap_file = "bta12_cgu.hap",
    map_file = "map.inp",
    allele_coding = "map",
    chr.name = 12
  )
  sink()
  
  scan <- scan_hh(hh, discard_integration_at_border = FALSE)
  
  requireNamespace("rehh.data", quietly = TRUE)
  data(wgscan.cgu)
  
  expect_equal(scan$IHH_A,
               wgscan.cgu[wgscan.cgu$CHR == 12, ]$iHH_A,
               tolerance = 0.02,
               scale = NULL)
  expect_equal(scan$IHH_D,
               wgscan.cgu[wgscan.cgu$CHR == 12, ]$iHH_D,
               tolerance = 0.01,
               scale = NULL)
  expect_equal(scan$INES, wgscan.cgu[wgscan.cgu$CHR == 12, ]$iES_Tang_et_al_2007)
  file.remove("test_scan_hh.log")
})

test_that("checked_calc_ehh(s)", {
  # compare calculation of calc_ehh(s) with scan_hh
  sink(file = "test_calc_ehh.log")
  hh <- data2haplohh(
    hap_file = "bta12_cgu.hap",
    map_file = "map.inp",
    allele_coding = "map",
    chr.name = 12
  )
  sink()
  # standard scan
  scan <- scan_hh(hh)
  
  # slow scan
  # create empty vectors of size nmrk
  NHAPLO_A <- NHAPLO_D <- vector("integer", nmrk(hh))
  FREQ_A <-
    FREQ_D <-
    IHH_A <- IHH_D <- IES <- INES <- vector("numeric", nmrk(hh))
  # invoke calc_ehh and calc_ehhs for each marker
  for (i in 1:nmrk(hh)) {
    ehh <- calc_ehh(hh, mrk = i, include_nhaplo = TRUE)
    FREQ_A[i] <- ehh$freq["FREQ_A"]
    FREQ_D[i] <- ehh$freq["FREQ_D"]
    NHAPLO_A[i] <- ehh$ehh[mrk.names(hh)[i], "NHAPLO_A"]
    NHAPLO_D[i] <- ehh$ehh[mrk.names(hh)[i], "NHAPLO_D"]
    IHH_A[i] <- ehh$ihh["IHH_A"]
    IHH_D[i] <- ehh$ihh["IHH_D"]
    
    ehhs <- calc_ehhs(hh, mrk = i)
    IES[i] <- ehhs$IES
    INES[i] <- ehhs$INES
  }
  # create data frame
  slow_scan <- data.frame(
    FREQ_A = FREQ_A,
    FREQ_D = FREQ_D,
    NHAPLO_A = NHAPLO_A,
    NHAPLO_D = NHAPLO_D,
    IHH_A = IHH_A,
    IHH_D = IHH_D,
    IES = IES,
    INES = INES,
    row.names = mrk.names(hh)
  )
  expect_identical(scan[3:10], slow_scan)
  file.remove("test_calc_ehh.log")
})
