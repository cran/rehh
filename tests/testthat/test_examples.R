#checkes various statistics of example 1
context("examples")

test_that("checked_example1", {
  sink(file = "test_example1.log")
  hh <-
    data2haplohh(
      hap_file = "example1.hap",
      map_file = "example1.map",
      allele_coding = "map"
    )
  
  res <- calc_ehh(
    hh,
    mrk = "rs6",
    limehh = 0,
    discard_integration_at_border = FALSE,
    include_nhaplo = TRUE
  )
  expect_equivalent(res$freq, c(0.5, 0.5))
  expected_ehh_ancestral <-
    c(0, 0, 0, 1 / 6, 1 / 2, 1, 1 / 3, 1 / 6, 0, 0, 0)
  expected_ehh_derived <- c(1 / 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 / 3)
  expect_identical(res$ehh$EHH_A, expected_ehh_ancestral)
  expect_identical(res$ehh$EHH_D, expected_ehh_derived)
  expect_identical(res$ehh$NHAPLO_A, c(0L, 0L, rep(4L, 7), 0L, 0L))
  expect_identical(res$ehh$NHAPLO_D, rep(4L, 11))
  expected_ihh_A <-
    sum(expected_ehh_ancestral[2:10] * 10000) + sum(expected_ehh_ancestral[c(1, 11)] *
                                                      5000)
  expected_ihh_D <-
    sum(expected_ehh_derived[2:10] * 10000) + sum(expected_ehh_derived[c(1, 11)] *
                                                    5000)
  expect_equivalent(res$ihh, cbind(expected_ihh_A, expected_ihh_D))
  scan <- scan_hh(hh, discard_integration_at_border = FALSE)
  expected <- data.frame(
    CHR = rep("chr1", 11),
    POSITION = c(
      10000,
      20000,
      30000,
      40000,
      50000,
      60000,
      70000,
      80000,
      90000,
      1e+05,
      110000
    ),
    FREQ_A = c(
      0.875,
      0.875,
      0.75,
      0.75,
      0.875,
      0.5,
      0.75,
      0.875,
      0.75,
      0.875,
      0.5
    )
    ,
    FREQ_D = c(
      0.125,
      0.125,
      0.25,
      0.25,
      0.125,
      0.5,
      0.25,
      0.125,
      0.25,
      0.125,
      0.5
    )
    ,
    NHAPLO_A = c(7L, 7L, 6L, 6L, 7L, 4L, 6L, 7L, 6L, 7L, 4L),
    NHAPLO_D = c(1L, 1L, 2L, 2L, 1L, 4L, 2L, 1L, 2L, 1L, 4L),
    IHH_A = c(
      21992.2619047619,
      36190.4761904762,
      45000,
      47666.6666666667,
      33333.3333333333,
      18816.6666666667,
      45333.3333333333,
      38571.4285714286,
      50666.6666666667,
      35000,
      25075
    ),
    IHH_D = c(
      0,
      0,
      23512.5,
      28025,
      0,
      89166.6666666667,
      28025,
      0,
      23512.5,
      0,
      25833.3333333333
    ),
    IES = c(
      15295.2380952381,
      25892.8571428571,
      22678.5714285714,
      24285.7142857143,
      23750,
      19821.4285714286,
      23035.7142857143,
      27678.5714285714,
      25714.2857142857,
      25000,
      8032.14285714286
    ),
    INES = c(
      21992.2619047619,
      36190.4761904762,
      43437.5,
      46250,
      33333.3333333333,
      52916.6666666667,
      44062.5,
      38571.4285714286,
      48750,
      35000,
      25416.6666666667
    )
  )
  row.names(expected) <- row.names <- c("rs1",
                                        "rs2",
                                        "rs3",
                                        "rs4",
                                        "rs5",
                                        "rs6",
                                        "rs7",
                                        "rs8",
                                        "rs9",
                                        "rs10",
                                        "rs11")
  expect_equal(scan, expected)
  
  # test extract_regions
  regions <- data.frame(
    CHR = c("chr1", "chr2", "chr1", "chr1"),
    START = c(100000, 10000, 50000, 55000),
    END = c(110000, 120000, 61000, 70000)
  )
  scan_regions <- extract_regions(scan, regions)
  
  expect_identical(scan_regions, scan[c(6, 7, 11), ])
  
  # check integration over stepwise constant rehh curve
  ehh <- calc_ehh(
    hh,
    mrk = "rs6",
    discard_integration_at_border = FALSE,
    phased = FALSE,
    limehh = 0.05,
    interpolate = FALSE
  )
  expected_ihh_A <- 25000 - 0.05 * 30000
  expected_ihh_D <- 100000 - 0.05 * 100000
  expect_equivalent(ehh$ihh, cbind(expected_ihh_A, expected_ihh_D))
  
  ehh <- calc_ehh(
    hh,
    mrk = "rs6",
    discard_integration_at_border = FALSE,
    phased = FALSE,
    limehh = 0,
    interpolate = FALSE
  )
  expected_ihh_A <- 25000
  expected_ihh_D <- 100000
  expect_equivalent(ehh$ihh, cbind(expected_ihh_A, expected_ihh_D))
  
  # check that min_maf excludes mrks
  hh <- data2haplohh(
    hap_file = "example1.hap",
    map_file = "example1.map",
    allele_coding = "map",
    min_maf = 0.2
  )
  
  res <- calc_ehh(
    hh,
    mrk = "rs6",
    limehh = 0,
    discard_integration_at_border = FALSE,
    include_nhaplo = TRUE
  )
  expect_equivalent(res$freq, c(0.5, 0.5))
  expected_ehh_ancestral <- c(0, 1 / 3, 1, 1 / 3, 0, 0)
  expected_ehh_derived <- c(1, 1, 1, 1, 1, 1 / 3)
  expect_identical(res$ehh$EHH_A, expected_ehh_ancestral)
  expect_identical(res$ehh$EHH_D, expected_ehh_derived)
  expect_identical(res$ehh$NHAPLO_A, c(rep(4L, 5), 0L))
  expect_identical(res$ehh$NHAPLO_D, rep(4L, 6))
  
  expected_ihh_A <-
    (1 / 3 * 10000 + (1 / 3 + 1) * 20000 + (1 + 1 / 3) * 10000 +
       1 / 3 * 20000) / 2
  expected_ihh_D <-
    ((1 + 1) * 10000 + (1 + 1) * 20000 + (1 + 1) * 30000 + (1 + 1 /
                                                              3) * 20000) / 2
  expect_equivalent(res$ihh, cbind(expected_ihh_A, expected_ihh_D))
  
  # check that scalegap reduces gaps
  ehh <- calc_ehh(
    hh,
    mrk = "rs6",
    limehh = 0,
    scalegap = 15000,
    discard_integration_at_border = FALSE,
    include_nhaplo = TRUE
  )
  expected_ihh_A <-
    (1 / 3 * 10000 + (1 / 3 + 1) * 15000 + (1 + 1 / 3) * 10000 +
       1 / 3 * 15000) / 2
  expected_ihh_D <-
    ((1 + 1) * 10000 + (1 + 1) * 15000 + (1 + 1) * 25000 + (1 + 1 /
                                                              3) * 15000) / 2
  expect_equivalent(ehh$ihh, cbind(expected_ihh_A, expected_ihh_D))
  
  # check that maxgap cuts integration
  ehh <- calc_ehh(
    hh,
    mrk = "rs6",
    limehh = 0,
    discard_integration_at_border = FALSE,
    maxgap = 15000
  )
  expected_ihh_A <- ((1 + 1 / 3) * 10000) / 2
  expected_ihh_D <- ((1 + 1) * 10000) / 2
  expect_equivalent(ehh$ihh, cbind(expected_ihh_A, expected_ihh_D))
  
  sink()
  file.remove("test_example1.log")
})


test_that("checked_example2", {
  sink(file = "test_example2.log")
  hh <- data2haplohh(
    hap_file = "example2.hap",
    map_file = "example2.map",
    allele_coding = "01",
    min_perc_geno.hap = 0,
    min_perc_geno.mrk = 1
  )
  
  res <- calc_ehh(
    hh,
    mrk = "rs6",
    limehh = 0,
    discard_integration_at_border = FALSE,
    include_nhaplo = TRUE
  )
  expect_equivalent(res$freq, c(4 / 11, 3 / 11, 4 / 11))
  expected_ehh_ancestral <-
    c(0, 0, 0, 1 / 6, 1 / 2, 1, 1 / 3, 1 / 6, 0, 0, 0)
  expected_ehh_derived1 <- c(rep(1, 10), 0)
  expected_ehh_derived2 <-
    c(0, 0, 1 / 3, 1 / 3, 1, 1, 1, 1, 0, 0, 0)
  expect_identical(res$ehh$EHH_A, expected_ehh_ancestral)
  expect_identical(res$ehh$EHH_D1, expected_ehh_derived1)
  expect_identical(res$ehh$EHH_D2, expected_ehh_derived2)
  expect_identical(res$ehh$NHAPLO_A, c(0L, 0L, 2L, rep(4L, 6), 0L, 0L))
  expect_identical(res$ehh$NHAPLO_D1, c(rep(3L, 9), 2L, 2L))
  expect_identical(res$ehh$NHAPLO_D2, c(0L, 2L, rep(3L, 3), rep(4L, 3), 3L, 0L, 0L))
  expected_ihh_A <-
    sum(expected_ehh_ancestral[2:10] * 10000) + sum(expected_ehh_ancestral[c(1, 11)] *
                                                      5000)
  expected_ihh_D1 <-
    sum(expected_ehh_derived1[2:10] * 10000) + sum(expected_ehh_derived1[c(1, 11)] *
                                                     5000)
  expected_ihh_D2 <-
    sum(expected_ehh_derived2[2:10] * 10000) + sum(expected_ehh_derived2[c(1, 11)] *
                                                     5000)
  expect_equivalent(res$ihh,
                    cbind(expected_ihh_A, expected_ihh_D1, expected_ihh_D2))
  
  ## check scan
  scan <- scan_hh(hh, discard_integration_at_border = FALSE)
  expected <- data.frame(
    CHR = rep("chr1", 11),
    POSITION = c(
      10000,
      20000,
      30000,
      40000,
      50000,
      60000,
      70000,
      80000,
      90000,
      1e+05,
      110000
    ),
    FREQ_A = c(
      0.916666666666667,
      0.909090909090909,
      0.9,
      0.75,
      0.909090909090909,
      0.363636363636364,
      0.833333333333333,
      0.916666666666667,
      0.7,
      0.777777777777778,
      0.666666666666667
    )
    ,
    FREQ_D = c(
      0.0833333333333333,
      0.0909090909090909,
      0.1,
      0.25,
      0.0909090909090909,
      0.363636363636364,
      0.166666666666667,
      0.0833333333333333,
      0.2,
      0.222222222222222,
      0.333333333333333
    ),
    NHAPLO_A = c(11L, 10L, 9L, 9L, 10L, 4L, 10L, 11L, 7L, 7L, 8L),
    NHAPLO_D = c(1L, 1L, 1L, 3L, 1L, 4L, 2L, 1L, 2L, 2L, 4L)
    ,
    IHH_A = c(
      25045.2380952381,
      32750.9920634921,
      39881.9444444445,
      38453.373015873,
      29680.1587301587,
      18816.6666666667,
      28960.5158730159,
      27370.1298701299,
      40142.8571428571,
      24452.380952381,
      14291.6666666667
    ),
    IHH_D = c(
      0,
      0,
      0,
      21383.3333333333,
      0,
      43216.6666666667,
      28025,
      0,
      33012.5,
      9025,
      13037.5
    ),
    IES = c(
      19362.1212121212,
      25023.5930735931,
      30057.1428571429,
      23654.6717171717,
      22727.2005772006,
      9582.1608946609,
      17869.696969697,
      21479.797979798,
      15730.1587301587,
      11326.9841269841,
      5890.83694083694
    ),
    INES = c(
      25045.2380952381,
      32750.9920634921,
      39881.9444444445,
      38145.5069124424,
      29680.1587301587,
      49005.9523809524,
      28743.6766567576,
      27370.1298701299,
      39857.9545454545,
      23125,
      14158.2264957265
    )
  )
  row.names(expected) <- row.names <- c("rs1",
                                        "rs2",
                                        "rs3",
                                        "rs4",
                                        "rs5",
                                        "rs6",
                                        "rs7",
                                        "rs8",
                                        "rs9",
                                        "rs10",
                                        "rs11")
  expect_equal(scan, expected)
  
  sink()
  file.remove("test_example2.log")
})

test_that("checked_unpolarized", {
  sink(file = "test_unpolarized.log")
  hh <- data2haplohh(hap_file = "example1.hap",
                     map_file = "example1.map", 
                     allele_coding = "map")
  
  ihh_1 <- scan_hh(hh,
                   discard_integration_at_border = FALSE,
                   polarized = TRUE)
  ihh_2 <- scan_hh(hh,
                   discard_integration_at_border = FALSE,
                   polarized = FALSE)
  
  ihs1 <- ihh2ihs(ihh_1, freqbin = 1)
  ihs2 <- ihh2ihs(ihh_2, freqbin = 1)
  
  expect_identical(ihs1, ihs2)
  
  # test that allele coding is irrelevant for unpolarized data,
  # (as long as the order is not changed, since for equal frequency,
  # the "major" allele is the one with lower encoding
  hh@haplo[hh@haplo == 0] <- 4L # arbitrary integer number
  hh@haplo[hh@haplo == 1] <- 7L # arbitrary, bigger integer number
  
  ihh_3 <- scan_hh(hh,
                   discard_integration_at_border = FALSE,
                   polarized = FALSE)
  
  expect_identical(ihh_2, ihh_3)
  sink()
  file.remove("test_unpolarized.log")
})

test_that("checked_limhaplo", {
  sink(file = "test_limhaplo.log")
  hh <- data2haplohh(
    hap_file = "example2.hap",
    map_file = "example2.map",
    allele_coding = "map",
    min_perc_geno.mrk = 1
  )
  
  expected_ihh <-
    c(IHH_A = 18816.6666666667,
      IHH_D1 = 80512.5,
      IHH_D2 = 43216.6666666667)
  ehh <-
    calc_ehh(
      hh,
      mrk = "rs6",
      limhaplo = 3,
      discard_integration_at_border = FALSE
    )
  expect_equal(ehh$ihh, expected_ihh)
  
  expected_ihh <-
    c(IHH_A = 18816.6666666667,
      IHH_D1 = 0,
      IHH_D2 = 28025)
  ehh <-
    calc_ehh(
      hh,
      mrk = "rs6",
      limhaplo = 4,
      discard_integration_at_border = FALSE
    )
  expect_equal(ehh$ihh, expected_ihh)
  
  sink()
  file.remove("test_limhaplo.log")
})

