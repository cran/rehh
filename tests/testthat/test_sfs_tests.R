#checks sfs statistics of example 1
context("sfs_test_statistics")

test_that("checked_example1", {
  sink(file = "test_sfs_example1.log")
  hh <-
    data2haplohh(hap_file = "example1.hap",
                 map_file = "example1.map",
                 allele_coding = "map")
  
  res <- calc_sfs_tests(hh, window_size = 110000)
  
  expect_equal(res$THETA_S, 4.24242424242424)
  expect_equal(res$TAJIMA_D,-0.159251163979525)
  expect_equal(res$FAY_WU_H, 0.928627986893831)
  expect_equal(res$ZENG_E,-1.07358846844319)
  
  res <- calc_sfs_tests(hh,
                        window_size = 50000,
                        overlap = 25000,
                        right = FALSE)
  
  expect_equal(
    res$THETA_S,
    c(
      1.54269972451791,
      1.92837465564738,
      1.92837465564738
    )
  )
  expect_equal(
    res$TAJIMA_D,
    c(
      -0.524737350335118,
      0.420384850289377,
      0.000462724105986242
    )
  )
  expect_equal(
    res$FAY_WU_H,
    c(
      0.898992052528498,
      0.822978681944425,
      0.768113436481463
    )
  )
  expect_equal(
    res$ZENG_E,
    c(
      -1.16431403321423,
      -0.527192758957502,
      -0.738186112321971
    )
  )
  
  sink()
  file.remove("test_sfs_example1.log")
})
