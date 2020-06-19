#compares reading in of the three different hap file formats
context("data2haplohh")

test_that("checked_example1", {
  sink(file = "test_example1.log")
  
  hh1 <-
    data2haplohh(hap_file = "example1.hap",
                 map_file = "example1.map",
                 allele_coding = "01")
  hh2 <-
    data2haplohh(hap_file = "example1.hap",
                 map_file = "example1.map",
                 allele_coding = "map")
  expect_identical(hh1, hh2)
  
  hh3 <-
    data2haplohh(hap_file = "example1.hap",
                 map_file = "example1.map",
                 allele_coding = "none")
  expect_identical(hh1, hh3)
  
  skip_if_not_installed("vcfR")
  
  hh2 <- data2haplohh(hap_file = "example1.vcf")
  
  expect_identical(hh1, hh2)
  sink()
  
  file.remove("test_example1.log")
})

test_that("checked_example2", {
  sink(file = "test_example2.log")
  
  hh1 <-
    data2haplohh(hap_file = "example2.hap",
                 map_file = "example2.map",
                 allele_coding = "01")
  hh2 <-
    data2haplohh(hap_file = "example2.hap",
                 map_file = "example2.map",
                 allele_coding = "map")
  expect_identical(hh1, hh2)
  
  skip_if_not_installed("vcfR")
  
  hh2 <- data2haplohh(hap_file = "example2.vcf")
  
  expect_identical(hh1, hh2)
  
  ## test with filtering
  hh1 <- data2haplohh(
    hap_file = "example2.hap",
    map_file = "example2.map",
    allele_coding = "01",
    min_perc_geno.hap = 50,
    min_perc_geno.mrk = 50
  )
  hh2 <-
    data2haplohh(
      hap_file = "example2.vcf",
      min_perc_geno.hap = 50,
      min_perc_geno.mrk = 50
    )
  
  expect_equivalent(hh1, hh2)
  
  sink()
  
  file.remove("test_example2.log")
})

test_that("checked_bta", {
  sink(file = "test_bta.log")
  
  hap1 <- data2haplohh(
    hap_file = "bta12_cgu.hap",
    map_file = "map.inp",
    allele_coding = "map",
    chr.name = 12
  )
  
  hap2 <-
    data2haplohh(
      hap_file = "bta12_cgu.thap",
      map_file = "map.inp",
      haplotype.in.columns = TRUE,
      allele_coding = "map",
      chr.name = 12
    )
  #reading standard format includes row-names, reading
  #transposed format does not, hence only equivalent, but not equal
  expect_equivalent(hap1, hap2)
  
  hap3 <-
    data2haplohh(
      hap_file = "bta12_hapguess_switch.out",
      map_file = "map.inp",
      allele_coding = "map",
      popsel = 7,
      chr.name = 12
    )
  
  expect_equivalent(hap1, hap3)
  
  skip_if_not_installed("vcfR")
  
  hap4 <-
    data2haplohh(hap_file = "bta12_cgu.vcf.gz",
                 polarize_vcf = FALSE)
  
  expect_equal(hap3, hap4)
  
  sink()
  
  file.remove("test_bta.log")
})
