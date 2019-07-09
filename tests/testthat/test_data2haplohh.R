#compares reading in of the three different hap file formats
context("data2haplohh")

test_that("checked_data2haplohh", {
  sink(file = "test_data2haplohh.log")
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
  hap3 <-
    data2haplohh(
      hap_file = "bta12_hapguess_switch.out",
      map_file = "map.inp",
      allele_coding = "map",
      popsel = 7,
      chr.name = 12
    )
  hap4 <-
    data2haplohh(hap_file = "bta12_cgu.vcf.gz",
                 polarize_vcf = FALSE)

  sink()
  #reading standard format includes row-names, hence not equality
  expect_equivalent(hap1, hap2)
  expect_equivalent(hap2, hap3)
  expect_equal(hap3, hap4)
  
  file.remove("test_data2haplohh.log")
})
