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
  
  skip_if_not_installed("data.table")
  skip_if_not_installed("R.utils")
  
  hap4 <-
    data2haplohh(hap_file = "bta12_cgu.vcf.gz",
                 polarize_vcf = FALSE,
                 vcf_reader = "data.table")
  
  expect_equal(hap3, hap4)
  
  skip_if_not_installed("vcfR")
  
  hap5 <-
    data2haplohh(hap_file = "bta12_cgu.vcf.gz",
                 polarize_vcf = FALSE,
                 vcf_reader = "vcfR")
  
  expect_equal(hap3, hap5)
  
  sink()
  
  file.remove("test_bta.log")
})

test_that("data2haplohh_AA_scan", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("vcfR")
  sink(file = "test_AA_scan.log")
  
  tmp <- tempfile()
  ## check if ancestral allele is correctly parsed
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles\">",
      "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">",
      "##INFO=<ID=VT,Number=.,Type=String,Description=\"Variant Type\">",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\trs1\tG\tT\t100\tPASS\tAA=G\tGT\t.\t0|0",
      "chr1\t20000\trs2\tG\tA\t100\tPASS\t.\tGT\t0\t1|0",
      "chr1\t30000\trs3\tA\tC\t100\tPASS\tAC=2;AF=0.5;AA=C;VT=SNP\tGT\t0\t1|0"
    ),
    tmp
  )
  
  hh1 <-
    suppressWarnings(data2haplohh(tmp, vcf_reader = "data.table", verbose = FALSE))
  hh2 <-
    suppressWarnings(data2haplohh(tmp, vcf_reader = "vcfR", verbose = FALSE))
  
  expect_equal(hh1, hh2)
  
  unlink(tmp)
  sink()
  file.remove("test_AA_scan.log")
})
