#checks error messages
context("error messages")

setup({
  sink(file = "test_error_messages.log")
})

teardown({
  sink()
  file.remove("test_error_messages.log")
})

test_that("errors_parameter", {
  tmp1 <- tempfile()
  writeLines(c("H1 0"), tmp1)
  tmp2 <- tempfile()
  writeLines(c("rs1 1 10000"), tmp2)
  expect_error(
    data2haplohh(tmp1, tmp2, min_perc_geno.hap = 101),
    "min_perc_geno.hap should lie in the interval"
  )
  expect_error(
    data2haplohh(tmp1, tmp2, min_perc_geno.mrk = 0),
    "min_perc_geno.mrk should lie in the interval"
  )
  expect_error(data2haplohh(tmp1, tmp2, min_maf = 0.6),
               "min_maf should lie in the interval")
  expect_error(
    data2haplohh(tmp1, tmp2, position_scaling_factor = 0),
    "position_scaling_factor must be a positive real number"
  )
  expect_error(data2haplohh(tmp1),
               "No map file specified.")
})

test_that("errors_map_file", {
  tmp1 <- tempfile()
  # non-unique haplotype identifiers
  writeLines(c("H1 1 1", "H1 1 1"), tmp1)
  
  tmp2 <- tempfile()
  writeLines(c("rs1 1 10000", "rs2 1 20000"), tmp2)
  expect_warning(data2haplohh(tmp1, tmp2),
                 "Haplotype identifiers were not unique")
  
  tmp1 <- tempfile()
  writeLines(c("H1 0 1", "H2 1 1"), tmp1)
  expect_error(data2haplohh(tmp1), "No map file specified.")
  
  tmp2 <- tempfile()
  # only one marker in map file
  writeLines(c("rs1 1 10000"), tmp2)
  expect_error(data2haplohh(tmp1, tmp2), "number of markers")
  
  tmp2 <- tempfile()
  writeLines(c("rs1 10000", "rs2 20000"), tmp2)
  # map files contains only two columns
  expect_error(data2haplohh(tmp1, tmp2), "Wrong format for map file.")
  
  tmp2 <- tempfile()
  writeLines(c("rs1 1 20000", "rs2 2 10000"), tmp2)
  # two chromsomes, none specified
  expect_error(data2haplohh(tmp1, tmp2), "specify a chromosome name")
  # wrong chromosome specified
  expect_error(data2haplohh(tmp1, tmp2, chr.name = 3),
               "specify one chromosome")
  
  tmp2 <- tempfile()
  writeLines(c("rs1 1 20000", "rs2 1 10000"), tmp2)
  # positions not ordered
  expect_error(data2haplohh(tmp1, tmp2),
               "Markers must be ordered numerically in the map file.")
  
  tmp2 <- tempfile()
  writeLines(c("rs1 1 10000", "rs2 1 10000"), tmp2)
  # multiple markers with same position
  expect_error(data2haplohh(tmp1, tmp2),
               "1 markers have non-unique positions")
  
  # warning if removed
  expect_warning(
    data2haplohh(tmp1, tmp2, remove_multiple_markers = TRUE),
    "Removed 1 markers with non-unique positions"
  )
})

test_that("errors_hap_file", {
  tmp1 <- tempfile()
  writeLines(c("H1 A G", "H2 C A"), tmp1)
  tmp2 <- tempfile()
  writeLines(c("rs1 1 10000", "rs2 1 20000"), tmp2)
  # wrong coding
  expect_error(
    data2haplohh(tmp1, tmp2, allele_coding = "12"),
    "Alleles are not coded in format \"12\""
  )
  expect_error(
    data2haplohh(tmp1, tmp2, allele_coding = "01"),
    "Alleles are not coded in format \"01\""
  )
  tmp1 <- tempfile()
  writeLines(c("H1 1 0", "H2 -1 2"), tmp1)
  expect_error(
    data2haplohh(tmp1, tmp2, allele_coding = "01"),
    "Found alleles coded by negative numbers"
  )
  
  expect_error(
    data2haplohh("bta12_hapguess_switch.out", "map.inp", chr.name = 12),
    "Please specify by 'popsel' one of the following population numbers"
  )
})

test_that("errors_ms", {
  skip_if_not_installed("gap")
  
  expect_error(data2haplohh("ms.out"),
               "Please select one by specifying its number in 'chr.name'.")
  expect_error(
    data2haplohh("ms.out", chr.name = 0),
    "Please select one by specifying its number in 'chr.name'."
  )
  expect_error(
    data2haplohh("ms.out", chr.name = "chr1"),
    "For ms output files 'chr.name' has to be an integer number"
  )
})


test_that("errors_vcf_vcfR", {
  skip_if_not_installed("vcfR")
  tmp <- tempfile()
  ## different ploidy at different markers
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=NS,Number=1,Type=Integer>",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\trs1\tG\tT\t100\tPASS\tNS=2\tGT\t.\t0",
      "chr1\t20000\trs2\tG\tA\t100\tPASS\tNS=2\tGT\t0\t1|0"
    ),
    tmp
  )
  expect_error(
    data2haplohh(tmp, polarize_vcf = FALSE, vcf_reader = "vcfR"),
    "1 individuals have different ploidy at different markers"
  )
  
  ## absent ancestral alleles
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=NS,Number=1,Type=Integer>",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\trs1\tG\tT\t100\tPASS\tNS=2\tGT\t.\t0|0",
      "chr1\t20000\trs2\tG\tA\t100\tPASS\tNS=2\tGT\t0\t1|0"
    ),
    tmp
  )
  
  expect_error(data2haplohh(tmp, vcf_reader = "vcfR"),
               "Key 'AA' not found")

  ## no polarization possible
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\trs1\tG\tT\t100\tPASS\tAA=C\tGT\t.\t0|0",
      "chr1\t20000\trs2\tG\tA\t100\tPASS\tAA=C\tGT\t0\t1|0"
    ),
    tmp
  )
  
  expect_error(data2haplohh(tmp, vcf_reader = "vcfR"),
               "No marker could be polarized")
  
  ## no marker identifiers
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=NS,Number=1,Type=Integer>",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\t.\tG\tT\t100\tPASS\tNS=2\tGT\t.\t0|0",
      "chr1\t20000\t.\tG\tA\t100\tPASS\tNS=2\tGT\t0\t1|0"
    ),
    tmp
  )
  expect_output(
    data2haplohh(tmp, polarize_vcf = FALSE, vcf_reader = "vcfR"),
    "No marker identifiers found in vcf file"
  )
  
  ## non unique identifiers
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=NS,Number=1,Type=Integer>",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\trs1\tG\tT\t100\tPASS\tNS=2\tGT\t.\t0|0",
      "chr1\t20000\trs1\tG\tA\t100\tPASS\tNS=2\tGT\t0\t1|0"
    ),
    tmp
  )
  expect_error(
    data2haplohh(tmp, polarize_vcf = FALSE, vcf_reader = "vcfR"),
    "ID column contains non-unique names"
  )
  
  unlink(tmp)
  
})

test_that("errors_vcf_data.frame", {
  skip_if_not_installed("data.table")
  tmp <- tempfile()
  ## different ploidy at different markers
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=NS,Number=1,Type=Integer>",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\trs1\tG\tT\t100\tPASS\tNS=2\tGT\t.\t0",
      "chr1\t20000\trs2\tG\tA\t100\tPASS\tNS=2\tGT\t0\t1|0"
    ),
    tmp
  )
  expect_error(
    data2haplohh(tmp, polarize_vcf = FALSE, vcf_reader = "data.table"),
    "1 individuals have different ploidy at different markers"
  )
  
  ## no polarization possible
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\trs1\tG\tT\t100\tPASS\tAA=C\tGT\t.\t0|0",
      "chr1\t20000\trs2\tG\tA\t100\tPASS\tAA=C\tGT\t0\t1|0"
    ),
    tmp
  )
  
  expect_error(data2haplohh(tmp, vcf_reader = "data.table"),
               "No marker could be polarized")
  
  ## no marker identifiers
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=NS,Number=1,Type=Integer>",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\t.\tG\tT\t100\tPASS\tNS=2\tGT\t.\t0|0",
      "chr1\t20000\t.\tG\tA\t100\tPASS\tNS=2\tGT\t0\t1|0"
    ),
    tmp
  )
  expect_output(
    data2haplohh(tmp, polarize_vcf = FALSE, vcf_reader = "data.table"),
    "No marker identifiers found in vcf file"
  )
  

  ## non unique identifiers
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=NS,Number=1,Type=Integer>",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2",
      "chr1\t10000\trs1\tG\tT\t100\tPASS\tNS=2\tGT\t.\t0|0",
      "chr1\t20000\trs1\tG\tA\t100\tPASS\tNS=2\tGT\t0\t1|0"
    ),
    tmp
  )
  expect_error(
    data2haplohh(tmp, polarize_vcf = FALSE, vcf_reader = "data.table"),
    "ID column contains non-unique names"
  )
  
  unlink(tmp)
  
})
