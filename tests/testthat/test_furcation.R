#compares iHH calculated from file with pre-calculated values in dataset
#compares iHH and iES from the functions calc_ehh(s) with scan_hh
context("furcation")

test_that("checked_furcation", {
  # compare calculation with pre-calculated values in object wgscan.cgu
  sink(file = "test_furcation.log")
  
  # example 1
  hh <- data2haplohh(
    hap_file = "example1.hap",
    map_file = "example1.map",
    allele_coding = "01",
    verbose = FALSE
  )
  
  f <- calc_furcation(hh, mrk = "rs6")
  
  expect_identical(
    as.newick(f, allele = 0, side = "left"),
    "(((4:0,(3:0,1:0):10000):10000,2:0):10000);"
  )
  expect_identical(
    as.newick(f, allele = 0, side = "right"),
    "(((1:0,4:0):10000,(3:0,2:0):20000):10000);"
  )
  expect_identical(as.newick(f, allele = 1, side = "left"),
                   "((5/6/8:0,7:0):50000);")
  expect_identical(as.newick(f, allele = 1, side = "right"),
                   "((7/8:0,5/6:0):50000);")
  
  # example 2
  hh <- data2haplohh(
    hap_file = "example2.hap",
    map_file = "example2.map",
    allele_coding = "01",
    verbose = FALSE,
    min_perc_geno.mrk = 50
  )
  
  f <- calc_furcation(hh, mrk = "rs6")
  
  expected_f <- dget("serialized_furcation_example2_rs6_phased.txt")
  expect_equal(f, expected_f, check.names = FALSE)
  
  expect_identical(
    as.newick(f, allele = 0, side = "left"),
    "(((4:0,(1/3:20000):10000):10000,2:0):10000);"
  )
  expect_identical(
    as.newick(f, allele = 0, side = "right"),
    "(((1:0,4:0):10000,(3:0,2:0):20000):10000);"
  )
  expect_identical(as.newick(f, allele = 1, side = "left"),
                   "(5/6/8:50000);")
  expect_identical(as.newick(f, allele = 1, side = "right"),
                   "(((8:0,5:0):10000,6:0):40000);")
  expect_identical(
    as.newick(f, allele = 2, side = "left"),
    "((((9:0,11:0):20000,10:0):10000,12:0):10000);"
  )
  expect_identical(as.newick(f, allele = 2, side = "right"),
                   "((12:0,11:0,9:0,10:0):30000);")
  
  # example 2, unphased
  f <- calc_furcation(hh, mrk = "rs6", phased = FALSE)
  expected_f <-
    dget("serialized_furcation_example2_rs6_unphased.txt")
  
  expect_equal(f, expected_f, check.names = FALSE)
  
  expect_identical(as.newick(f, allele = 0, side = "left"),
                   "(((1:0,2:0):10000,(4:0,3:0):20000):0);")
  expect_identical(as.newick(f, allele = 0, side = "right"),
                   "(((1:0,2:0):10000,(4:0,3:0):10000):0);")
  expect_identical(as.newick(f, allele = 1, side = "left"),
                   "((5/6:50000):0);")
  expect_identical(as.newick(f, allele = 1, side = "right"),
                   "(((5:0,6:0):40000):0);")
  expect_identical(
    as.newick(f, allele = 2, side = "left"),
    "(((9:0,10:0):20000,(11:0,12:0):10000):0);"
  )
  expect_identical(
    as.newick(f, allele = 2, side = "right"),
    "(((9:0,10:0):30000,(12:0,11:0):30000):0);"
  )
  
  # bta data
  hh <- data2haplohh(
    hap_file = "bta12_cgu.hap",
    map_file = "map.inp",
    chr.name = "12",
    allele_coding = "map",
    verbose = FALSE
  )
  
  f <- calc_furcation(hh, mrk = "F1205400")
  expected_f <- dget("serialized_furcation_F1205400.txt")
  
  expect_equal(f, expected_f, check.names = FALSE)
  
  expected_newick <- readLines("furcation_F1205400.newick")[-1]
  expect_identical(as.newick(f, allele = 0, side = "left"), expected_newick[1])
  expect_identical(as.newick(f, allele = 0, side = "right"), expected_newick[2])
  expect_identical(as.newick(f, allele = 1, side = "left"), expected_newick[3])
  expect_identical(as.newick(f, allele = 1, side = "right"), expected_newick[4])
  
  ## haplen calculates longest shared haplotype per chromosome,
  ## for left and right side independently
  h <- calc_haplen(f)
  lengths_haplen_left <- positions(hh)["F1205400"] - h$haplen$MIN
  lengths_haplen_right <- h$haplen$MAX - positions(hh)["F1205400"]
  
  ## pairwise haplen calculates length of all mutually shared haplotypes
  ph <- calc_pairwise_haplen(hh, mrk = "F1205400", side = "left")
  lengths_longest_pairwise_haplen_left <- apply(ph, 1, max)
  expect_equivalent(lengths_haplen_left, lengths_longest_pairwise_haplen_left)
  
  ph <- calc_pairwise_haplen(hh, mrk = "F1205400", side = "right")
  lengths_longest_pairwise_haplen_right <- apply(ph, 1, max)
  expect_equivalent(lengths_haplen_right,
                    lengths_longest_pairwise_haplen_right)
  
  sink()
  file.remove("test_furcation.log")
})
