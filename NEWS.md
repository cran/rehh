### rehh 3.2.1 (November 03, 2020)

* bug-fix: memory error in function calc_sfs_tests()

### rehh 3.2.0 (October 30, 2020)

* corrections and extensions of the vignettes
* added inset between chromosomes in function manhattanplot()
* added option to highlight specific markers in plot.haplohh()
* function subset() can filter markers on maximum number of alleles
* option to set maximal haplotype extension (in base pairs) in function scan_hh_full()
* added function haplohh2sweepfinder() that extracts allele frequencies which can serve as input for the programs SweepFinder or SweeD
* added function calc_sfs_tests() to calculate Tajima's D and Fay & Wu's H, allowing for data with missing values

### rehh 3.1.2 (July 17, 2020)

* added option to parse vcf files using package 'data.table' in order to avoid reliance on package 'vcfR'

### rehh 3.1.1 (June 19, 2020)

* option to use stepwise continuous EHH(S) curves for integration and plotting
* refactored C-code for calculation of furcation trees
  + bug fix for false rendering of missing values in furcations of unphased data
* bug fix: removed reliance of examples and tests on suggested packages
* updated section on calc_ehh() in main vignette
* added functions calc_pairwise_haplen() and scan_hh_full()

### rehh 3.1.0 (March 03, 2020)

* functions ies2xpehh() and ines2rsb():
  + bug fix for unphased data and markers with no homozygous individuals
* function manhattanplot(): 
  + added option to plot absolute values
  + added option to highlight individual markers
* function manhattanplot() and distribplot(qqplot=TRUE):
  + added option to rasterize plot (reduce number of plotted data points)
* function calc_ehh(), calc_furcation() and plot.haplen():
  + changed default order of alleles from their frequency to internal coding
* all plot functions: 
  + allow user to set plot margins
  + allow user to overwrite some plot parameters such as "xlab" and "ylab"
* new plot function for (small) haplohh-objects
* minor optimization of parallelization
* additional comments about unphased data in the main vignette

### rehh 3.0.1 (July 11, 2019)

* bug-fix: memory allocation in C

### rehh 3.0.0 (July 08, 2019)

* Major refactoring of underlying C code and its interface to R
* The package can now handle multi-allelic markers
* Enabled data input from vcf files and the simulation program 'ms'
* Calculation and visualization functions have been separated
* Added support for handling unphased haplotypes and unpolarized markers
* Output in form of matrices has been replaced by data frames
* Simplified column names in data output (use only capital letters and underscores)
* Improved visualization of furcation trees
* Added visualization of shared haplotype lengths
* Added simple function to delineate candidate regions of selection
* Added options to yield virtually identical results with program hapbin
* Made all plots more customizable
* Added new vignette with in-depth discussion of example files
* Disagreement between code and vignette about uni-lateral p-value for Rsb/XP-EHH has been cleared

### rehh 2.0.4 (February 13, 2019)

* Corrected implementation of maxgap in combination with discard_integration_at_border=FALSE
* Added parameter scalegap
* Included section about gaps in vignette, added details in comparison with hapbin, corrected rendering
* Corrected paper title in CITATIONS

### rehh 2.0.3 (June 12, 2017)

* Bug corrected in the ihh2ihsplot (it lead to a crash when the freqbin argument was set to 0)
* Inactivation of graphical display in multiple windows (i.e., plot.new calls) in the calc_ehh, calc_ehhs, ihsplot, rsbplot, xpehhplot and distribplot functions
* Minor corrections in the manual pages

### rehh 2.0.2 (November 15, 2016)

* Option maxgap modified (functions return NA value if not satisfied)
* Example data sets removed from the rehh package and included in the newly developed rehh.data package. The function make.example.files() was modified accordingly.

### rehh 2.0.1 (October 24, 2016)

* Option maxgap added to the calc_ehh, calc_ehhs and scan_hh functions

### rehh 2.0.0 (July 30, 2016)

* Major modification of the algorithm to explore haplotype variability (more than one order of magnitude faster), including parallelization
* Major modification of the data2haplohh function to improve reading efficiency and allele recoding. In addition, a new input data file format (haplotype.in.columns) is now available
* Computation of xp-EHH
* A vignette that details how to use the package is now included in the package
* Several other minor modifications in other functions

### rehh 1.13 (May 13, 2015)

* Removing smartlegend() calls (deprecated in new version of gplot) in calc_ehh.R and fistribplot.R

### rehh 1.12 (April 7, 2015)

* Minor modification of the CITATION file

### rehh 1.11 (August 25, 2013)

* Minor modification of the distribplot.R file (call to plot.density was disabled)

### rehh 1.1 (April 25, 2013)

* Minor modification of the rehh-package.Rd file

* Printing of the advance of each site out of the total SNPs scanned when using scan_hh() function was disabled to save time and simplify screen output (particularly relevant for large data sets).

* Minor bugs corrected in the ehh_utils.c and r_scan_hh.c codes. These bugs had no effect on the results.

### rehh 1.0 (March 15, 2012)