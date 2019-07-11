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