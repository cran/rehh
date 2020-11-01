## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(comment = ">",fig.height = 4.5, fig.width = 4.5, fig.show = "hold")

## ----library, message = FALSE-------------------------------------------------
library(rehh)

## ----make_examples, results = 'hide'------------------------------------------
make.example.files()

## ----minimalcodeexample, fig.align = 'center', results = "hide"---------------
hh <-                  # data input
  data2haplohh(
    hap_file = "bta12_cgu.hap",
    map_file = "map.inp",
    chr.name = "12",
    allele_coding = "map"
  )
scan <- scan_hh(hh)    # calculation of EHH and integration
                       # (combine results from different chromosomes)
ihs <- ihh2ihs(scan)   # log ratio for alleles and standardization
manhattanplot(ihs)     # plot of the statistics

## ----map_inp------------------------------------------------------------------
# show first 6 lines of file
cat(readLines("map.inp", n = 6), sep = "\n")

## ----eval = FALSE-------------------------------------------------------------
#  ?data2haplohh

## ----example1-----------------------------------------------------------------
hh <- data2haplohh(hap_file = "bta12_cgu.hap",
                   map_file = "map.inp",
                   chr.name = 12,
                   allele_coding = "map")

## ----error = TRUE-------------------------------------------------------------
hh <- data2haplohh(hap_file = "bta12_cgu.hap",
                   map_file = "map.inp",
                   allele_coding = "map") 

## ----error = TRUE-------------------------------------------------------------
hh <- data2haplohh(hap_file = "bta12_cgu.hap",
                   map_file = "map.inp",
                   chr.name = 18,
                   allele_coding = "map")

## ----example2-----------------------------------------------------------------
hh <- data2haplohh(hap_file = "bta12_cgu.thap",
                   map_file = "map.inp",
                   chr.name = 12,
                   allele_coding = "map",
                   haplotype.in.columns = TRUE)

## ----example3-----------------------------------------------------------------
hh <- data2haplohh(hap_file = "bta12_hapguess_switch.out",
                   map_file = "map.inp",
                   chr.name = 12,
                   popsel = 7,
                   allele_coding = "map")

## ----error = TRUE-------------------------------------------------------------
hh <- data2haplohh(hap_file = "bta12_hapguess_switch.out",
                   map_file = "map.inp",
                   chr.name = 12,
                   allele_coding = "map")

## ----vcf_example, eval = requireNamespace("data.table", quietly = TRUE) & requireNamespace("R.utils", quietly = TRUE)----
hh <- data2haplohh(hap_file = "bta12_cgu.vcf.gz",
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")

## ----ms_example, eval = requireNamespace("gap", quietly = TRUE)---------------
hh <- data2haplohh(hap_file = "ms.out",
                   chr.name = 2,
                   position_scaling_factor = 1000)

## -----------------------------------------------------------------------------
hh_subset = subset(hh, select.hap = 1:5, min_maf = 0)

## -----------------------------------------------------------------------------
hh_subset = subset(hh, select.mrk = -1)

## ----eval=FALSE---------------------------------------------------------------
#  ?calc_ehh

## -----------------------------------------------------------------------------
#example haplohh object (280 haplotypes, 1424 SNPs) see ?haplohh_cgu_bta12 for details
data(haplohh_cgu_bta12)
#computing EHH statistics for the focal SNP with name "F1205400"
#which displays a strong signal of selection
res <- calc_ehh(haplohh_cgu_bta12, 
                mrk = "F1205400", 
                include_nhaplo = TRUE)

## -----------------------------------------------------------------------------
res$mrk.name

## -----------------------------------------------------------------------------
res$freq

## -----------------------------------------------------------------------------
res$ehh

## -----------------------------------------------------------------------------
res$ihh

## ---- ehh, fig.align = 'center', fig.cap = "Graphical output of the plot.ehh() function", fig.pos = '!h', fig.lp = 'fig:'----
plot(res)

## ----ehh2, fig.align = 'center', fig.cap = "Graphical output of the plot.ehh() function", fig.pos = '!h', fig.lp = 'fig:'----
plot(calc_ehh(haplohh_cgu_bta12, 
              mrk = "F1205400", 
              phased = FALSE),
     xlim = c(2.55E7, 3.05E7))

## ---- eval=FALSE--------------------------------------------------------------
#  ?calc_ehhs

## -----------------------------------------------------------------------------
data(haplohh_cgu_bta12)
res <- calc_ehhs(haplohh_cgu_bta12, 
                 mrk = "F1205400", 
                 include_nhaplo = TRUE)

## -----------------------------------------------------------------------------
res$mrk.name

## -----------------------------------------------------------------------------
res$ehhs 

## -----------------------------------------------------------------------------
res$IES

## -----------------------------------------------------------------------------
res$INES

## ----ehhs, fig.align = 'center', fig.cap = 'Graphical output of the plot.ehhs() function', fig.pos = '!h', fig.lp = 'fig:'----
plot(res)

## ----nehhs, fig.align = 'center', fig.cap = 'Graphical output of the plot.ehhs() function', fig.pos = '!h', fig.lp = 'fig:'----
plot(res, nehhs = TRUE)

## ----ehhs2, fig.align = 'center', fig.cap = 'Graphical output of the plot.ehhs() function', fig.pos = '!h', fig.lp = 'fig:'----
plot(calc_ehhs(haplohh_cgu_bta12,
               mrk = "F1205400",
               phased = FALSE))

## ----eval=FALSE---------------------------------------------------------------
#  ?scan_hh

## -----------------------------------------------------------------------------
data(haplohh_cgu_bta12)
scan <- scan_hh(haplohh_cgu_bta12)

## -----------------------------------------------------------------------------
scan[453:459,]

## -----------------------------------------------------------------------------
# perform scan using scan_hh
system.time(scan <- scan_hh(haplohh_cgu_bta12))

## ----slow_scan----------------------------------------------------------------
# perform scan applying calc_ehh and calc_ehhs to each marker
slow_scan_hh <- function(haplohh) {
  # create empty vectors of size nmrk
  IHH_A <- IHH_D <- IES <- INES <- vector("numeric", nmrk(haplohh))
  # invoke calc_ehh and calc_ehhs for each marker
  for (i in 1:nmrk(haplohh)) {
    res <- calc_ehh(haplohh, mrk = i)
    IHH_A[i] <- res$ihh["IHH_A"]
    IHH_D[i] <- res$ihh["IHH_D"]
    res <- calc_ehhs(haplohh, mrk = i)
    IES[i] <- res$IES
    INES[i] <- res$INES
  }
  # create data frame (the return value of this function)
  data.frame(IHH_A = IHH_A, 
             IHH_D = IHH_D,
             IES = IES,
             INES = INES)
}
system.time(slow_scan <- slow_scan_hh(haplohh_cgu_bta12))

## ----scancomp-----------------------------------------------------------------
identical(slow_scan[, "IHH_A"], scan[, "IHH_A"])
identical(slow_scan[, "IHH_D"], scan[, "IHH_D"])
identical(slow_scan[, "IES"], scan[, "IES"])
identical(slow_scan[, "INES"], scan[, "INES"])

## ----ihh2ihs1, eval = FALSE---------------------------------------------------
#  ## demo code - no data files for all chromosomes provided
#  for(i in 1:29) {
#    # haplotype file name for each chromosome
#    hap_file = paste("hap_chr_", i, ".cgu", sep = "")
#    # create internal representation
#    hh <- data2haplohh(hap_file = hap_file,
#                       map_file = "map.inp",
#                       chr.name = i,
#                       allele_coding = "map")
#    # perform scan on a single chromosome (calculate iHH values)
#    scan <- scan_hh(hh)
#    # concatenate chromosome-wise data frames to
#    # a data frame for the whole genome
#    # (more efficient ways certainly exist...)
#    if (i == 1) {
#      wgscan <- scan
#    } else {
#      wgscan <- rbind(wgscan, scan)
#    }
#  }
#  # calculate genome-wide iHS values
#  wgscan.ihs <- ihh2ihs(wgscan)

## ----ihh2ihs2-----------------------------------------------------------------
library(rehh.data)
data(wgscan.cgu)
wgscan.ihs.cgu <- ihh2ihs(wgscan.cgu)

## -----------------------------------------------------------------------------
head(wgscan.ihs.cgu$ihs)

## -----------------------------------------------------------------------------
head(wgscan.ihs.cgu$frequency.class)

## ----ines2rsb2----------------------------------------------------------------
data(wgscan.cgu) ; data(wgscan.eut)
## results from a genome scan (44,057 SNPs) see ?wgscan.eut and ?wgscan.cgu for details
rsb.cgu_eut <- ines2rsb(scan_pop1 = wgscan.cgu,
                        scan_pop2 = wgscan.eut,
                        popname1 = "CGU",
                        popname2 = "EUT")

## -----------------------------------------------------------------------------
head(rsb.cgu_eut)

## ----ies2xpehh2---------------------------------------------------------------
data(wgscan.cgu) ; data(wgscan.eut)
## results from a genome scan (44,057 SNPs) see ?wgscan.eut and ?wgscan.cgu for details
xpehh.cgu_eut <- ies2xpehh(scan_pop1 =  wgscan.cgu,
                           scan_pop2 =  wgscan.eut,
                           popname1 = "CGU",
                           popname2 = "EUT")

## -----------------------------------------------------------------------------
head(xpehh.cgu_eut)

## ----candidate_regions--------------------------------------------------------
cr.cgu <- calc_candidate_regions(wgscan.ihs.cgu,
                                 threshold = 4,
                                 pval = TRUE,
                                 window_size = 1E6,
                                 overlap = 1E5,
                                 min_n_extr_mrk = 2)
cr.cgu

## ----freqbin, fig.align='center', fig.lp='fig:', fig.cap='Graphical output of the freqbinplot() function', fig.pos='!h'----
freqbinplot(wgscan.ihs.cgu)

## ----comp, echo = TRUE, fig.align = 'center', fig.lp = 'fig:', fig.cap = 'Comparison of Rsb and XP-EHH values across the CGU and EUT populations', fig.pos = "!h"----
plot(rsb.cgu_eut[, "RSB_CGU_EUT"],
     xpehh.cgu_eut[, "XPEHH_CGU_EUT"],
     xlab = "Rsb",
     ylab = "XP-EHH",
     pch = ".",
     xlim = c(-7.5, 7.5),
     ylim = c(-7.5, 7.5))
# add red circle for marker with name "F1205400"
points(rsb.cgu_eut["F1205400", "RSB_CGU_EUT"],
       xpehh.cgu_eut["F1205400", "XPEHH_CGU_EUT"],
       col = "red")
# add dashed diagonal
abline(a = 0, b = 1, lty = 2)

## ----distribplot, fig.align = 'center', fig.lp = 'fig:', fig.cap = 'Graphical output of the distribplot() function', fig.pos="!h"----
distribplot(wgscan.ihs.cgu$ihs$IHS, xlab = "iHS")

## ----qqplot, fig.align = 'center', fig.lp='fig:', fig.cap = 'Graphical output of the distribplot() function', fig.pos = "!h"----
distribplot(wgscan.ihs.cgu$ihs$IHS, 
            xlab = "iHS", 
            qqplot = TRUE)

## ----manhattanplot, fig.align = 'center', fig.width = 7, fig.lp = 'fig:', fig.cap = 'Graphical output of the manhattanplot() function', fig.pos = '!h'----
manhattanplot(wgscan.ihs.cgu,
              main = "iHS (CGU cattle breed)")

## ----pval, fig.align = 'center', fig.width = 7, fig.lp = 'fig:', fig.cap = "Graphical output of the manhattanplot() function", fig.pos = '!h'----
manhattanplot(wgscan.ihs.cgu,
              pval = TRUE,
              threshold = 4,
              main = "p-value of iHS (CGU cattle breed)")

## ----manhattanplotsub, fig.align = 'center', fig.width = 7, fig.lp = 'fig:', fig.cap = 'Graphical output of the manhattanplot() function', fig.pos = '!h'----
# re-define colors
palette(c("red", "green"))
manhattanplot(wgscan.ihs.cgu, 
              pval = TRUE,
              threshold = 4, 
              chr.name = c("1", "4", "5", "12"), 
              main = "iHS (CGU cattle breed)", 
              cr = cr.cgu,
              mrk = "F1205400",
              inset = 1E+7,
              resolution = c(200000, 0.05))
# set back to default colors
palette("default")

## ----rehh2qqman---------------------------------------------------------------
# extract data frame from result list
ihs <- wgscan.ihs.cgu$ihs
# create new data frame
wgscan.cgu.ihs.qqman <- data.frame(
    CHR = as.integer(factor(ihs$CHR, 
                            levels = unique(ihs$CHR))),
                               # chromosomes as integers
    BP = ihs$POSITION,         # base pairs
    P = 10**(-ihs$LOGPVALUE),  # transform back to p-values
    SNP = row.names(ihs)       # SNP names
    )

## ----qqman, eval = requireNamespace("qqman", quietly = TRUE), fig.align = 'center', fig.width = 7, fig.lp = 'fig:', fig.cap = 'Graphical output of the qqman::manhattan() function', fig.pos = '!h', message = FALSE----
library(qqman)
qqman::manhattan(wgscan.cgu.ihs.qqman,
                 col = c("red","green"),
                 chrlabs = unique(ihs$CHR),
                 suggestiveline = 4,
                 highlight = "F1205400",
                 annotatePval = 0.0001)

## ----plothaplo, results = "hide", fig.align = 'center', fig.width = 7, fig.height = 5, fig.lp = 'fig:', fig.cap = 'Graphical output of the plot.haplohh() function', fig.pos = '!h'----
hh_subset <- subset(haplohh_cgu_bta12, select.mrk = 350:550)
oldpar <- par(mar = c(3, 2, 2, 2) + 0.1)
plot(
  hh_subset,
  mrk = "F1205400",
  group_by_allele = TRUE,
  ignore.distance = TRUE,
  col = c(NA, "red"),
  linecol = c("lightblue", "lightpink"),
  mrk.col = "black",
  cex = 0.1,
  pos.lab.hap = "none",
  pos.lab.mrk = "none"
)
par(oldpar)

## -----------------------------------------------------------------------------
data(haplohh_cgu_bta12)
furcation <- calc_furcation(haplohh_cgu_bta12,
                            mrk = "F1205400")

## ----furc, fig.align = 'center', fig.width = 7, fig.height = 7.5, fig.lp = 'fig:', fig.cap = 'Graphical output of the plot.furcation() function', fig.pos = "!h"----
plot(furcation)

## ----furczoom, fig.align = 'center', fig.width = 7, fig.height = 7.5, fig.lp = 'fig:', fig.cap = 'Graphical output of the plot.furcation() function', fig.pos = "!h"----
plot(furcation,
     xlim = c(2.8E+7, 3E+7),
     lwd = 0.05,
     hap.names = hap.names(haplohh_cgu_bta12),
     cex.lab = 0.3)

## -----------------------------------------------------------------------------
newick <- as.newick(furcation,
                    allele = 0,
                    side = "left",
                    hap.names = hap.names(haplohh_cgu_bta12))

## ----newick, eval = requireNamespace("ape", quietly = TRUE), fig.align = 'center', fig.width = 6, fig.height = 6, fig.lp = 'fig:', fig.cap = 'Graphical output of the ape::plot.phylo() function', fig.pos = "!h"----
library(ape)
tree <- ape::read.tree(text = newick)
plot(tree, 
     cex = 0.5, 
     direction = "leftwards", 
     edge.color = "blue",
     underscore = TRUE,
     no.margin = TRUE)

## -----------------------------------------------------------------------------
haplen <- calc_haplen(furcation)

## -----------------------------------------------------------------------------
haplen$mrk.name

## -----------------------------------------------------------------------------
haplen$position

## -----------------------------------------------------------------------------
haplen$xlim

## -----------------------------------------------------------------------------
head(haplen$haplen)

## ----haplo, fig.align = 'center', fig.width = 7, fig.height = 7.5, fig.lp = 'fig:', fig.cap = 'Graphical output of the plot.haplen() function', fig.pos = "!h"----
plot(haplen)

## ----haplozoom, fig.align = 'center', fig.width = 7, fig.height = 7.5, fig.lp = 'fig:', fig.cap = 'Graphical output of the plot.haplen() function', fig.pos = "!h"----
plot(haplen,
     allele = 0,
     xlim = c(haplen$xlim[1], haplen$position),
     hap.names = hap.names(haplohh_cgu_bta12),
     cex.lab = 0.3,
     legend.xy.coords = "none")

## -----------------------------------------------------------------------------
# finding the index number of marker "F1205400"
mrk.nr = which(mrk.names(haplohh_cgu_bta12) == "F1205400")
# subset of all markers on the "left" of the focal one
hh_left = subset(haplohh_cgu_bta12, select.mrk = 1:mrk.nr)
# check for duplicated rows
which(duplicated(haplo(hh_left)))
# row 248 is identical to a previous row, but which one?
# get the other row by a search in opposite direction
which(duplicated(haplo(hh_left), fromLast = TRUE))
# extract the corresponding haplotype names
hap.names(hh_left)[c(236, 248)]

## -----------------------------------------------------------------------------
remove.example.files()

