## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(comment = ">", fig.height = 4.5, fig.width = 4.5, fig.show = "hold")

## ----library, message = FALSE-------------------------------------------------
library(rehh)

## -----------------------------------------------------------------------------
make.example.files()

## -----------------------------------------------------------------------------
cat(readLines("example1.hap"), sep = "\n")

## -----------------------------------------------------------------------------
cat(readLines("example1.map"), sep = "\n")

## -----------------------------------------------------------------------------
cat(readLines("example1.vcf"), sep = "\n")

## -----------------------------------------------------------------------------
hh <- data2haplohh(hap_file = "example1.hap",
                   map_file = "example1.map",
                   allele_coding = "01")

## -----------------------------------------------------------------------------
hh_map <- data2haplohh(hap_file = "example1.hap",
                       map_file = "example1.map",
                       allele_coding = "map",
                       verbose = FALSE)
identical(hh, hh_map)

hh_none <- data2haplohh(hap_file = "example1.hap",
                        map_file = "example1.map",
                        allele_coding = "none",
                        verbose = FALSE)
identical(hh, hh_none)

## ---- eval = requireNamespace("data.table", quietly = TRUE)-------------------
hh_vcf <- data2haplohh(hap_file = "example1.vcf",
                       vcf_reader = "data.table",
                       verbose = FALSE)
identical(hh, hh_vcf)

## -----------------------------------------------------------------------------
plot(hh)

## -----------------------------------------------------------------------------
res <- calc_ehh(hh, mrk = "rs6", include_nhaplo = TRUE)
res

## -----------------------------------------------------------------------------
haplo(hh)

## ----ehh, fig.align = 'center', fig.cap = "Graphical output of the plot.ehh() function", fig.pos = '!h', fig.lp = 'fig:'----
plot(res)

## -----------------------------------------------------------------------------
res <- calc_ehh(hh, 
                mrk = "rs6", 
                include_nhaplo = TRUE, 
                phased = FALSE)
res

## ----ehhunphased, fig.align = 'center', fig.cap = "Graphical output of the plot.ehh() function", fig.pos = '!h', fig.lp = 'fig:'----
plot(res)

## -----------------------------------------------------------------------------
res <- calc_ehhs(hh,
                 mrk = "rs6", 
                 include_nhaplo = TRUE,
                 discard_integration_at_border = FALSE)
res

## ----ehhs1, fig.align = 'center', fig.cap = "Graphical output of the plot.ehhs() function", fig.pos = '!h', fig.lp = 'fig:'----
plot(res)

## -----------------------------------------------------------------------------
res.scan <- scan_hh(hh, discard_integration_at_border = FALSE)
res.scan

## -----------------------------------------------------------------------------
ihs <- ihh2ihs(res.scan, freqbin = 1, verbose = FALSE)
ihs

## -----------------------------------------------------------------------------
cr <- calc_candidate_regions(ihs, threshold = 1.5, ignore_sign = TRUE, window_size = 1)
cr

## ----manhattan11, fig.align = 'center', fig.cap = "Graphical output of the manhattanplot() function", fig.pos = '!h', fig.lp = 'fig:'----
manhattanplot(ihs, threshold = c(-1.5,1.5), cr = cr, ylim = c(-2.5,2.5), pch = 20)

## ----furcation11, fig.align = 'center', fig.cap = "Graphical output of the plot.furcation() function", fig.pos = '!h', fig.lp = 'fig:'----
f <- calc_furcation(hh, mrk = "rs6")
# set equal plot margins on left and right side and save old ones
oldpar <- par(mar = (c(5, 3, 4, 3) + 0.1))
plot(f,
     # increase line width
     lwd = 1.5,
     # set habplotype identifiers as labels
     hap.names = hap.names(hh),
     # find a place for the legend ...
     legend.xy.coords = c(60000, 0.2)
)
# reset old margins
par(oldpar)

## ----newick1, , eval = requireNamespace("ape", quietly = TRUE), fig.align = 'center', fig.cap = "Graphical output of the plot.phylo() function of package ape", fig.pos = '!h', fig.lp = 'fig:'----
newick <- as.newick(f, 
                    allele = 0, 
                    side = "left", 
                    hap.names = hap.names(hh))
newick
library(ape)
tree <- ape::read.tree(text = newick)
plot(tree, 
     direction = "leftwards", 
     edge.color = "blue",
     underscore = TRUE,
     no.margin = TRUE)

## ----haplen11, fig.align = 'center', fig.cap = "Graphical output of the plot.haplen() function", fig.pos = '!h', fig.lp = 'fig:'----
h <- calc_haplen(f)
plot(h, 
     hap.names = hap.names(hh), 
     legend.xy.coords = c(70000,1.15))

## ----furcation12, fig.align = 'center', fig.cap = "Graphical output of the plot.haplen() function", fig.pos = '!h', fig.lp = 'fig:'----
f <- calc_furcation(hh, mrk = "rs6", phased = FALSE)
# set equal plot margins on left and right side and save old ones
oldpar <- par(mar = (c(5, 3, 4, 3) + 0.1))
plot(f,
     # increase line width
     lwd = 1.5,
     # set haplotype identifiers as labels
     hap.names = hap.names(hh),
     # no place for a legend inside the plot
     legend.xy.coords = "none")
# reset old margins
par(oldpar)

## ----haplen13, fig.align = 'center', fig.cap = "Graphical output of the plot.haplen() function", fig.pos = '!h', fig.lp = 'fig:'----
h <- calc_haplen(f)
plot(h, hap.names = hap.names(hh))

## -----------------------------------------------------------------------------
cat(readLines("example2.hap"),sep = "\n")

## -----------------------------------------------------------------------------
cat(readLines("example2.map"),sep = "\n")

## -----------------------------------------------------------------------------
cat(readLines("example2.vcf"), sep = "\n")

## -----------------------------------------------------------------------------
hh <- data2haplohh(hap_file = "example2.hap",
                   map_file = "example2.map",
                   allele_coding = "01")

## -----------------------------------------------------------------------------
hh <- data2haplohh(hap_file = "example2.hap",
                   map_file = "example2.map",
                   allele_coding = "01",
                   min_perc_geno.mrk = 50)

## -----------------------------------------------------------------------------
hh_map <- data2haplohh(hap_file = "example2.hap",
                       map_file = "example2.map",
                       allele_coding = "map",
                       min_perc_geno.mrk = 50,
                       verbose = FALSE)
identical(hh, hh_map)

hh_none <- data2haplohh(hap_file = "example2.hap",
                        map_file = "example2.map",
                        allele_coding = "none",
                        min_perc_geno.mrk = 50,
                        verbose = FALSE)
identical(hh, hh_none)

## ---- eval = requireNamespace("data.table", quietly = TRUE)-------------------
hh_vcf <- data2haplohh(hap_file = "example2.vcf",
                       min_perc_geno.mrk = 50,
                       vcf_reader = "data.table",
                       verbose = FALSE)
identical(hh, hh_vcf)

## -----------------------------------------------------------------------------
plot(hh)

## -----------------------------------------------------------------------------
res <- calc_ehh(hh, 
                mrk = "rs6", 
                include_nhaplo = TRUE,
                include_zero_values = TRUE,
                discard_integration_at_border = FALSE )
res

## ----ehh21, fig.align = 'center', fig.cap = "Graphical output of the plot.ehh() function", fig.pos = '!h', fig.lp = 'fig:'----
plot(res)

## -----------------------------------------------------------------------------
res <- calc_ehhs(hh,
                 mrk = "rs6", 
                 include_nhaplo = TRUE,
                 include_zero_values = TRUE,
                 discard_integration_at_border = FALSE)
res

## ----ehhs21, fig.align = 'center', fig.cap = "Graphical output of the plot.ehh() function", fig.pos = '!h', fig.lp = 'fig:'----
plot(res)

## -----------------------------------------------------------------------------
scan <- scan_hh(hh, discard_integration_at_border = FALSE)
scan

## -----------------------------------------------------------------------------
scan_unpol <- scan_hh(hh, 
                      discard_integration_at_border = FALSE,
                      polarized = FALSE)
scan_unpol

## -----------------------------------------------------------------------------
ihs <- ihh2ihs(scan, freqbin = 1)
ihs

## ----manhattan21, fig.align = 'center', fig.cap = "Graphical output of the manhattanplot() function", fig.pos = '!h', fig.lp = 'fig:'----
manhattanplot(ihs, threshold = c(-1.5, 1.5), ylim = c(-2.5,2.5), pch = 20)

## -----------------------------------------------------------------------------
hh_subset = subset(hh,
                   # exclude haplotypes with allele 2 at "rs6"
                   select.hap = haplo(hh)[ ,"rs6"] != 2, 
                   min_perc_geno.mrk = 50)
scan <- scan_hh(hh_subset, discard_integration_at_border = FALSE)
scan

## ----manhattan22, fig.align = 'center', fig.cap = "Graphical output of the manhattanplot() function", fig.pos = '!h', fig.lp = 'fig:'----
ihs <- ihh2ihs(scan, freqbin = 1, verbose = FALSE)
manhattanplot(ihs, threshold = c(-1.5, 1.5), ylim = c(-2.5,2.5), pch = 20)

## ----furcation21, fig.align = 'center', fig.cap = "Graphical output of the plot.furcation() function", fig.pos = '!h', fig.lp = 'fig:'----
f <- calc_furcation(hh, mrk = "rs6")
# set equal plot margins on left and right side and save old ones
oldpar <- par(mar = (c(5, 3, 4, 3) + 0.1))
plot(f,
     lwd = 1.5,
     hap.names = hap.names(hh),
     # no place for a legend inside the plot
     legend.xy.coords = "none")
par(oldpar)

## ----newick2, , eval = requireNamespace("ape", quietly = TRUE), fig.align = 'center', fig.cap = "Graphical output of the plot.phylo() function of package ape", fig.pos = '!h', fig.lp = 'fig:'----
newick <- as.newick(f, 
                    allele = 0, 
                    side = "left", 
                    hap.names = hap.names(hh))
tree <- ape::read.tree(text = newick)
plot(tree, 
     direction = "leftwards", 
     edge.color = "blue",
     underscore = TRUE,
     no.margin = TRUE)


## ----haplen21, fig.align = 'center', fig.cap = "Graphical output of the plot.haplen() function", fig.pos = '!h', fig.lp = 'fig:'----
h <- calc_haplen(f)
plot(h, hap.names = hap.names(hh), legend.xy.coords = "none")

## -----------------------------------------------------------------------------
remove.example.files()

