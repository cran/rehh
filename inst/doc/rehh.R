## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(comment=">")

## ----message=FALSE-------------------------------------------------------
library(rehh)

## ----results='hide'------------------------------------------------------
make.example.files()

## ------------------------------------------------------------------------
head(read.table("map.inp"))

## ------------------------------------------------------------------------
?data2haplohh

## ------------------------------------------------------------------------
hap<-data2haplohh(hap_file="bta12_cgu.hap",map_file="map.inp",
                  recode.allele=TRUE,chr.name=12)

## ----eval=FALSE----------------------------------------------------------
#  hap<-data2haplohh(hap_file="bta12_cgu.hap",map_file="map.inp",
#                    recode.allele=TRUE)

## ----echo=FALSE----------------------------------------------------------
cat("More than one chromosome name in Map file: map.inp\nWhich chromosome should be considered among:\n1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29\n1:")

## ----eval=FALSE----------------------------------------------------------
#  12

## ----echo=FALSE----------------------------------------------------------
hap<-data2haplohh(hap_file="bta12_cgu.hap",map_file="map.inp",recode.allele=TRUE,chr.name=12)

## ----error=TRUE----------------------------------------------------------
hap<-data2haplohh(hap_file="bta12_cgu.hap",map_file="map.inp",
                  recode.allele=TRUE,chr.name=18)

## ------------------------------------------------------------------------
hap<-data2haplohh(hap_file="bta12_cgu.thap",map_file="map.inp",haplotype.in.columns=TRUE,
                  recode.allele=TRUE,chr.name=12)

## ------------------------------------------------------------------------
hap<-data2haplohh(hap_file="bta12_hapguess_switch.out",map_file="map.inp",
                  recode.allele=TRUE,popsel=7,chr.name=12)

## ----eval=FALSE----------------------------------------------------------
#  hap<-data2haplohh(hap_file="bta12_hapguess_switch.out",map_file="map.inp",
#                    recode.allele=TRUE,chr.name=12)

## ----echo=FALSE----------------------------------------------------------
cat("Map file seems OK: 1424  SNPs declared for chromosome 12\nLooks like a FastPHASE haplotype file\nHaplotypes originate from  8  different populations in the fastPhase output file\nChosen pop. is not in the list of pop. number: 1 2 3 4 5 6 7 8\nWhich population should be considered among: 1 2 3 4 5 6 7 8\n1: 
")

## ----eval=FALSE----------------------------------------------------------
#  7

## ----echo=FALSE----------------------------------------------------------
hap<-data2haplohh(hap_file="bta12_hapguess_switch.out",map_file="map.inp",recode.allele=TRUE,popsel=7,chr.name=12)

## ----eval=FALSE----------------------------------------------------------
#  ?calc_ehh

## ---- ehhplot,fig.align='center',out.height="7cm",fig.cap='Graphical output of the \\texttt{calc\\_ehh()} function',fig.pos='!h',fig.lp='fig:',fig.keep='all'----
#example haplohh object (280 haplotypes, 1424 SNPs) see ?haplohh_cgu_bta12 for details
data(haplohh_cgu_bta12)
#computing EHH statistics for the focal SNP with name "F1205400" 
#which displays a strong signal of selection
res.ehh<-calc_ehh(haplohh_cgu_bta12,mrk="F1205400") 

## ------------------------------------------------------------------------
res.ehh$ehh[1:2,454:458]

## ------------------------------------------------------------------------
res.ehh$nhaplo_eval[1:2,454:458]

## ------------------------------------------------------------------------
res.ehh$freq_all1

## ------------------------------------------------------------------------
res.ehh$ihh

## ---- eval=FALSE---------------------------------------------------------
#  ?calc_ehhs

## ----ehhsplot,fig.align='center',out.height="7cm",fig.cap='Graphical output of the \\texttt{calc\\_ehhs()} function',fig.pos='!h',fig.lp='fig:'----
data(haplohh_cgu_bta12)
res.ehhs<-calc_ehhs(haplohh_cgu_bta12,mrk="F1205400")

## ------------------------------------------------------------------------
res.ehhs$EHHS_Sabeti_et_al_2007[453:459] 

## ------------------------------------------------------------------------
res.ehhs$EHHS_Tang_et_al_2007[453:459] 

## ------------------------------------------------------------------------
res.ehhs$nhaplo_eval[453:459] 

## ------------------------------------------------------------------------
res.ehhs$IES_Tang_et_al_2007

## ------------------------------------------------------------------------
res.ehhs$IES_Sabeti_et_al_2007

## ----eval=TRUE-----------------------------------------------------------
data(haplohh_cgu_bta12)
res.scan<-scan_hh(haplohh_cgu_bta12)

## ------------------------------------------------------------------------
dim(res.scan)

## ------------------------------------------------------------------------
res.scan[453:459,]

## ------------------------------------------------------------------------
system.time(res.scan<-scan_hh(haplohh_cgu_bta12))

## ------------------------------------------------------------------------
foo<-function(haplo){
  res.ihh=res.ies=matrix(0,haplo@nsnp,2)
  for(i in 1:haplo@nsnp){
    res.ihh[i,]=calc_ehh(haplo,mrk=haplo@snp.name[i],plotehh=FALSE)$ihh
    tmp=calc_ehhs(haplo,mrk=haplo@snp.name[i],plotehhs=FALSE)
    res.ies[i,1]=tmp$IES_Tang_et_al_2007
    res.ies[i,2]=tmp$IES_Sabeti_et_al_2007  
  }
  list(res.ies=res.ies,res.ihh=res.ihh)
}
system.time(res.scan2<-foo(haplohh_cgu_bta12))

## ------------------------------------------------------------------------
identical(res.scan2$res.ihh[,1],res.scan[,4])
identical(res.scan2$res.ihh[,2],res.scan[,5])
identical(res.scan2$res.ies[,1],res.scan[,6])
identical(res.scan2$res.ies[,2],res.scan[,7])

## ----eval=FALSE----------------------------------------------------------
#  for(i in 1:29){
#    hap_file=paste("hap_chr_",i,".pop1",sep="")
#    data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
#    res<-scan_hh(data)
#    if(i==1){wg.res<-res}else{wg.res<-rbind(wg.res,res)}
#  }
#  wg.ihs<-ihh2ihs(wg.res)

## ------------------------------------------------------------------------
data(wgscan.cgu)
## results from a genome scan (44,057 SNPs) see ?wgscan.eut and ?wgscan.cgu for details
ihs.cgu<-ihh2ihs(wgscan.cgu)

## ------------------------------------------------------------------------
head(ihs.cgu$iHS)

## ------------------------------------------------------------------------
head(ihs.cgu$frequency.class)

## ----ihsplot,fig.align='center',fig.width=16,fig.height=12,fig.lp='fig:',fig.cap='Graphical output of the \\texttt{ihsplot()} function',fig.pos='!h'----
layout(matrix(1:2,2,1))
ihsplot(ihs.cgu,plot.pval=TRUE,ylim.scan=2,main="iHS (CGU cattle breed)")

## ----eval=FALSE----------------------------------------------------------
#  for(i in 1:29){
#    hap_file=paste("hap_chr_",i,".pop1",sep="")
#    data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
#    res<-scan_hh(data)
#    if(i==1){wg.res.pop1<-res}else{wg.res.pop1<-rbind(wg.res.pop1,res)}
#    hap_file=paste("hap_chr_",i,".pop2",sep="")
#    data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
#    res<-scan_hh(data)
#    if(i==1){wg.res.pop2<-res}else{wg.res.pop2<-rbind(wg.res.pop2,res)}
#  }
#  wg.rsb<-ies2rsb(wg.res.pop1,wg.res.pop2)

## ------------------------------------------------------------------------
data(wgscan.cgu) ; data(wgscan.eut)
## results from a genome scan (44,057 SNPs) see ?wgscan.eut and ?wgscan.cgu for details
cguVSeut.rsb<-ies2rsb(wgscan.cgu,wgscan.eut,"CGU","EUT")

## ------------------------------------------------------------------------
head(cguVSeut.rsb)

## ----rsbplot,fig.align='center',fig.width=16,fig.height=12,fig.lp='fig:',fig.cap='Graphical output of the \\texttt{rsbplot()} function',fig.pos='!h'----
layout(matrix(1:2,2,1))
rsbplot(cguVSeut.rsb,plot.pval=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  for(i in 1:29){
#    hap_file=paste("hap_chr_",i,".pop1",sep="")
#    data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
#    res<-scan_hh(data)
#    if(i==1){wg.res.pop1<-res}else{wg.res.pop1<-rbind(wg.res.pop1,res)}
#    hap_file=paste("hap_chr_",i,".pop2",sep="")
#    data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
#    res<-scan_hh(data)
#    if(i==1){wg.res.pop2<-res}else{wg.res.pop2<-rbind(wg.res.pop2,res)}
#  }
#  wg.xpehh<-ies2xpehh(wg.res.pop1,wg.res.pop2)

## ------------------------------------------------------------------------
data(wgscan.cgu) ; data(wgscan.eut)
## results from a genome scan (44,057 SNPs) see ?wgscan.eut and ?wgscan.cgu for details
cguVSeut.xpehh<-ies2xpehh(wgscan.cgu,wgscan.eut,"CGU","EUT")

## ------------------------------------------------------------------------
head(cguVSeut.xpehh)

## ----xpehhplot,fig.align='center',fig.width=16,fig.height=12,fig.lp='fig:',fig.cap='Graphical output of the \\texttt{xpehhplot()} function',fig.pos="!h"----
layout(matrix(1:2,2,1))
xpehhplot(cguVSeut.xpehh,plot.pval=TRUE)

## ----comp,echo=TRUE,fig.align='center',fig.width=6,fig.height=6,out.height="8cm",fig.lp='fig:',fig.cap='Comparison of Rsb and XP-EHH values across the CGU and EUT populations',fig.pos="!h"----
plot(cguVSeut.rsb[,3],cguVSeut.xpehh[,3],xlab="Rsb",ylab="XP-EHH",pch=".",
     xlim=c(-7.5,7.5),ylim=c(-7.5,7.5))
points(cguVSeut.rsb["F1205400",3],cguVSeut.xpehh["F1205400",3],col="red")
abline(a=0,b=1,lty=2)

## ----distribplot,fig.align='center',fig.width=6,fig.height=12,out.height="14cm",fig.lp='fig:',fig.cap='Graphical output of the function \\texttt{distribplot}',fig.pos="!h"----
layout(matrix(1:2,2,1))
distribplot(ihs.cgu$iHS[,3],xlab="iHS")

## ----bifdia,eval=TRUE,fig.align='center',out.height="12cm",fig.width=16,fig.height=12,fig.lp='fig:',fig.cap='Graphical output of the function \\texttt{bifurcation.diagram()}',fig.pos="!h"----
data(haplohh_cgu_bta12)
layout(matrix(1:2,2,1))
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc="F1205400",all_foc=1,nmrk_l=20,nmrk_r=20,
                    main="Bifurcation diagram (RXFP2 SNP on BTA12): Ancestral Allele")
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc="F1205400",all_foc=2,nmrk_l=20,nmrk_r=20,
                    main="Bifurcation diagram (RXFP2 SNP on BTA12): Derived Allele")

