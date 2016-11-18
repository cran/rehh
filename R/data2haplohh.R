setClass(Class = "haplohh",
	representation(haplo = "matrix",position = "numeric",snp.name="character",chr.name = "character",nhap = "numeric",nsnp = "numeric")
)


is.haplohh <- function (x) 
{
    res <- (is(x,"haplohh") & validObject(x))
    return(res)
}


make.example.files <- function() {
     file.copy(system.file('bta12_hapguess_switch.out.zip',package='rehh.data'),'bta12_hapguess_switch.out.zip')
     unzip('bta12_hapguess_switch.out.zip')
     file.copy(system.file('map.inp.zip',package='rehh.data'),'map.inp.zip')
     unzip('map.inp.zip')
     file.copy(system.file('bta12_cgu.hap.zip',package='rehh.data'),'bta12_cgu.hap.zip')
     unzip('bta12_cgu.hap.zip')
     file.copy(system.file('bta12_cgu.thap.zip',package='rehh.data'),'bta12_cgu.thap.zip')
     unzip('bta12_cgu.thap.zip')
}

data2haplohh<-function(hap_file,map_file,min_maf=0,min_perc_geno.hap=100,min_perc_geno.snp=100,chr.name=NA,popsel=NA,recode.allele=FALSE,haplotype.in.columns=FALSE){
res<-new("haplohh")
if(min_perc_geno.hap<0 | min_perc_geno.hap>100){stop("The value for the argument min_perc_geno.hap should lie between 0 and 100")}
if(min_perc_geno.snp<0 | min_perc_geno.snp>100){stop("The value for the argument min_perc_geno.snp should lie between 0 and 100")}
#####Fichier map (verification)
map<-read.table(map_file,row.names=1,colClasses="character") #snp name, chromosome, position, allele ancestral, allele derive
 if(ncol(map)!=4){
   cat("Wrong format for map data file: ",map_file,"\n Should contain 5 columns (Snp name, Chromosome Number, Position, Ancestral Allele and Derived Allele) and no header\n")
  stop("Conversion stopped")
 }else{
   tmp_chr=unique(as.character(map[,1]))
   if(length(tmp_chr)!=1 & is.na(chr.name)){
    cat("More than one chromosome name in Map file:",map_file,"\n")
       repeat {
         cat('Which chromosome should be considered among:\n',tmp_chr,"\n")
         chr.name <- scan(file = '',n = 1,what = character(),quiet = TRUE)
         if (chr.name %in% tmp_chr) break
         }
    }
    
    if(is.na(chr.name)){chr.name=tmp_chr[1]}
     map=map[as.character(map[,1])==chr.name,]
     tmp_nom<-rownames(map)
     if(nrow(map)==0){
      cat("No SNPs mapping to chromosome ",chr.name," are found in the map file\n")
      tmp_chr=unique(as.character(map[,1]))
      cat("Here is the list of chromosome names in the map file:\n",tmp_chr,"\n")
      stop("Conversion stopped")
      }
   #premier checks OK
    res@chr.name<-as.character(chr.name)
    res@snp.name<-tmp_nom
    tmp_pos<-as.numeric(map[,2])

    if(sum(diff(tmp_pos)<0)>0){
      stop("SNP should be ordered on the map, check also that both haplotypes (column of haplo)\nand map (row of map) are ordered in the same way")}
    if(sum(diff(tmp_pos)==0)>0){
      warning("Some SNPs map to the same position")}
    res@position<-tmp_pos
  }

res@nsnp=nrow(map)
cat("Map file seems OK:",res@nsnp," SNPs declared for chromosome",res@chr.name,"\n")

###Fichier haplo
if(haplotype.in.columns){
  cat("Haplotype are in columns with no header\n")
  tmp.nhap=length(unlist(strsplit(readLines(hap_file,n=1),split="\t|\\s+"))) 
  tmp_haplos=matrix(scan(hap_file,what="character",quiet=TRUE),nrow=tmp.nhap)
#  tmp_haplos=t(as.matrix(read.table(hap_file,colClasses="numeric")))
#  tmp_haplos=readChar(hap_file,file.info(hap_file)$size,useBytes=T)
#  tmp_haplos=matrix(strsplit(gsub(tmp_haplos,pattern="\n|\t|\\s+",replacement=" "),split=" ",useBytes=T,fixed=T)[[1]],nrow=182)
  tmp.nsnp=ncol(tmp_haplos)
   if(tmp.nsnp!=res@nsnp){
    cat("The number of snp in the haplotypes",tmp.nsnp," is not equal\nto the number of snps declared in the map file",res@nsnp,"\n")
    stop("Conversion stopped")
   }
}else{
 out_fphase<-scan(hap_file,what="character",sep="\n",quiet=TRUE,nlines=15)
 test_fphase_1=grep("fastPHASE",out_fphase)
 test_fphase_2=grep("BEGIN COMMAND_LINE",out_fphase)

 if(length(test_fphase_1)>0 & length(test_fphase_2)>0){ #fichier fastphase
  cat("Looks like a FastPHASE haplotype file\n")
  out_fphase<-scan(hap_file,what="character",sep="\n",quiet=TRUE)
  deb_geno=grep("BEGIN GENOTYPES",out_fphase)[1] + 1
  fin_geno=grep("END GENOTYPES",out_fphase)[1] - 1
  out_fphase<-out_fphase[deb_geno:fin_geno]

 test_poplabel=grep("subpop. label:",out_fphase)
 if(length(test_poplabel)>1){
   nom_hap_cplet=out_fphase[test_poplabel]
   nhap_tot=length(nom_hap_cplet)
   pop_label=numeric(nhap_tot)
   tmp_poplab=(strsplit(nom_hap_cplet,split="subpop. label:"))
   for(i in 1:nhap_tot){
     pop_label[i]=as.numeric(unlist(strsplit(tmp_poplab[[i]][2],split="\\(internally"))[1])
    }
    liste_pop=unique(pop_label)
    cat("Haplotypes originate from ",length(liste_pop)," different populations in the fastPhase output file\n")
   
    if(!(popsel %in% pop_label)){
     cat("Chosen pop. is not in the list of pop. number:\n",liste_pop,"\n")
     popsel=NA
    }
    if(is.na(popsel)){ #on demande interactivement la pop a selectionner
       repeat {
         cat('Which population should be considered among:',liste_pop,"\n")
         popsel <- scan(file = '',n = 1,what = integer(0),quiet = TRUE)
         if (popsel %in% pop_label) break
         }
    }
   hapsel=(which(pop_label==popsel) -1)*3 + 1
   hapsel=sort(as.numeric(cbind(hapsel,hapsel+1,hapsel+2)))
   out_fphase=out_fphase[hapsel]
   }else{
    cat("Haplotypes seems to originate from only one population\n")
    if(!(is.na(popsel))){ #bizarre car pas de sous population observe
         cat('No popsel are thus considered:',liste_pop,"\n")
    }
 }

 nlignes=length(out_fphase) ; nind=nlignes/3 ; nhap=2*nind
 tmp_haplos=matrix(,nhap,res@nsnp)
   for(i in 1:nind){
    id_line=3*(i-1) + 1 ; hap1_line=id_line+1 ; hap2_line=id_line+2
    hap1_index=2*(i-1) +1 ; hap2_index=hap1_index + 1
  
    hap1=unlist(strsplit(out_fphase[hap1_line],split=" "))
    if(length(hap1)!=res@nsnp){
      stop("Number of snp in haplotypes differ from the number declared in the map file")
      }else{tmp_haplos[hap1_index,]=hap1}
  
    hap2=unlist(strsplit(out_fphase[hap2_line],split=" "))
    if(length(hap2)!=res@nsnp){
      stop("Number of snp in haplotypes differ from the number declared in the map file")
      }else{tmp_haplos[hap2_index,]=hap2}
   }

 }else{ #fichier au format standard
   cat("Standard rehh input file assumed\n")
   tmp.ncol=length(unlist(strsplit(readLines(hap_file,n=1),split="\t|\\s+"))) #retrocompatibilite: ca marchait avec la tabulation avant
   tmp.nsnp=tmp.ncol-1 #retrocompatibilite: les haplos contiennent un ID
   if(tmp.nsnp!=res@nsnp){
    cat("The number of snp in the haplotypes",tmp.nsnp," is not equal\nto the number of snps declared in the map file",res@nsnp,"\n")
    stop("Conversion stopped")
   }
   tmp_haplos=matrix(scan(hap_file,what="character",quiet=TRUE),ncol=tmp.ncol,byrow=TRUE)[,-1]
 }
 }
 res@nhap=nrow(tmp_haplos)
##recodage des alleles si necessaire
 if(recode.allele){
  cat("Alleles are being recoded according to map file as:\n\t0 (missing data), 1 (ancestral allele) or 2 (derived allele)\n")
  anc_allele=matrix(rep(map[,3],res@nhap),res@nhap,res@nsnp,byrow=TRUE)
  all_anc=tmp_haplos==anc_allele
  rm(anc_allele) 
  der_allele=matrix(rep(map[,4],res@nhap),res@nhap,res@nsnp,byrow=TRUE)
  all_der=tmp_haplos==der_allele
  rm(der_allele)
  rm(tmp_haplos) 
  res@haplo=matrix(0,res@nhap,res@nsnp)
  res@haplo[all_anc]=1 ; res@haplo[all_der]=2 
  }else{
   if(sum(!(tmp_haplos=="0" | tmp_haplos=="1" | tmp_haplos=="2" ))>0){
    cat("Alleles are not coded in the appropriate format:\n\t0 (missing data), 1 (ancestral allele) or 2 (derived allele)\n")
    cat("Check your data or use recode.allele=TRUE option to recode according to the map information\n")
    stop("Conversion stopped")
   }
  res@haplo=matrix(as.numeric(tmp_haplos),res@nhap,res@nsnp)
  rm(tmp_haplos) 
 }
#selection des haplos d'apres les donnees manquantes
# if( !(is.na(min_perc_geno.hap)) ){
   cat("Discard Haplotype with less than ",min_perc_geno.hap,"% of genotyped SNPs\n")
   hap_sel=(100*rowSums(res@haplo!=0)/res@nsnp)>=min_perc_geno.hap
   if(sum(hap_sel)==res@nhap){
    cat("No haplotype discarded\n")
   }else{
    cat(res@nhap-sum(hap_sel)," Haplos discarded\n")
    res@haplo=res@haplo[hap_sel,]
    res@nhap=sum(hap_sel)
    cat(res@nhap," Haplos remaining\n")
   }
   if(res@nhap==0){stop("No haplotype left after filtering of missing data: check allele coding (e.g., correspondance of allele coding between the map_inp and haplotype input file) data or relax the value of min_perc_geno.hap to allow for more missing data")}
# }
#selection des snps d'apres les donnees manquantes
# if( !(is.na(min_perc_geno.snp)) ){
   cat("Discard SNPs genotyped on less than ",min_perc_geno.snp,"% of haplotypes\n")
   snp_sel=(100*colSums(res@haplo!=0)/res@nhap)>=min_perc_geno.snp
   if(sum(snp_sel)==res@nsnp){
    cat("No SNP discarded\n")
   }else{
    cat(res@nsnp-sum(snp_sel)," SNPs discarded\n")
    res@haplo=res@haplo[,snp_sel]
    res@nsnp=sum(snp_sel)
    res@position=res@position[snp_sel]
    res@snp.name=res@snp.name[snp_sel]
    cat(res@nsnp," SNPs remaining\n")
   }
   if(res@nsnp==0){stop("No SNP left after filtering of missing data: check allele coding (e.g., correspondance  ofallele coding between the map_inp and haplotype input file) data or relax the value of min_perc_geno.snp to allow for more missing data")}
   # }

#selection des snps sur MAF
 if(min_maf>0){
   cat("Discard SNPs with MAF below ",min_maf,"\n")
   tmp_n1=colSums(res@haplo==1) ; tmp_n=colSums(res@haplo!=0)
   tmp_maf=0.5-abs(0.5-tmp_n1/tmp_n)
   snp_sel=tmp_maf>min_maf
   if(sum(snp_sel)==res@nsnp){
    cat("No SNP discarded\n")
   }else{
    cat(res@nsnp-sum(snp_sel)," SNPs discarded\n")
    res@haplo=res@haplo[,snp_sel]
    res@nsnp=sum(snp_sel)
    res@position=res@position[snp_sel]
    res@snp.name=res@snp.name[snp_sel]
    cat(res@nsnp," SNPs remaining\n")
   }
   rm(tmp_n1) ; rm(tmp_n) ; rm(tmp_maf)
 }

    cat("Data consists of",res@nhap,"haplotypes and",res@nsnp,"SNPs\n")
    return(res)
}
