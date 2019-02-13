calc_ehh <- function(haplohh,mrk,limhaplo = 2,limehh = 0.05,scalegap = NA,maxgap = NA,discard_integration_at_border = TRUE,plotehh = TRUE,lty = 1,lwd = 1.5,col = c("blue","red"),xlab = "Position",ylab = expression(Extended~haplotype~homozygosity~(italic(EHH))),cex.lab = 1.25,main = NA,cex.main = 1.5) {

	if (!(is.haplohh(haplohh))) {stop("The data are not formatted as a valid haplohh object... (see the data2haplohh() function)")}
  if(is.numeric(mrk)){
    mrk=as.integer(mrk)
    if(mrk<1){
      stop(paste0("No marker numbers smaller than 1 allowed"))
    }
    if(mrk>haplohh@nsnp){
      stop(paste0("The marker number ",mrk," is bigger than the number of SNPs in the data set (",haplohh@nsnp,")"))
    }
  }else{
    mrk = as.character(mrk)
    if (!(mrk %in% haplohh@snp.name)) {stop(paste0("A marker with name '",mrk,"' is not contained in the data set"))}
    mrk = which(haplohh@snp.name == mrk)
  }
	if (limhaplo < 2) {stop("limhaplo must be larger than 1")}
	if (limehh < 0 | limehh > 1) {stop("limehh must lie between 0 and 1")}
	if (is.na(maxgap)) {maxgap = (max(haplohh@position) + 1)}
  if (is.na(scalegap)){
    scalegap = (max(haplohh@position) + 1)
  } else{
    if (scalegap > maxgap) {
      stop("scalegap has to be smaller than maxgap in order to have an effect")
    }
  }
 	ehh <- nhaplo_eval <- matrix(0,nrow = haplohh@nsnp,ncol = 2)
	ihh <- rep(0,2)
	res.ehh <- .C("CALL_EHH",
  				Rdata = as.integer(haplohh@haplo),
  				focal_SNP = as.integer(mrk),
  				number_SNPs  = as.integer(haplohh@nsnp),
  				number_chromosomes = as.integer(haplohh@nhap),
  				number_haplotypes = as.integer(nhaplo_eval),
  				min_number_haplotypes = as.integer(limhaplo),
		  		discard_integration_at_border = as.integer(discard_integration_at_border),
  				min_EHH = as.double(limehh),
  				scale_gap = as.double(scalegap),
	  			max_gap = as.double(maxgap),
  				map = as.double(haplohh@position),
  				EHH = as.double(ehh),
  				IHH = as.double(ihh)
  				)
	nhaplo_eval = matrix(res.ehh$number_haplotypes,2,haplohh@nsnp,byrow = T)
	ehh = matrix(res.ehh$EHH,2,haplohh@nsnp,byrow = T)
	rownames(ehh) = rownames(nhaplo_eval) = c("Ancestral allele","Derived allele")
	colnames(ehh) = colnames(nhaplo_eval) = haplohh@snp.name
	ihh = replace(res.ehh$IHH,which(res.ehh$IHH == -1),NA)
	names(ihh) = c("Ancestral allele","Derived allele")
	if (plotehh) {
		sel_reg <- (colSums(nhaplo_eval) > 0)
		if (sum(sel_reg) > 0) {
			if (max(haplohh@position[sel_reg]) < 1e3) {
				scale <- 1
				unit = "(bp)"
			}
			else if (max(haplohh@position[sel_reg]) < 1e6) {
				scale <- 1e3
				unit = "(kb)"
			}
			else if (max(haplohh@position[sel_reg]) < 1e9) {
				scale <- 1e6
				unit = "(Mb)"
			}
			else {
				scale <- 1e9
				unit = "(Gb)"
			}
		    if (!is.null(names(dev.list())) && ((names(dev.cur()) == "windows") | (names(dev.cur()) == "X11") | (names(dev.cur()) == "quartz"))) {
   		    		dev.new()
    		}
			par(mar = c(5,5,4,2) + 0.1)
			matplot(haplohh@position[sel_reg] / scale,t(ehh[,sel_reg]),col = col,type = "l",lty = lty,lwd = lwd,main = main,bty = "n",xlab = paste(xlab,unit),ylab = ylab,cex.lab = cex.lab, cex.main = cex.main)
			abline(v = haplohh@position[mrk] / scale,lty = 2)
			if (haplohh@position[mrk] > sum(range(haplohh@position[sel_reg])) / 2) {
				legend("topleft",c("Ancestral allele","Derived allele"),col = col,bty = "n",lty = lty,lwd = lwd)
			}
			else {
				legend("topright",c("Ancestral allele","Derived allele"),col = col,bty = "n",lty = lty,lwd = lwd)
			}
		}
	}
	return(list(ehh = ehh,nhaplo_eval = nhaplo_eval,freq_all1 = nhaplo_eval[1,mrk] / sum(nhaplo_eval[,mrk]),ihh = ihh))
}
