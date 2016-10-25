calc_ehhs <- function(haplohh,mrk,limhaplo = 2,limehhs = 0.05,maxgap = NA,plotehhs = TRUE,lty = 1,lwd = 1.5,col = c("blue","red"),xlab = "Position",ylab = expression(Site~specific~italic(EHH)~(italic(EHHS))),cex.lab = 1.25,main = NA,cex.main = 1.5) {

	if (!(is.haplohh(haplohh))) {stop("The data are not formatted as a valid haplohh object... (see the data2haplohh() function)")} 
	if (mrk < 1 | mrk > haplohh@nsnp) {stop(paste("The focal SNP index must lie between",1,"and",haplohh@nsnp))}
	if (limhaplo < 2) {stop("limhaplo must be larger than 1")}
	if (limehhs < 0 | limehhs > 1) {stop("limehhs must lie between 0 and 1")}
	if (is.na(maxgap)) {maxgap = (max(haplohh@position) + 1)}

  nhaplo_eval <- ehhs_tang <- ehhs_sabeti <- rep(0,haplohh@nsnp)
  ies_tang <- ies_sabeti <- 0
  res.ehhs <- .C("CALL_EHHS",
  				Rdata = as.integer(haplohh@haplo),
  				focal_SNP = as.integer(mrk),
  				number_SNPs  = as.integer(haplohh@nsnp),
  				number_chromosomes = as.integer(haplohh@nhap),
  				number_haplotypes = as.integer(nhaplo_eval),
  				min_number_haplotypes = as.integer(limhaplo),
  				min_EHHS = as.double(limehhs),
  				max_gap = as.double(maxgap),
  				map = as.double(haplohh@position),
  				EHHS_TANG = as.double(ehhs_tang),
  				IES_TANG = as.double(ies_tang),
  				EHHS_SABETI = as.double(ehhs_sabeti),
  				IES_SABETI = as.double(ies_sabeti)
  				)
  				
	nhaplo_eval = res.ehhs$number_haplotypes 
	ehhs_tang = res.ehhs$EHHS_TANG
	ies_tang = res.ehhs$IES_TANG
	ehhs_sabeti = res.ehhs$EHHS_SABETI
	ies_sabeti = res.ehhs$IES_SABETI
	names(ehhs_tang) = names(ehhs_sabeti) = names(nhaplo_eval) = haplohh@snp.name
	
	if (plotehhs) {
		sel_reg <- (nhaplo_eval > 0)
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
			dev.new()
			plot.new()
			par(mar = c(5,5,4,2) + 0.1)
			plot(haplohh@position[sel_reg] / scale,ehhs_sabeti[sel_reg],col = col[1],type = "l",lty = lty,lwd = lwd,main = main,bty = "n",xlab = paste(xlab,unit),ylab = ylab,cex.lab = cex.lab, cex.lab = cex.main)
			lines(haplohh@position[sel_reg] / scale,ehhs_tang[sel_reg],col = col[2],lty = lty,lwd = lwd)
			abline(v = haplohh@position[mrk] / scale,lty = 2)
			if (haplohh@position[mrk] > sum(range(haplohh@position[sel_reg])) / 2) {
				legend("topleft",c("Sabeti et al. (2007)","Tang et al. (2007)"),col = col,bty = "n",lty = lty,lwd = lwd)
				
			}
			else {
				legend("topright",c("Sabeti et al. (2007)","Tang et al. (2007)"),col = col,bty = "n",lty = lty,lwd = lwd)
			}
		}
	}

	return(list(nhaplo_eval = nhaplo_eval,EHHS_Tang_et_al_2007 = ehhs_tang,IES_Tang_et_al_2007 = ies_tang,EHHS_Sabeti_et_al_2007 = ehhs_sabeti,IES_Sabeti_et_al_2007 = ies_sabeti))
}

