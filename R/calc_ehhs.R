calc_ehhs <- function(haplohh,mrk,limhaplo = 2,limehhs = 0.05,maxgap = NA,discard_integration_at_border = TRUE,plotehhs = TRUE,lty = 1,lwd = 1.5,col = c("blue","red"),xlab = "Position",ylab = expression(Site~specific~italic(EHH)~(italic(EHHS))),cex.lab = 1.25,main = NA,cex.main = 1.5) {

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
				discard_integration_at_border = as.integer(discard_integration_at_border),
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
	ies_tang = replace(res.ehhs$IES_TANG,which(res.ehhs$IES_TANG == -1),NA)
	ehhs_sabeti = res.ehhs$EHHS_SABETI
	ies_sabeti = replace(res.ehhs$IES_SABETI,which(res.ehhs$IES_SABETI == -1),NA)
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
		    if (!is.null(names(dev.list())) && ((names(dev.cur()) == "windows") | (names(dev.cur()) == "X11") | (names(dev.cur()) == "quartz"))) {
   		    		dev.new()
    		}
			par(mar = c(5,5,4,2) + 0.1)
			plot(haplohh@position[sel_reg] / scale,ehhs_sabeti[sel_reg],col = col[1],type = "l",lty = lty,lwd = lwd,main = main,bty = "n",xlab = paste(xlab,unit),ylab = ylab,cex.lab = cex.lab, cex.lab = cex.main, ylim=c(0,1))
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

