scan_hh <- function(haplohh,limhaplo = 2,limehh = 0.05,limehhs = 0.05,maxgap = NA,threads = 1) {
	
	if (!(is.haplohh(haplohh))) {stop("The data are not formatted as a valid haplohh object... (see the data2haplohh() function)")} 
	if (limhaplo < 2) {stop("limhaplo must be larger than 1")}
	if (limehh < 0 | limehh > 1) {stop("limehh must lie between 0 and 1")}
	if (limehhs < 0 | limehhs > 1) {stop("limehhs must lie between 0 and 1")}
	if (is.na(maxgap)) {maxgap = (max(haplohh@position) + 1)}

	ihh <- matrix(0.0,nrow = haplohh@nsnp,ncol = 2)
	ies_tang <- ies_sabeti <- vector(mode = "numeric",length = haplohh@nsnp)
	res_scan <- .C("CALL_SCAN_HH",
				Rdata = as.integer(haplohh@haplo),
				number_SNPs  = as.integer(haplohh@nsnp),
				number_chromosomes = as.integer(haplohh@nhap),
				min_number_haplotypes = as.integer(limhaplo),
				min_EHH = as.double(limehh),
				min_EHHS = as.double(limehhs),
  				max_gap = as.double(maxgap),
				map = as.double(haplohh@position),
				IHH = as.double(ihh),
				IES_TANG = as.double(ies_tang),
				IES_SABETI = as.double(ies_sabeti),
				number_threads = as.integer(threads)
				)
	tmp_n1 = colSums(haplohh@haplo == 1)
	freq_A = tmp_n1 / (tmp_n1 + colSums(haplohh@haplo == 2))
	CHR = rep(haplohh@chr.name,haplohh@nsnp)
	POSITION = haplohh@position
	tmp.ihh = matrix(res_scan$IHH,haplohh@nsnp,2)
	tmp.ihh[tmp.ihh == -1]=NA
	iHH_A = tmp.ihh[,1]
	iHH_D = tmp.ihh[,2]
	iES_Tang_et_al_2007 = replace(res_scan$IES_TANG,which(res_scan$IES_TANG == -1),NA)
	iES_Sabeti_et_al_2007 = replace(res_scan$IES_SABETI,which(res_scan$IES_SABETI == -1),NA)
	return(data.frame(CHR,POSITION,freq_A,iHH_A,iHH_D,iES_Tang_et_al_2007,iES_Sabeti_et_al_2007,row.names=haplohh@snp.name,stringsAsFactors=FALSE))
}
