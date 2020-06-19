#include "calc_ehh.h"
#include "haplotypes.h"
#include "homozygosity.h"

/**
 * Calculates ehh values in increasing distance (either in positive or negative direction) form the focal marker.
 */
void extend_ehh(const int* const data, const int nbr_chr, const int foc_mrk, const int end_mrk, const int lim_haplo, const int lim_homo_haplo,
		const double lim_ehh, const bool phased, int* const hap, int* const nbr_hap, int* const nbr_chr_with_hap,
		int* const tot_nbr_chr_in_hap, double* const ehh) {

	int increment = foc_mrk <= end_mrk ? 1 : -1;

	for (int mrk = foc_mrk + increment; mrk != end_mrk + increment; mrk += increment) { //walk along the chromosome, away from the focal SNP
		if (update_hap(data, nbr_chr, mrk, hap, nbr_hap, nbr_chr_with_hap)) { //if there is a change in haplotypes (i.e. at least one furcation)
			for (int i = 0; i < *nbr_hap; i++) {
				tot_nbr_chr_in_hap[mrk] += nbr_chr_with_hap[i]; //total number of chromosomes; changes only by missing data
			}
			if (tot_nbr_chr_in_hap[mrk] < lim_haplo) {
				break; //stop, if too few chromosomes left due to missing values
			}
			
			//if unphased, sequences from different individuals are different by definition
	    int nbr_homo_haplo = phased ? tot_nbr_chr_in_hap[mrk] - *nbr_hap + 1 : (tot_nbr_chr_in_hap[mrk] - *nbr_hap) * 2;
			if(nbr_homo_haplo < lim_homo_haplo){
		    break; //stop, if too few homozygous chromosomes left
		  }
			
			ehh[mrk] = homozygosity(tot_nbr_chr_in_hap[mrk], *nbr_hap, nbr_chr_with_hap, phased);

			if (ehh[mrk] <= lim_ehh) { //if homozygosity is small (nearly all haplotypes are different)
				ehh[mrk] = 0.0; // reset to zero
				break; //stop calculation
			}
		} else { //if no change in the number of haplotypes, then re-use previous values
			tot_nbr_chr_in_hap[mrk] = tot_nbr_chr_in_hap[mrk - increment];
			ehh[mrk] = ehh[mrk - increment];
		}
	}
}

/**
 * Calculates ehh values in increasing distance to the focal marker
 */
void calc_ehh(const int* const data, const int nbr_chr, const int nbr_mrk, const int foc_mrk, const int foc_all,
		const int lim_haplo, const int lim_homo_haplo, const double lim_ehh, const bool phased, int* const tot_nbr_chr_in_hap, double* const ehh) {

	int nbr_hap;                                                   //number of distinct haplotypes

	int *hap = (int*) malloc(nbr_chr * sizeof(int));               //vector index to chromosomes, ordered by haplotype
	int *nbr_chr_with_hap = (int *) malloc(nbr_chr * sizeof(int)); //for each haplotype gives the number of chromosomes sharing it

	//initialization of return values
	for (int mrk = 0; mrk < nbr_mrk; mrk++) {
		ehh[mrk] = 0.0;
		tot_nbr_chr_in_hap[mrk] = 0;
	}

	init_allele_hap(data, nbr_chr, foc_mrk, foc_all, phased, hap, &nbr_hap, nbr_chr_with_hap); //initialize the array of haplotypes (for the focal SNP) that contain the `allele'

	for (int i = 0; i < nbr_hap; i++) {
	  tot_nbr_chr_in_hap[foc_mrk] += nbr_chr_with_hap[i]; //total number of chromosomes; changes only by missing data
	}

	//if unphased, sequences from different individuals are different by definition
	int nbr_homo_haplo = phased ? tot_nbr_chr_in_hap[foc_mrk] - nbr_hap + 1 : (tot_nbr_chr_in_hap[foc_mrk] - nbr_hap) * 2;
	if(nbr_homo_haplo < lim_homo_haplo){
	  //do nothing, if too few homozygous chromosomes
	} else if (tot_nbr_chr_in_hap[foc_mrk] >= lim_haplo) { //if the total number of chromosomes is greater than minimum
		ehh[foc_mrk] = 1.0; //EHH at the focal marker is 1.0 for each allele, by definition

		extend_ehh(data, nbr_chr, foc_mrk, 0, lim_haplo, lim_homo_haplo, lim_ehh, phased, hap, &nbr_hap, nbr_chr_with_hap,
				tot_nbr_chr_in_hap, ehh);

		init_allele_hap(data, nbr_chr, foc_mrk, foc_all, phased, hap, &nbr_hap, nbr_chr_with_hap); // initialize again

		extend_ehh(data, nbr_chr, foc_mrk, nbr_mrk - 1, lim_haplo, lim_homo_haplo, lim_ehh, phased, hap, &nbr_hap, nbr_chr_with_hap,
				tot_nbr_chr_in_hap, ehh);
	}
	
	free(hap);
	free(nbr_chr_with_hap);
}
