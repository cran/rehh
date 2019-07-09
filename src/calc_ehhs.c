#include "calc_ehhs.h"
#include "haplotypes.h"
#include "homozygosity.h"

/**
 * Calculates ehhs values in increasing distance (either in positive or negative direction) form the focal marker.
 * EHHS Tang is equal to EHHS Sabeti normalized by the homozygosity at the focal marker
 */
void extend_ehhs(const int* const data, const int nbr_chr, const int foc_mrk, const int end_mrk, const int lim_haplo,
		const double lim_ehhs, const int phased, int* const hap, int* const nbr_hap, int* const nbr_chr_with_hap,
		int* const tot_nbr_chr_in_hap, double* const ehhs, double* const nehhs) {

	int increment = foc_mrk <= end_mrk ? 1 : -1;

	for (int mrk = foc_mrk + increment; mrk != end_mrk + increment; mrk += increment) { //walk along the chromosome, away from the focal SNP
		if (update_hap(data, nbr_chr, mrk, hap, nbr_hap, nbr_chr_with_hap)) { //if there is a change in haplotypes (i.e. at least one furcation)
			for (int i = 0; i < *nbr_hap; i++) {
				tot_nbr_chr_in_hap[mrk] += nbr_chr_with_hap[i]; //total number of chromosomes; changes only by missing data
			}
			if (tot_nbr_chr_in_hap[mrk] < lim_haplo) {
				break; //stop, if too few chromosomes left due to missing values
			}
			ehhs[mrk] = homozygosity(tot_nbr_chr_in_hap[mrk], *nbr_hap, nbr_chr_with_hap, phased); // Compute the standardized hapotype homozygosity at the jth SNP, following Sabeti et al's (2007)

			if (phased) {
				/* calculate number of alleles at the focal marker in order to obtain EHHS Tang
				 * this considers only remaining chromosomes (i.e. without missing values up to marker 'mrk')
				 * the numbers may hence change at any newly considered marker with missing values!
				 */
				int *nbr_chr_with_foc_all = (int *) malloc(*nbr_hap * sizeof(int));
				int *foc_all = (int *) malloc(*nbr_hap * sizeof(int));
				int nbr_dstnct_foc_all = 0;
				int index = 0;
				for (int i = 0; i < *nbr_hap; i++) {                       //loop over all distinct haplotypes
					int allele_found = 0;
					for (int j = 0; j < nbr_dstnct_foc_all; j++) {
						//read as data[hap[index]][start_mrk]
						if (data[foc_mrk * nbr_chr + hap[index]] == foc_all[j]) {
							nbr_chr_with_foc_all[j] += nbr_chr_with_hap[i];
							allele_found = 1;
							break;
						}
					}
					if (!allele_found) {
						//read as data[hap[index]][start_mrk]
						foc_all[nbr_dstnct_foc_all] = data[foc_mrk * nbr_chr + hap[index]];
						nbr_chr_with_foc_all[nbr_dstnct_foc_all] = nbr_chr_with_hap[i];
						nbr_dstnct_foc_all++;
					}
					index += nbr_chr_with_hap[i];
				}
				free(foc_all);

				nehhs[mrk] = ehhs[mrk]
						/ homozygosity(tot_nbr_chr_in_hap[mrk], nbr_dstnct_foc_all, nbr_chr_with_foc_all, phased);

				free(nbr_chr_with_foc_all);

			} else {
				//for un-phased data, the homozygosity at the focal marker is always 1
				nehhs[mrk] = ehhs[mrk];
			}

			if (ehhs[mrk] <= lim_ehhs) { //if homozygosity is small
				ehhs[mrk] = 0.0;
				//note that EHHS Tang is always equal or greater than EHHS Sabeti
				if (nehhs[mrk] <= lim_ehhs) {
					nehhs[mrk] = 0.0;
					break;
				}
			}
		} else { //if no change in the number of haplotypes, then re-use previous values
			tot_nbr_chr_in_hap[mrk] = tot_nbr_chr_in_hap[mrk - increment];
			nehhs[mrk] = nehhs[mrk - increment];
			ehhs[mrk] = ehhs[mrk - increment];
		}
	}
}

/**
 * Calculates ehhs values in increasing distance to the focal marker
 */
void calc_ehhs(const int* const data, const int nbr_chr, const int nbr_mrk, const int foc_mrk, const int lim_haplo,
		const double lim_ehhs, const int phased, int* const tot_nbr_chr_in_hap, double* const ehhs,
		double* const nehhs) {

	int nbr_hap;                                                   //number of distinct haplotypes

	int *hap = (int*) malloc(nbr_chr * sizeof(int));               //vector index to chromosomes, ordered by haplotype
	int *nbr_chr_with_hap = (int *) malloc(nbr_chr * sizeof(int)); //for each haplotype gives the number of chromosomes sharing it

	//initialization of return values
	for (int mrk = 0; mrk < nbr_mrk; mrk++) {
		ehhs[mrk] = 0.0;
		nehhs[mrk] = 0.0;
		tot_nbr_chr_in_hap[mrk] = 0;
	}

	init_site_hap(data, nbr_chr, foc_mrk, phased, hap, &nbr_hap, nbr_chr_with_hap); // initialize the array of haplotypes

	for (int i = 0; i < nbr_hap; i++) {
		tot_nbr_chr_in_hap[foc_mrk] += nbr_chr_with_hap[i];
	}

	if (tot_nbr_chr_in_hap[foc_mrk] >= lim_haplo) { //if the total number of chromosomes is greater than minimum
		nehhs[foc_mrk] = 1.0; //normalized EHHS is 1.0 at the focal marker by definition
		ehhs[foc_mrk] = homozygosity(tot_nbr_chr_in_hap[foc_mrk], nbr_hap, nbr_chr_with_hap, phased);

		extend_ehhs(data, nbr_chr, foc_mrk, 0, lim_haplo, lim_ehhs, phased, hap, &nbr_hap, nbr_chr_with_hap,
				tot_nbr_chr_in_hap, ehhs, nehhs);

		init_site_hap(data, nbr_chr, foc_mrk, phased, hap, &nbr_hap, nbr_chr_with_hap); // initialize again

		extend_ehhs(data, nbr_chr, foc_mrk, nbr_mrk - 1, lim_haplo, lim_ehhs, phased, hap, &nbr_hap, nbr_chr_with_hap,
				tot_nbr_chr_in_hap, ehhs, nehhs);
	}
	free(hap);
	free(nbr_chr_with_hap);
}
