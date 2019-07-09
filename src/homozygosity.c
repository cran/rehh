#include "homozygosity.h"

/**
 * tot_nbr_chr_in_hap refers to the total number of chromosomes considered
 * nbr_hap is the number of distinct haplotypes (classes of identical chromosomes)
 * nbr_chr_in_hap is a vector which gives for each distinct haplotype the number of
 * chromosomes belonging to this haplotype (i.e. which form a group of identical chromosomes)
 * the sum over all nbr_chr_in_hap should be equal to tot_nbr_chr_in_hap; since the latter
 * is calculated already elsewhere (because it is evaluated as stopping condition), it is given
 * as parameter although it could be calculated from the other two parameters
 *
 * phased should be boolean
 */
double homozygosity(int tot_nbr_chr_in_hap, int nbr_hap, int *nbr_chr_in_hap, int phased) {
	double homzgsty = 0.0;

	if (tot_nbr_chr_in_hap > 1) {
		if (phased) {
			for (int i = 0; i < nbr_hap; i++) {  // Loop over all distinct haplotypes
				homzgsty += (double) nbr_chr_in_hap[i] * (nbr_chr_in_hap[i] - 1.0);
			}
			homzgsty /= (double) (tot_nbr_chr_in_hap * (tot_nbr_chr_in_hap - 1.0));

		} else {
			homzgsty = (tot_nbr_chr_in_hap - nbr_hap) * 2 / (double) tot_nbr_chr_in_hap;
		}
	}
	return (homzgsty);
}
