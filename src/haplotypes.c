#include "haplotypes.h"

/**
 * The data structure that describes haplotypes has been reversed with respect to the previous release.
 * Instead of each chromosom referring to a haplotype (via a certain chromosome designating as "reference"),
 * now the haplotypes refer to chromosomes. Logically, each haplotype is associated with a vector of
 * variable size containing the indices of the chromosomes. For performance reasons this information is
 * stored consecutively in a single vector and the number of chromosomes belonging to each haplotype is
 * stored in another vector of.
 * Example:
 * Assume chromosomes 1,2,5 are identical, chromosomes 0,3 are identical and chromosome 4 is different
 * from the other two haplotypes. Then the number of distinct haplotypes is 3 and their
 * respective indices are 0, 1 and 2.
 * This might be represented by the vector 'hap' as (1,2,5,4,0,3) and the vector nbr_chr_with_hap yields (3,1,2).
 * The haplotypes are ordered, roughly speaking, by their amount of ancestral and derived alleles; details are
 * given below.
 */

#ifdef DEBUG
/**
 * Only for debugging. Prints the index(+1) of the marker and the indices(+1) of chromosomes, ordered by haplotype.
 */
void print_haplotypes(const int mrk, const int *hap, const int *nbr_hap, const int *nbr_chr_with_hap) {
	printf("Haplotypes at marker %i: ", mrk + 1);
  printf("Number of distinct haplotypes : %i\n",*nbr_hap);

	int counter = 0;
	for (int i = 0; i < *nbr_hap; i++) {
		for (int j = counter; j < counter + nbr_chr_with_hap[i]; j++) {
			printf("%i ", hap[j] + 1);
		}
		counter += nbr_chr_with_hap[i];
		if(i < *nbr_hap - 1) {
			printf("| ");
		} else {
			printf("\n");
		}
	}
}
#endif

/**
 * Consider only chromosomes with the specified allele at the focal marker. This subset of all chromsomes
 * is thus be definition homozygous at the focal marker.
 * Example:
 * If chromosomes 0,3,4 have the allele, while chromosomes 1,2,5 don't, we get a single haplotype
 * with index 0, the vector 'hap' yields (0,3,4,x,x,x), nbr_chr_with_hap=(3) and nbr_hap=1.
 */
void init_allele_hap(const int *data, const int nbr_chr, const int foc_mrk, const int foc_all, const int phased,
		int* const hap, int* const nbr_hap, int* const nbr_chr_with_hap) {

	*nbr_hap = 0;
	nbr_chr_with_hap[0] = 0;                                           // Initialize the haplotype counts

#ifdef DEBUG
	printf("### Allele: %i\n", foc_all);
#endif

	if (!phased) {
		for (int i = 0; i < nbr_chr - 1; i += 2) {
			//read as data[i][foc_mrk], data[i+1][foc_mrk]
			if (data[foc_mrk * nbr_chr + i] == foc_all && data[foc_mrk * nbr_chr + i + 1] == foc_all) {
        hap[(*nbr_hap) * 2] = i;
			  hap[(*nbr_hap) * 2 + 1] = i + 1;
				nbr_chr_with_hap[*nbr_hap] = 2;
				(*nbr_hap)++;
			}
		}
	} else {
		for (int i = 0; i < nbr_chr; i++) {								// Loop over all chromosomes in the sample
			//read as data[i][foc_mrk]
			if (data[foc_mrk * nbr_chr + i] == foc_all) {
				hap[nbr_chr_with_hap[0]] = i;
				nbr_chr_with_hap[0]++;
			}
		}

		if (nbr_chr_with_hap[0] > 0) {
			*nbr_hap = 1;
		}
	}

#ifdef DEBUG
	print_haplotypes(foc_mrk, hap, nbr_hap, nbr_chr_with_hap);
#endif
}

/**
 * Consider all chromosomes which have no missing alleles at the focal marker. If the marker is a SNP,
 * thus polymorphic, by definition the chromosomes are not homozygous. Note that in comparison
 * among populations also monomorphic sites are considered!
 * The algorithm for phased data starts with a fictitious monomorphic site on which the update
 * algorithm is applied for the focal site.
 */
void init_site_hap(const int *data, const int nbr_chr, const int foc_mrk, const int phased, int* const hap,
		int* const nbr_hap, int* const nbr_chr_with_hap) {

	if (!phased) {
		*nbr_hap = 0;
		nbr_chr_with_hap[0] = 0;                                           // Initialize the haplotype counts
		for (int i = 0; i < nbr_chr - 1; i += 2) {
			//read as data[i][foc_mrk]==data[i+1][foc_mrk]
			if (data[foc_mrk * nbr_chr + i] != MISSING_VALUE
					&& data[foc_mrk * nbr_chr + i] == data[foc_mrk * nbr_chr + i + 1]) {
				hap[(*nbr_hap) * 2] = i;
				hap[(*nbr_hap) * 2 + 1] = i + 1;
				nbr_chr_with_hap[*nbr_hap] = 2;
				(*nbr_hap)++;
			}
		}
	} else {
		for (int i = 0; i < nbr_chr; i++) {
			hap[i] = i;
		}
		nbr_chr_with_hap[0] = nbr_chr;
		*nbr_hap = 1;
		update_hap(data, nbr_chr, foc_mrk, hap, nbr_hap, nbr_chr_with_hap);
	}
}

/**
 * Updates given haplotypes after inspecting alleles at a certain given marker.
 * If chromosomes of the same given haplotype have different alleles at the marker, the haplotype is
 * devided into new sub-haplotypes. These are stored in the same positions within the vector 'hap'
 * as their "parent".
 * Example:
 * Given a vector hap=(1,2,5,4,0,3), nbr_hapl=3 and nbr_chr_with_hap (3,1,2).
 * Hence we have three different haplotypes, the first referring to chromosomes 1,2,5
 * the second to chromosome 4 and the last two to chromosomes 0 and 3.
 * Assume chromosomes 1 and 5 have one allele at the marker and chromosome 2 another and that chromosomes
 * 0 and 3 have different alleles a that site, too. Then the update may yield
 * hap=(1,5,2,4,0,3), nbr_chr_with_hap=(2,1,1,1,1) and obviously we have five distinct haplotypes.
 * The haplotypes are ordered with respect to the allele, more precisely, the respective representation,
 * i.e. 0 for ancestral and 1,2,3,... for derived variants.
 *
 * Chromosomes with missing data at the marker are removed. Note, that for the normalized EHHS,
 * this removal affects the normalization at the focal marker, too.
 *
 * Already existing "singleton"-haplotypes are not changed any further, i.e. neither subdivided nor removed.
 */
int update_hap(const int *data, const int nbr_chr, const int mrk, int* const hap, int* const nbr_hap,
		int* const nbr_chr_with_hap) {

	int index_in_hap = 0;  // current index within hap
	int *sub_nbr_chr_with_hap = (int *) malloc(nbr_chr * sizeof(int)); // the number of identical hapotypes within each group
	int *dstnct_alleles = (int *) malloc(nbr_chr * sizeof(int));	// the distinct alleles at the current mrk
	const int mrk_index = mrk * nbr_chr; // position of the marker column in the matrix regarded as linear array

	short int tree_changed = 0;

	// Look over all distinct haplotypes
	for (int i = 0; i < *nbr_hap; i++) {
		// Loop over the chromosomes with the i-th haplotype
		// If the i-th haplotype is present only on a single chromosome, nothing can change any more; skip it.
		if (nbr_chr_with_hap[i] > 1) {
			// remove all chromosomes with missing alleles at the marker
			for (int j = index_in_hap; j < index_in_hap + nbr_chr_with_hap[i]; j++) {
				if (data[mrk_index + hap[j]] == MISSING_VALUE) {
					// delete chromosomes by shifting the remaineder leftwards
					for (int k = j; k < nbr_chr - 1; k++) {
						hap[k] = hap[k + 1];
					}
					nbr_chr_with_hap[i]--;
					j--;
					tree_changed = 1;
				}
			}
			if (nbr_chr_with_hap[i] == 0) { // implies missing values found
				// empty set, hence remove haplotype group
				for (int k = i; k < *nbr_hap - 1; k++) {
					nbr_chr_with_hap[k] = nbr_chr_with_hap[k + 1];
				}
				(*nbr_hap)--;
	    } else {
				// Sort the chromosomes with respect to the allele at the current marker ("Insertion sort")
				// important: this sort is conservative, it does not change the order of equal elements
				for (int j = index_in_hap + 1; j < index_in_hap + nbr_chr_with_hap[i]; j++) {
					int tmp = hap[j];
					int k = j;
					//read as data[hap[k-1]][mrk] > data[hap[j]][mrk]
					while (k > index_in_hap && data[mrk_index + hap[k - 1]] > data[mrk_index + tmp]) {
						hap[k] = hap[k - 1];
						k--;
					}
					hap[k] = tmp;
				}

				// now lets partition the sorted haplotypes
				int sub_nbr_hap = 1;

				sub_nbr_chr_with_hap[0] = 1;
				// loop over the same chromosomes as before
				for (int j = index_in_hap + 1; j < index_in_hap + nbr_chr_with_hap[i]; j++) {
					//read as data[hap[j]][mrk] != data[hap[j-1]][mrk]
					//if *subsequent* chromosomes have different alleles at the marker,
					//define a new subgroup
					if (data[mrk_index + hap[j]] != data[mrk_index + hap[j - 1]]) {
						sub_nbr_hap++;
						sub_nbr_chr_with_hap[sub_nbr_hap - 1] = 1;
					} else {
						sub_nbr_chr_with_hap[sub_nbr_hap - 1]++;
					}
				}

				int old_nbr_chr_with_hap = nbr_chr_with_hap[i];

				//if there is a furcation, insert new haps
				if (sub_nbr_hap > 1) {
					//signal tree structure change
					tree_changed = 1;
					//shift all groups to the right and make place for new subgroups
					for (int j = *nbr_hap - 1; j > i; j--) {
						nbr_chr_with_hap[j + sub_nbr_hap - 1] = nbr_chr_with_hap[j];
					}
					//insert new subgroups
					for (int j = 0; j < sub_nbr_hap; j++) {
						nbr_chr_with_hap[i + j] = sub_nbr_chr_with_hap[j];
					}
					//shift index
					i += sub_nbr_hap - 1;
					//add number of new sub groups to total number of groups
					(*nbr_hap) += sub_nbr_hap - 1;
				}

				index_in_hap += old_nbr_chr_with_hap;
			}
		} else {
			index_in_hap++;
		}
	}

	free(dstnct_alleles);
	free(sub_nbr_chr_with_hap);
#ifdef DEBUG
	print_haplotypes(mrk, hap, nbr_hap, nbr_chr_with_hap);
#endif
	return (tree_changed);
}
