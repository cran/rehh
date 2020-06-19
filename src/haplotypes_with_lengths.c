#include "definitions.h"
#include "haplotypes_with_lengths.h"


/**
 */
int update_hap_with_lengths(const int* const data, const int nbr_chr, const int mrk, int* const hap, int* const nbr_hap,
               int* const nbr_chr_with_hap, const double length, double* const hap_length) {
  
  int index_in_hap = 0;  // position within hap array
  const int mrk_index = mrk * nbr_chr; // position of the marker column in the matrix regarded as linear array
  short int tree_changed = 0; //"boolean" value that notifies if tree structure has changed
  
  // Loop over all distinct haplotypes
  for (int i = 0; i < *nbr_hap; i++) {
    if (nbr_chr_with_hap[i] == 1) {
      //haplogroup cannot change any more; set index to next haplogroup
      index_in_hap++;
    } else {
      // Loop over the chromosomes with the i-th haplotype
      // remove all chromosomes with a missing value at the current marker
      for (int j = index_in_hap; j < index_in_hap + nbr_chr_with_hap[i]; j++) {
        if(data[mrk_index + hap[j]] == MISSING_VALUE) {
          
          //set shared lengths with all other chromosomes of the same haplogroup
          for(int k = index_in_hap; k < index_in_hap + nbr_chr_with_hap[i]; k++){
            if(k != j){
              int chr1 = hap[j];
              int chr2 = hap[k];
              hap_length[chr1 * nbr_chr + chr2] += length;
              hap_length[chr2 * nbr_chr + chr1] += length;
            }
          }
          
          //remove chromosome from the hap array by shifting the remainder indices leftwards
          for (int k = j; k < nbr_chr - 1; k++) {
            hap[k] = hap[k + 1];
          }
          nbr_chr_with_hap[i]--;
          j--;
          //"pseudo" furcation (new "node" with missing value)
          tree_changed = 1;
        }
      }
      if (nbr_chr_with_hap[i] == 0) { //whole haplogroup has missing values
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
        
        // loop over the sorted chromosomes
        int j = 1;
        while(j < nbr_chr_with_hap[i]) {
          //do *subsequent* chromosomes within the haplotype group have different alleles at the marker?
          if (data[mrk_index + hap[index_in_hap + j - 1]] != data[mrk_index + hap[index_in_hap + j]]) {
            //new furcation
            tree_changed = 1;
            //shift upper indices to the right
            for (int k = *nbr_hap; k > i + 1; k--) {
              nbr_chr_with_hap[k] = nbr_chr_with_hap[k - 1];
            }
            //split current number of haplotypes to the two "daughter nodes"
            nbr_chr_with_hap[i + 1] = nbr_chr_with_hap[i] - j;
            nbr_chr_with_hap[i] = j;
            
            for(int k=0; k < nbr_chr_with_hap[i]; k++){
              for(int l=0; l < nbr_chr_with_hap[i + 1]; l++){
                int chr1 = hap[index_in_hap + k];
                int chr2 = hap[index_in_hap + nbr_chr_with_hap[i] + l];
                hap_length[chr1 * nbr_chr + chr2] += length;
                hap_length[chr2 * nbr_chr + chr1] += length;
              }
            }
            //increase number of haplogroups
            (*nbr_hap)++;
            //set indices to next haplogroup
            index_in_hap += nbr_chr_with_hap[i];	
            i++;
            //reset index within new haplogroup
            j = 1;
          }else{
            //next index in same haplogroup
            j++;
          }
        }
        //set index to next haplogroup
        index_in_hap += nbr_chr_with_hap[i];
      }
    }
  }
  
#ifdef DEBUG
  print_haplotypes(mrk, hap, nbr_hap, nbr_chr_with_hap);
#endif
  return (tree_changed);
}
