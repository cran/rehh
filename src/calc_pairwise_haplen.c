#include "definitions.h"
#include "haplotypes.h"
#include "haplotypes_with_lengths.h"
#include "calc_pairwise_haplen.h"

/**
 * Extends haplength on either left or right side of the focal marker.
 */
int extend_haplen(const int* const data, const int nbr_chr, const double* map, const int foc_mrk, const int end_mrk,
                  int* const hap, int* const nbr_hap, int* const nbr_chr_with_hap,
                  const int max_gap, const bool discard_integration_at_border, double* const pairwise_haplen) {
  
  int increment = foc_mrk <= end_mrk ? 1 : -1;
  int mrk;
  for (mrk = foc_mrk + increment; mrk != end_mrk + increment; mrk += increment) { // walk along the chromosome, away from the focal SNP
    double gap = increment * (map[mrk] - map[mrk - increment]);
    if(gap > max_gap){
      if(discard_integration_at_border){
        return(1);
      }
      break;
    }
    double length = increment * (map[mrk] - map[foc_mrk]);
    if (update_hap_with_lengths(data, nbr_chr, mrk, hap, nbr_hap, nbr_chr_with_hap, length, pairwise_haplen)) {
      
      int tot_nbr_chr_in_hap = 0;
      for (int i = 0; i < *nbr_hap; i++) {
        tot_nbr_chr_in_hap += nbr_chr_with_hap[i];
      }
      
      //stop if all haplotypes are distinct
      if (tot_nbr_chr_in_hap == *nbr_hap) {
        break;
      }
    }
  }
  
  mrk -= increment; //set back to last index in loop
  if(mrk == end_mrk && discard_integration_at_border){
    return(1);
  }
  
  int index_in_hap = 0;
  for(int i=0; i < *nbr_hap; i++){
    if(nbr_chr_with_hap[i] > 1){
      if(discard_integration_at_border){
        return(1);
      }
      for(int k = 1; k < nbr_chr_with_hap[i]; k++){
        for(int l = 0; l < k; l++){
          int chr1 = hap[index_in_hap + k];
          int chr2 = hap[index_in_hap + l];
          pairwise_haplen[chr1 * nbr_chr + chr2] += increment * (map[mrk] - map[foc_mrk]);
          pairwise_haplen[chr2 * nbr_chr + chr1] += increment * (map[mrk] - map[foc_mrk]);
        }
      }
    }
    index_in_hap += nbr_chr_with_hap[i];
  }
  
  return(0);
}

/**
 * Computes furcation trees on both sides of the focal marker.
 */
int calc_pairwise_haplen(const int* const data, const int nbr_chr, const int nbr_mrk, 
                         const double* map, const int foc_mrk, const int max_gap,
                         const bool phased, const bool discard_integration_at_border, 
                         double* const pairwise_haplen) {
  
  int nbr_hap;                                                  //number of distinct haplotypes
  bool discard = false; 
  
  int *hap = (int*) malloc(nbr_chr * sizeof(int));              //vector index to chromosomes, ordered by haplotype
  int *nbr_chr_with_hap = (int*) malloc(nbr_chr * sizeof(int)); //for each haplotype gives the number of chromosomes sharing it
  
  init_site_hap(data, nbr_chr, foc_mrk, phased, hap, &nbr_hap, nbr_chr_with_hap); 
  
  discard = extend_haplen(data, nbr_chr, map, foc_mrk, 0, hap, &nbr_hap, nbr_chr_with_hap,
                             max_gap, discard_integration_at_border, pairwise_haplen);
  if(!discard){
    init_site_hap(data, nbr_chr, foc_mrk, phased, hap, &nbr_hap, nbr_chr_with_hap); 
    discard = extend_haplen(data, nbr_chr, map, foc_mrk, nbr_mrk - 1, hap, &nbr_hap, nbr_chr_with_hap,
                               max_gap, discard_integration_at_border, pairwise_haplen);
  }
  
  free(hap);
  free(nbr_chr_with_hap);
  return(discard ? 1 : 0);
}
