#include "definitions.h"
#include "haplotypes_with_nodes.h"
#include "calc_furcation.h"

/**
 * Extends furcation tree on either left or right side of the focal marker.
 */
void extend_furcation(const int* const data, const int nbr_chr, const int foc_mrk, const int end_mrk,
                      const int lim_haplo, int* const hap, int* const nbr_hap, int* const nbr_chr_with_hap, 
                      int* const hap_node, int* const node_size, int* const node_mrk, int* const node_parent, 
                      int* const node_with_missing_data, int* const nbr_node, int* const label_parent) {
  
  int increment = foc_mrk <= end_mrk ? 1 : -1;
  for (int mrk = (foc_mrk + increment); mrk != end_mrk + increment; mrk += increment) { // walk along the chromosome, away from the focal SNP
    if (update_hap_with_nodes(data, nbr_chr, mrk, hap, nbr_hap, nbr_chr_with_hap, 
                              hap_node, node_size, node_mrk, node_parent, node_with_missing_data, nbr_node, label_parent)) {
      int tot_nbr_chr_in_hap = 0;
      for (int i = 0; i < *nbr_hap; i++) {
        tot_nbr_chr_in_hap += nbr_chr_with_hap[i];
      }
      
      //stop if two few haplotypes left or all haplotypes are distinct
      if ((tot_nbr_chr_in_hap < lim_haplo) | (tot_nbr_chr_in_hap == *nbr_hap)) {
        break;
      }
    }
  }
}

/**
 * Computes furcation trees on both sides of the focal marker.
 */
void calc_furcation(const int* const data, const int nbr_chr, const int foc_mrk, const int end_mrk, const int foc_all,
                    const int lim_haplo, const int phased, int* const node_mrk, int* const node_parent,
                    int* const node_with_missing_data, int* const nbr_node, int* const label_parent) {
  
  int nbr_hap;                                                  //number of distinct haplotypes
  
  int *hap = (int*) malloc(nbr_chr * sizeof(int));              //vector index to chromosomes, ordered by haplotype
  int *nbr_chr_with_hap = (int*) malloc(nbr_chr * sizeof(int)); //for each haplotype gives the number of chromosomes sharing it
  int *node_size = (int*) malloc((nbr_chr * 2 - 1) * sizeof(int));//number of children nodes; at most 2n-1 nodes in a tree of n leaves
  int *hap_node = (int*) malloc(nbr_chr * sizeof(int));         //remember for each haplo-group the corresponding node
  
  //initialization of return values
  for (int chr = 0; chr < nbr_chr; chr++) {
    label_parent[chr] = NA_INTEGER; // no label node for chromosome 'chr'
  }
  for (int i = 0; i < 2 * nbr_chr - 1; i++) {
    node_mrk[i] = NA_INTEGER; // no marker associated with node i
    node_parent[i] = NA_INTEGER; // no parent associated with node i
    node_with_missing_data[i] = 0; // no missing data associated with node i (0 = FALSE)
  }
  
  init_allele_hap_with_nodes(data, nbr_chr, foc_mrk, foc_all, phased, hap, &nbr_hap, nbr_chr_with_hap,
                             hap_node, node_size, node_mrk, node_parent, nbr_node); 
  
  if (nbr_chr_with_hap[0] >= lim_haplo) {
    extend_furcation(data, nbr_chr, foc_mrk, end_mrk, lim_haplo, hap, &nbr_hap, nbr_chr_with_hap,
                     hap_node, node_size, node_mrk, node_parent, node_with_missing_data, nbr_node, label_parent);
    
    int index = 0;
    for (int i = 0; i < nbr_hap; i++) {
      for(int j = index; j < index + nbr_chr_with_hap[i]; j++){
        label_parent[hap[j]] = hap_node[i];
      }
      index += nbr_chr_with_hap[i];
    }
    
  } else {
    *nbr_node = 0;
  }
  free(hap);
  free(nbr_chr_with_hap);
  free(node_size);
  free(hap_node);
}
