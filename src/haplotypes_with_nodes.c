#include "definitions.h"
#include "haplotypes_with_nodes.h"

/**
 *Extensions of the corresponding functions in haplotypes.c to build the tree structure of furcations.
 *For performance reasons this code is used only when the furcation structure is actually needed.
 */


void init_allele_hap_with_nodes(const int *data, const int nbr_chr, const int foc_mrk, const int foc_all, const int phased,
                                int* const hap, int* const nbr_hap, int* const nbr_chr_with_hap,
                                int* const hap_node, int* const node_size, 
                                int* const node_mrk, int* const node_parent, int* const nbr_node) {
  
  *nbr_hap = 0;  // current haplogroup has number 0
  nbr_chr_with_hap[0] = 0;  // haplogroup 0 has 0 elements
  
#ifdef DEBUG
  printf("### Allele: %i\n", foc_all);
#endif
  
  if (!phased) {
    // run through all (assumed) diploid individuals
    for (int i = 0; i < nbr_chr - 1; i += 2) {
      //check if individual is homozygous for the focal allele
      //read as data[i][foc_mrk], data[i+1][foc_mrk]
      if (data[foc_mrk * nbr_chr + i] == foc_all && data[foc_mrk * nbr_chr + i + 1] == foc_all) {
        //current haplogroup points to the two chromosomes of the individual
        hap[(*nbr_hap) * 2] = i;
        hap[(*nbr_hap) * 2 + 1] = i + 1;
        //size of haplogroup
        nbr_chr_with_hap[*nbr_hap] = 2;
        //next haplogroup
        (*nbr_hap)++;
      }
    }
  } else {
    for (int i = 0; i < nbr_chr; i++) {								// Loop over all chromosomes in the sample
      //read as data[i][foc_mrk]
      if (data[foc_mrk * nbr_chr + i] == foc_all) {
        //add all chromosomes carrying the focal allele to a single haplogroup
        hap[nbr_chr_with_hap[0]] = i;
        nbr_chr_with_hap[0]++;
      }
    }
    
    if (nbr_chr_with_hap[0] > 0) {
      *nbr_hap = 1;
    }
  }
  
  node_mrk[0] = foc_mrk; // root node marker points to focal marker
  *nbr_node = 1; // start with one node (the root)
  if (phased) {
    hap_node[0] = 0;
    node_size[0] = nbr_chr_with_hap[0];
  } else {
    int index = 0;
    // the root node gets a daughter node for each (homozygous) individual
    for (int i = 0; i < *nbr_hap; i++) { //run through all haplogroups (=chromosomes of individuals)
      index += nbr_chr_with_hap[i]; // goto next haplogroup
      node_mrk[*nbr_node] = foc_mrk; // all node markers point to focal marker
      node_parent[*nbr_node] = 0; // all nodes have the root as parent
      hap_node[i] = *nbr_node;
      node_size[*nbr_node] = nbr_chr_with_hap[i];
      node_size[0] += nbr_chr_with_hap[i]; //root node size
      (*nbr_node)++; // create new node
    }
  }
  
#ifdef DEBUG
  print_haplotypes(foc_mrk, hap, nbr_hap, nbr_chr_with_hap);
#endif
}

int update_hap_with_nodes(const int *data, const int nbr_chr, const int mrk, int* const hap, int* const nbr_hap,
                          int* const nbr_chr_with_hap, int* const hap_node, int* const node_size,  int* const node_mrk, 
                          int* const node_parent, int* const node_with_missing_data, int* const nbr_node, int* const label_parent) {
  
  int index_in_hap = 0;  // position within hap array
  const int mrk_index = mrk * nbr_chr; // position of the marker column in the matrix regarded as linear array
  short int tree_changed = 0; //"boolean" value that notifies if tree structure has changed
  
  // Loop over all distinct haplotypes
  for (int i = 0; i < *nbr_hap; i++) {
    if (nbr_chr_with_hap[i] == 1) {
      //haplogroup cannot change any more; set index to next haplogroup
      index_in_hap++;
    } else {
      int nbr_missing = 0;
      
      // Loop over the chromosomes with the i-th haplotype
      // remove all chromosomes with a missing value at the current marker
      for (int j = index_in_hap; j < index_in_hap + nbr_chr_with_hap[i]; j++) {
        if (data[mrk_index + hap[j]] == MISSING_VALUE) {
          label_parent[hap[j]] = *nbr_node;
          
          //remove chromosome from the hap array by shifting the remainder indices leftwards
          for (int k = j; k < nbr_chr - 1; k++) {
            hap[k] = hap[k + 1];
          }
          
          nbr_chr_with_hap[i]--;
          j--;
          nbr_missing++;
          //"pseudo" furcation (new "node" with missing value)
          tree_changed = 1;
        }
      }
      if (nbr_chr_with_hap[i] == 0) { //whole haplogroup has missing values
        //add a single new node
        node_mrk[*nbr_node] = mrk;
        node_parent[*nbr_node] = hap_node[i];
        node_size[*nbr_node] = nbr_chr_with_hap[i];
        node_with_missing_data[*nbr_node] = 1;
        hap_node[i] = *nbr_node;
        (*nbr_node)++;
        
        //remove haplotype group
        for (int k = i; k < *nbr_hap - 1; k++) {
          nbr_chr_with_hap[k] = nbr_chr_with_hap[k + 1];
          hap_node[k] = hap_node[k + 1];
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
        int nbr_new_branches = 0;
        
        while(j < nbr_chr_with_hap[i]) {
          //do *subsequent* chromosomes within the haplotype group have different alleles at the marker?
          if (data[mrk_index + hap[index_in_hap + j - 1]] != data[mrk_index + hap[index_in_hap + j]]) {
            //new furcation
            tree_changed = 1;
            //shift upper indices to the right
            for (int k = *nbr_hap; k > i + 1; k--) {
              nbr_chr_with_hap[k] = nbr_chr_with_hap[k - 1];
              hap_node[k] = hap_node[k - 1];
            }
            //split current number of haplotypes to the two "daughter nodes"
            nbr_chr_with_hap[i + 1] = nbr_chr_with_hap[i] - j;
            nbr_chr_with_hap[i] = j;
            //increase number of haplogroups
            (*nbr_hap)++;
            //set indices to next haplogroup
            index_in_hap += nbr_chr_with_hap[i];	
            i++;
            nbr_new_branches++;
            //reset index within new haplogroup
            j = 1;
          }else{
            //next index in same haplogroup
            j++;
          }
        }
        
        if(nbr_new_branches + nbr_missing > 0){
          int parent = hap_node[i - nbr_new_branches];
          
          for(int j = 0; j <= nbr_new_branches; j++){  
            node_size[(*nbr_node) + j] = nbr_chr_with_hap[i + j];
            node_parent[(*nbr_node) + j] = parent;
            node_mrk[(*nbr_node) + j] = mrk;
            hap_node[i - nbr_new_branches + j] = (*nbr_node) + j;
          }
          (*nbr_node) += nbr_new_branches + 1;
          
          if(nbr_missing > 0){
            //shift label parents created above to next new node
            for(int k = 0; k < nbr_chr; k++){
              if(label_parent[k] == (*nbr_node) - nbr_new_branches - 1){
                label_parent[k] += nbr_new_branches + 1;
              }
            }
            //add a node for the haplotypes with missing values
            node_size[(*nbr_node)] = nbr_missing;
            node_parent[(*nbr_node)] = parent;
            node_mrk[(*nbr_node)] = mrk;
            node_with_missing_data[*nbr_node] = 1;
            (*nbr_node)++;
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
