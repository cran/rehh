#include "haplotypes.h"
#include "calc_furcation.h"

/**
 * Extends furcation tree on either left or right side of the focal marker.
 */
void extend_furcation(const int* const data, const int nbr_chr, const int foc_mrk, const int end_mrk,
                      const int lim_haplo, int* const hap, int* const nbr_hap, int* const nbr_chr_with_hap, int* const node_mrk,
                      int* const node_parent, int* const nbr_node, int* const label_parent) {
  
  int increment = foc_mrk <= end_mrk ? 1 : -1;
  
  int* node_size = (int*) malloc((nbr_chr * 2 - 1) * sizeof(int));
  int* inner_nodes = (int*) malloc((nbr_chr) * sizeof(int));
  
  node_size[0] = nbr_chr_with_hap[0];
  
  for (int i = 1; i < *nbr_hap; i++) {
    node_size[i] = nbr_chr_with_hap[i];
    node_size[0] += nbr_chr_with_hap[i];
  }
  
  for (int mrk = (foc_mrk + increment); mrk != end_mrk + increment; mrk += increment) { // walk along the chromosome, away from the focal SNP
    if (update_hap(data, nbr_chr, mrk, hap, nbr_hap, nbr_chr_with_hap)) {
      int tot_nbr_chr_in_hap = 0;
      for (int i = 0; i < *nbr_hap; i++) {
        tot_nbr_chr_in_hap += nbr_chr_with_hap[i];
      }
      // the (sub)partitioning of extended haplotypes is translated into a tree-structure
      // this is surprisingly awkward and could be done more elegantly in the "update_hap" function
      // however this would harm the performance of its more important goal of calculating EHH(S)
      
      // essentially this code looks for sub-partitions and adds nodes if one is found
      // if sub-partitions don't cover the "parent" partition, the remaining chromosomes
      // get their own node and are discarded as missing
      
      // unfortunately, if there are no sub-partitions at all of a partition, then that haplotype
      // has been entirely lost by missing values. Such "lost" haplotypes are checked by
      // the presence of inner nodes (=haplotypes of old partition) which have no sub-partitions
      
      // build a list of all inner nodes (node size > 1)
      int nbr_inner_nodes = 0;
      for (int i = 0; i < *nbr_node; i++) {
        if (node_size[i] > 1) {
          inner_nodes[nbr_inner_nodes] = i;
          nbr_inner_nodes++;
        }
      }
      int index_in_hap = 0;
      for (int i = 0; i < *nbr_hap; i++) {
        //get previously set "last node", i.e. label_parent of one chromosomes (the first) of the haplo group
        int hap_grp_parent = label_parent[hap[index_in_hap]];
        // remove that node from the list of inner nodes
        for (int j = 0; j < nbr_inner_nodes; j++) {
          if (inner_nodes[j] == hap_grp_parent) {
            for (int k = j; k < nbr_inner_nodes - 1; k++) {
              inner_nodes[k] = inner_nodes[k + 1];
            }
            nbr_inner_nodes--;
            break;
          }
        }
        //go to next haplotype
        int next_index_in_hap = index_in_hap + nbr_chr_with_hap[i];
        // if previous parent node size is different from current hap group size
        if (node_size[hap_grp_parent] != nbr_chr_with_hap[i]) {
          // link to label_parent is "re-grafted" to a new node (with number nbr_node)
          for (int k = index_in_hap; k < next_index_in_hap; k++) {
            label_parent[hap[k]] = *nbr_node;
          }
          // the new node gets the old node as parent
          node_parent[*nbr_node] = hap_grp_parent;
          // assign the marker number to the new node
          node_mrk[*nbr_node] = mrk;
          // assign a size (= size of current haplo group) to the new node
          node_size[*nbr_node] = nbr_chr_with_hap[i];
          // update counter ("ready" for the next node)
          (*nbr_node)++;
          
          //if last hap group has been reached or next hap group has a different parent node
          if (i == (*nbr_hap) - 1 || label_parent[hap[next_index_in_hap]] != hap_grp_parent) {
            //check if there are still label_parents referring to the "old" parent node
            //these are hence absent from current hap groups and must have a missing value
            //create a single new node for all chromosomes with the same parent node
            int already_found_missing = 0;
            // go through all chromosomes to check "dangling" links
            for (int k = 0; k < nbr_chr; k++) {
              // still a link to the "old" node?
              if (label_parent[k] == hap_grp_parent) {
                // if a dangling node was found...
                if (already_found_missing == 0) {
                  //create new node
                  node_parent[*nbr_node] = hap_grp_parent;
                  // assign marker number to node
                  node_mrk[*nbr_node] = mrk;
                  // assign a node size of 1
                  node_size[*nbr_node] = 1;
                  // "ready" for the next node
                  (*nbr_node)++;
                  already_found_missing = 1;
                } else {
                  //update node size
                  node_size[(*nbr_node) - 1]++;
                }
                // regraft label from old parent to the new node
                label_parent[k] = (*nbr_node) - 1;
              }
            }
          }
        }
        index_in_hap += nbr_chr_with_hap[i];
      }
      // tackle on inner nodes which have not been handled, i.e. completely lost by missing values
      for (int j = 0; j < nbr_inner_nodes; j++) {
        int already_found_missing = 0;
        for (int k = 0; k < nbr_chr; k++) {
          if (label_parent[k] == inner_nodes[j]) {
            if (already_found_missing == 0) {
              //create new node
              node_parent[*nbr_node] = inner_nodes[j];
              // assign marker number to node
              node_mrk[*nbr_node] = mrk;
              // assign a node size of 1
              node_size[*nbr_node] = 1;
              // "ready" for the next node
              (*nbr_node)++;
              already_found_missing = 1;
            } else {
              //update node size
              node_size[(*nbr_node) - 1]++;
            }
            // regraft label from old parent to new node
            label_parent[k] = (*nbr_node) - 1;
          }
        }
      }
      
      //stop if two few haplotypes left or all haplotypes are distinct
      if ((tot_nbr_chr_in_hap < lim_haplo) | (tot_nbr_chr_in_hap == *nbr_hap)) {
        break;
      }
    }
  }
  free(node_size);
}

/**
 * Computes furcation trees on both sides of the focal marker.
 */
void calc_furcation(const int* const data, const int nbr_chr, const int foc_mrk, const int end_mrk, const int foc_all,
                    const int lim_haplo, const int phased, int* const node_mrk, int* const node_parent,
                    int* const node_with_missing_data, int* const nbr_node, int* const label_parent) {
  
  int nbr_hap;                                                  //number of distinct haplotypes
  
  int *hap = (int*) malloc(nbr_chr * sizeof(int));             //vector index to chromosomes, ordered by haplotype
  int *nbr_chr_with_hap = (int*) malloc(nbr_chr * sizeof(int)); //for each haplotype gives the number of chromosomes sharing it
  
  //initialization of return values
  for (int chr = 0; chr < nbr_chr; chr++) {
    label_parent[chr] = NA_INTEGER;
  }
  for (int i = 0; i < 2 * nbr_chr - 1; i++) {
    node_mrk[i] = NA_INTEGER;
    node_parent[i] = NA_INTEGER;
    node_with_missing_data[i] = 0;
  }
  
  init_allele_hap(data, nbr_chr, foc_mrk, foc_all, phased, hap, &nbr_hap, nbr_chr_with_hap); // initialize the array of hapotypes (for the focal SNP) that contain the `allele'
  
  if (nbr_chr_with_hap[0] >= lim_haplo) {
    node_mrk[0] = foc_mrk;
    *nbr_node = 1;
    if (phased) {
      for (int j = 0; j < nbr_chr_with_hap[0]; j++) {
        label_parent[hap[j]] = 0;
      }
    } else {
      int index = 0;
      for (int i = 0; i < nbr_hap; i++) {
        for (int j = index; j < index + nbr_chr_with_hap[i]; j++) {
          label_parent[hap[j]] = *nbr_node;
        }
        index += nbr_chr_with_hap[i];
        node_mrk[*nbr_node] = foc_mrk;
        node_parent[*nbr_node] = 0;
        (*nbr_node)++;
      }
    }
    extend_furcation(data, nbr_chr, foc_mrk, end_mrk, lim_haplo, hap, &nbr_hap, nbr_chr_with_hap, node_mrk,
                     node_parent, nbr_node, label_parent);
    
    int tot_nbr_chr_in_hap = 0;
    for (int i = 0; i < nbr_hap; i++) {
      tot_nbr_chr_in_hap += nbr_chr_with_hap[i];
    }
    
    for (int i = 0; i < nbr_chr; i++) {
      if (label_parent[i] != NA_INTEGER) {
        int found = 0;
        for (int j = 0; j < tot_nbr_chr_in_hap; j++) {
          if (hap[j] == i) {
            found = 1;
            break;
          }
        }
        if (!found) {
          node_with_missing_data[label_parent[i]] = 1;
        }
      }
    }
  } else {
    *nbr_node = 0;
  }
  free(hap);
  free(nbr_chr_with_hap);
}
