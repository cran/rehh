void init_allele_hap_with_nodes(const int *data, const int nbr_chr, const int foc_mrk, const int foc_all, const int phased,
                                int* const hap, int* const nbr_hap, int* const nbr_chr_with_hap,
                                int* const hap_node, int* const node_size, 
                                int* const node_mrk, int* const node_parent, int* const nbr_node);

int update_hap_with_nodes(const int *data, const int nbr_chr, const int mrk, int* const hap, int* const nbr_hap,
                          int* const nbr_chr_with_hap, int* const hap_node, int* const node_size,  int* const node_mrk, 
                          int* const node_parent, int* const node_with_missing_data, int* const nbr_node, int* const label_parent);
