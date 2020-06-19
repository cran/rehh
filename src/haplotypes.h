#include "definitions.h"

void init_allele_hap(const int* const data, const int nbr_chr, const int foc_mrk, const int foc_all, const bool phased,
                     int* const hap, int* const nbr_hap, int* const hap_count);
void init_site_hap(const int* const data, const int nbr_chr, const int foc_mrk, const bool phased, int* const hap,
                   int* const nbr_hap, int* const hap_count);
int update_hap(const int* const data, const int nbr_chr, const int mrk, int* const hap, int* const nbr_hap,
               int* const hap_count);
