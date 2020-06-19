#include "definitions.h"

int calc_pairwise_haplen(const int* const data, const int nbr_chr, const int nbr_mrk, const double* map,
                         const int foc_mrk, const int max_gap, 
                         const bool phased, const bool discard_integration_at_border, double* const haplength);
