#include <stdio.h>
#include <stdlib.h>

#ifdef Conly
#define MISSING_VALUE -1
#else
#include <Rinternals.h>
#define MISSING_VALUE NA_INTEGER
#endif

void init_allele_hap(const int* const data, const int nbr_chr, const int foc_mrk, const int foc_all, const int phased,
		int* const hap, int* const nbr_hap, int* const hap_count);
void init_site_hap(const int* const data, const int nbr_chr, const int foc_mrk, const int phased, int* const hap,
		int* const nbr_hap, int* const hap_count);
int update_hap(const int* const data, const int nbr_chr, const int mrk, int* const hap, int* const nbr_hap,
		int* const hap_count);
