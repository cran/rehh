#include "definitions.h"

void calc_ehhs(const int* const data, const int nbr_chr, const int nbr_mrk, const int foc_mrk, const int lim_haplo,
		const int lim_homo_haplo, const double lim_ehhs, const bool phased, int* const tot_nbr_chr_in_hap, double* const ehhs,
		double* const nehhs);
