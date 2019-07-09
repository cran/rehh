#ifdef Conly
#define NA_INTEGER -1
#else
#include <Rinternals.h>
#endif

void calc_furcation(const int* const data, const int nbr_chr, const int foc_mrk, const int end_mrk, const int foc_all,
		const int lim_haplo, const int phased, int* const node_mrk, int* const node_parent,
		int* const node_with_missing_data, int* const nbr_node, int* const label_parent);
