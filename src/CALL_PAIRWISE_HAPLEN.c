#include <R.h>
#include "definitions.h"
#include "calc_pairwise_haplen.h"

/**
 * Interface between R and C.
 * R objects are marked by a trailing underscore.
 */
SEXP CALL_PAIRWISE_HAPLEN(SEXP data_, SEXP nbr_chr_, SEXP nbr_mrk_, SEXP map_,
                          SEXP foc_mrk_, SEXP max_gap_, SEXP phased_, SEXP pairwise_haplen_) {

	//get pointer to R data vector
	int* data = INTEGER(data_);

	//translate R vectors of size 1 to numbers
	int nbr_chr = asInteger(nbr_chr_);
	int nbr_mrk = asInteger(nbr_mrk_);
	int foc_mrk = asInteger(foc_mrk_) - 1; //change to C indexing!
	double max_gap = asReal(max_gap_);
	double* map = REAL(map_);
	bool phased = asLogical(phased_);
	
	double* pairwise_haplen = REAL(pairwise_haplen_);

	//perform calculation
	calc_pairwise_haplen(data, nbr_chr, nbr_mrk, map, foc_mrk, max_gap, phased, 0, pairwise_haplen);
	
	return ScalarLogical(1);
}
