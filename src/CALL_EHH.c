#include <R.h>
#include <Rinternals.h>
#include "calc_ehh.h"

/**
 * Interface between R and C.
 * R objects are marked by a trailing underscore.
 */
SEXP CALL_EHH(SEXP data_, SEXP nbr_chr_, SEXP nbr_mrk_, SEXP foc_mrk_, SEXP foc_all_, SEXP lim_haplo_, SEXP lim_homo_haplo_,
              SEXP lim_ehh_, SEXP phased_) {

	//get pointer to R data vector
	int *data = INTEGER(data_);

	//translate R vectors of size 1 to numbers
	int nbr_chr = asInteger(nbr_chr_);
	int nbr_mrk = asInteger(nbr_mrk_);
	int foc_mrk = asInteger(foc_mrk_) - 1; //change to C indexing!
	int foc_all = asInteger(foc_all_);
	int lim_haplo = asInteger(lim_haplo_);
	int lim_homo_haplo = asInteger(lim_homo_haplo_);
	double lim_ehh = asReal(lim_ehh_);
	int phased = asInteger(phased_);

	//create R integer vector
	SEXP nhaplo_eval_ = PROTECT(allocVector(INTSXP, nbr_mrk));
	//create R numeric vector
	SEXP ehh_ = PROTECT(allocVector(REALSXP, nbr_mrk));

	//get C references for R vectors
	int* nhaplo_eval = INTEGER(nhaplo_eval_);
	double* ehh = REAL(ehh_);

	//perform calculation
	calc_ehh(data, nbr_chr, nbr_mrk, foc_mrk, foc_all, lim_haplo, lim_homo_haplo, lim_ehh, phased, nhaplo_eval, ehh);

	//create R list of length 2
	SEXP list_ = PROTECT(allocVector(VECSXP, 2));

	//add R vectors to list
	SET_VECTOR_ELT(list_, 0, nhaplo_eval_);
	SET_VECTOR_ELT(list_, 1, ehh_);

	//unprotect all three created R objects
	UNPROTECT(3);

	//return R list
	return list_;
}
