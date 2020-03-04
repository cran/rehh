#include <R.h>
#include <Rinternals.h>
#include "calc_ehhs.h"

/**
 * Interface between R and C.
 * R objects are marked by a trailing underscore.
 */
SEXP CALL_EHHS(SEXP data_, SEXP nbr_chr_, SEXP nbr_mrk_, SEXP foc_mrk_, SEXP lim_haplo_, SEXP lim_homo_haplo_,
               SEXP lim_ehhs_, SEXP phased_) {

	//get pointer to R data vector
	int* data = INTEGER(data_);

	//translate R vectors of size 1 to numbers
	int nbr_chr = asInteger(nbr_chr_);
	int nbr_mrk = asInteger(nbr_mrk_);
	int foc_mrk = asInteger(foc_mrk_) - 1; //change to C indexing!
	int lim_haplo = asInteger(lim_haplo_);
	int lim_homo_haplo = asInteger(lim_homo_haplo_);
	double lim_ehhs = asReal(lim_ehhs_);
	int phased = asInteger(phased_);

	//create R vectors
	SEXP nhaplo_eval_ = PROTECT(allocVector(INTSXP, nbr_mrk));
	SEXP ehhs_ = PROTECT(allocVector(REALSXP, nbr_mrk));
	SEXP nehhs_ = PROTECT(allocVector(REALSXP, nbr_mrk));

	//get C references for R vectors
	int* nhaplo_eval = INTEGER(nhaplo_eval_);
	double* ehhs = REAL(ehhs_);
	double* nehhs = REAL(nehhs_);

	//perform calculation
	calc_ehhs(data, nbr_chr, nbr_mrk, foc_mrk, lim_haplo, lim_homo_haplo, lim_ehhs, phased, nhaplo_eval, ehhs, nehhs);

	//create R list of length 3
	SEXP list_ = PROTECT(allocVector(VECSXP, 3));

	//add R vectors to list
	SET_VECTOR_ELT(list_, 0, nhaplo_eval_);
	SET_VECTOR_ELT(list_, 1, ehhs_);
	SET_VECTOR_ELT(list_, 2, nehhs_);

	//unprotect all created R objects
	UNPROTECT(4);

	//return R list
	return list_;
}
