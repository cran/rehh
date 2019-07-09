#include <R.h>
#include <Rinternals.h>
#include "calc_furcation.h"

SEXP CALL_FURCATION(SEXP data_, SEXP nbr_chr_, SEXP foc_mrk_, SEXP end_mrk_, SEXP foc_all_, SEXP lim_haplo_,
		SEXP phased_) {
	//get pointer to R data vector
	int *data = INTEGER(data_);

	//translate R vectors of size 1 to numbers
	int nbr_chr = asInteger(nbr_chr_);
	int end_mrk = asInteger(end_mrk_) - 1; //change to C indexing
	int foc_mrk = asInteger(foc_mrk_) - 1; //change to C indexing
	int foc_all = asInteger(foc_all_);
	int lim_haplo = asInteger(lim_haplo_);
	int phased = asInteger(phased_);

	//create R vectors
	SEXP node_mrk_ = PROTECT(allocVector(INTSXP, 2 * nbr_chr - 1));
	SEXP node_parent_ = PROTECT(allocVector(INTSXP, 2 * nbr_chr - 1));
	SEXP node_with_missing_data_ = PROTECT(allocVector(INTSXP, 2 * nbr_chr - 1));
	SEXP label_parent_ = PROTECT(allocVector(INTSXP, nbr_chr));

	//get C references for R vectors
	int* node_mrk = INTEGER(node_mrk_);
	int* node_parent = INTEGER(node_parent_);
	int* node_with_missing_data = INTEGER(node_with_missing_data_);
	int* label_parent = INTEGER(label_parent_);

	int nbr_nodes = 0;

	calc_furcation(data, nbr_chr, foc_mrk, end_mrk, foc_all, lim_haplo, phased, node_mrk, node_parent,
			node_with_missing_data, &nbr_nodes, label_parent);

	//reset length from maximal length to actual length
	SETLENGTH(node_mrk_, nbr_nodes);
	SETLENGTH(node_parent_, nbr_nodes);
	SETLENGTH(node_with_missing_data_, nbr_nodes);

	//create R list of length 8
	SEXP list_ = PROTECT(allocVector(VECSXP, 4));

	//add R vectors to list
	SET_VECTOR_ELT(list_, 0, node_mrk_);
	SET_VECTOR_ELT(list_, 1, node_parent_);
	SET_VECTOR_ELT(list_, 2, node_with_missing_data_);
	SET_VECTOR_ELT(list_, 3, label_parent_);

	//unprotect all created R objects
	UNPROTECT(5);

	//return R list
	return list_;
}
