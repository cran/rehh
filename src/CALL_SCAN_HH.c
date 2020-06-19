#include <R.h>
#include "definitions.h"
#include "calc_ehh.h"
#include "calc_ehhs.h"
#include "integrate.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * Interface between R and C.
 * R objects are marked by a trailing underscore.
 */
SEXP CALL_SCAN_HH(SEXP data_, SEXP nbr_chr_, SEXP nbr_mrk_, SEXP first_allele_, SEXP second_allele_, SEXP map_,
                  SEXP lim_haplo_, SEXP lim_homo_haplo_, SEXP lim_ehh_, SEXP lim_ehhs_, 
                  SEXP lower_ehh_y_bound_, SEXP lower_ehhs_y_bound_, SEXP scale_gap_, SEXP max_gap_,  
                  SEXP phased_, SEXP discard_integration_at_border_, SEXP interpolate_, SEXP nbr_threads_) {

	//get pointer to R data vectors
	int* data = INTEGER(data_);
	double *map = REAL(map_);
	int* first_allele = INTEGER(first_allele_);
	int* second_allele = INTEGER(second_allele_);

	//translate R vectors of size 1 to numbers
	int nbr_chr = asInteger(nbr_chr_);
	int nbr_mrk = asInteger(nbr_mrk_);
	int lim_haplo = asInteger(lim_haplo_);
	int lim_homo_haplo = asInteger(lim_homo_haplo_);
	double lim_ehh = asReal(lim_ehh_);
	double lim_ehhs = asReal(lim_ehhs_);
	double lower_ehh_y_bound = asReal(lower_ehh_y_bound_);
	double lower_ehhs_y_bound = asReal(lower_ehhs_y_bound_);
	int max_gap = asInteger(max_gap_);
	int scale_gap = asInteger(scale_gap_);
	bool phased = asLogical(phased_);
	bool discard_integration_at_border = asLogical(discard_integration_at_border_);
	bool interpolate = asLogical(interpolate_);
	int nbr_threads = asInteger(nbr_threads_);

	//create R numerical vectors
	SEXP nhaplo_A_ = PROTECT(allocVector(INTSXP, nbr_mrk));
	SEXP nhaplo_D_ = PROTECT(allocVector(INTSXP, nbr_mrk));
	SEXP ihhA_ = PROTECT(allocVector(REALSXP, nbr_mrk));
	SEXP ihhD_ = PROTECT(allocVector(REALSXP, nbr_mrk));
	SEXP ies_ = PROTECT(allocVector(REALSXP, nbr_mrk));
	SEXP ines_ = PROTECT(allocVector(REALSXP, nbr_mrk));

	//get C references for R vectors
	int* nhaplo_A = INTEGER(nhaplo_A_);
	int* nhaplo_D = INTEGER(nhaplo_D_);
	double* ihhA = REAL(ihhA_);
	double* ihhD = REAL(ihhD_);
	double* ies = REAL(ies_);
	double* ines = REAL(ines_);

	//init output vectors by zeros
	for (int i = 0; i < nbr_mrk; i++) {
	  nhaplo_A[i] = 0;
	  nhaplo_D[i] = 0;
		ihhA[i] = 0.0;
		ihhD[i] = 0.0;
		ies[i] = 0.0;
		ines[i] = 0.0;
	}

	int j;
	
#ifdef _OPENMP
#pragma omp parallel num_threads(nbr_threads) private(j)
{
#endif
  
  int* nhaplo = (int *) malloc(nbr_mrk * sizeof(int));
  double* ehh = (double *) malloc(nbr_mrk * sizeof(double));
  double* ehhs = (double *) malloc(nbr_mrk * sizeof(double));
  double* nehhs = (double *) malloc(nbr_mrk * sizeof(double));
  
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 256)
#endif
  
	for (j = 0; j < nbr_mrk; j++) {

		// compute EHH for the ancestral allele
		calc_ehh(data, nbr_chr, nbr_mrk, j, first_allele[j], lim_haplo, lim_homo_haplo, lim_ehh, phased, nhaplo,
				ehh);
		
		// store nhaplo for ancestral allele at focal marker
		nhaplo_A[j] = nhaplo[j];
		
		// compute IHH for the ancestral allele
		ihhA[j] = integrate(map, ehh, nbr_mrk, j, lim_ehh, scale_gap, max_gap, discard_integration_at_border,
                      lower_ehh_y_bound, interpolate);

		// compute EHH for the derived allele
		calc_ehh(data, nbr_chr, nbr_mrk, j, second_allele[j], lim_haplo, lim_homo_haplo, lim_ehh, phased, nhaplo,
				ehh);
		
		// store nhaplo for derived allele at focal marker
		nhaplo_D[j] = nhaplo[j];
		
		// compute IHH for the derived allele
		ihhD[j] = integrate(map, ehh, nbr_mrk, j, lim_ehh, scale_gap, max_gap, discard_integration_at_border,
                      lower_ehh_y_bound, interpolate);

		//compute EHHS for both Tang et al.'s (2007) and Sabeti et al.'s (2007) definitions
		calc_ehhs(data, nbr_chr, nbr_mrk, j, lim_haplo, lim_homo_haplo, lim_ehhs, phased, nhaplo,
				ehhs, nehhs);
		//compute IES, using Tang et al.'s (2007) definition of EHHS
		ies[j] = integrate(map, ehhs, nbr_mrk, j, lim_ehhs, scale_gap, max_gap,
				discard_integration_at_border, lower_ehhs_y_bound, interpolate);
		//compute IES, using Sabeti et al.'s (2007) definition of EHHS
		ines[j] = integrate(map, nehhs, nbr_mrk, j, lim_ehhs, scale_gap, max_gap,
				discard_integration_at_border, lower_ehhs_y_bound, interpolate);
	}


	free(ehh);
	free(ehhs);
	free(nehhs);
	free(nhaplo);
	
#ifdef _OPENMP
}
#endif
  
	//create R list of length 6
	SEXP list_ = PROTECT(allocVector(VECSXP, 6));

	//add R vectors to list
	SET_VECTOR_ELT(list_, 0, nhaplo_A_);
	SET_VECTOR_ELT(list_, 1, nhaplo_D_);
	SET_VECTOR_ELT(list_, 2, ihhA_);
	SET_VECTOR_ELT(list_, 3, ihhD_);
	SET_VECTOR_ELT(list_, 4, ies_);
	SET_VECTOR_ELT(list_, 5, ines_);

	//unprotect all created R objects
	UNPROTECT(7);

	//return R list
	return list_;
}
