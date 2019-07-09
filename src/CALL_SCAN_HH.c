#include <R.h>
#include <Rinternals.h>
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
SEXP CALL_SCAN_HH(SEXP data_, SEXP nbr_chr_, SEXP nbr_mrk_, SEXP first_allele_, SEXP second_allele_,
                  SEXP lim_haplo_, SEXP lim_ehh_, SEXP lim_ehhs_,
                  SEXP scale_gap_, SEXP max_gap_, SEXP map_, SEXP phased_, SEXP discard_integration_at_border_,
                  SEXP lower_ehh_y_bound_, SEXP lower_ehhs_y_bound_, SEXP nbr_threads_) {

	//get pointer to R data vectors
	int* data = INTEGER(data_);
	double *map = REAL(map_);
	int* first_allele = INTEGER(first_allele_);
	int* second_allele = INTEGER(second_allele_);

	//translate R vectors of size 1 to numbers
	int nbr_chr = asInteger(nbr_chr_);
	int nbr_mrk = asInteger(nbr_mrk_);
	int lim_haplo = asInteger(lim_haplo_);
	double lim_ehh = asReal(lim_ehh_);
	double lim_ehhs = asReal(lim_ehhs_);
	int max_gap = asInteger(max_gap_);
	int scale_gap = asInteger(scale_gap_);
	int phased = asInteger(phased_);
	int discard_integration_at_border = asInteger(discard_integration_at_border_);
	double lower_ehh_y_bound = asReal(lower_ehh_y_bound_);
	double lower_ehhs_y_bound = asReal(lower_ehhs_y_bound_);
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

	int idx_thread;
	int **nhaplo; // re-used for each marker
	double **ehh, **ehhs, **nehhs;

	nhaplo = (int **) malloc(nbr_threads * sizeof(int *));          // Allocate memory for thread-specific vectors
	ehh = (double **) malloc(nbr_threads * sizeof(double *));
	ehhs = (double **) malloc(nbr_threads * sizeof(double *));
	nehhs = (double **) malloc(nbr_threads * sizeof(double *));

	for (int i = 0; i < nbr_threads; i++) {
		nhaplo[i] = (int *) malloc(nbr_mrk * sizeof(int));
		ehh[i] = (double *) malloc(nbr_mrk * sizeof(double));
		ehhs[i] = (double *) malloc(nbr_mrk * sizeof(double));
		nehhs[i] = (double *) malloc(nbr_mrk * sizeof(double));
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(nbr_threads) private(idx_thread)
#endif

	for (int j = 0; j < nbr_mrk; j++) {

#ifdef _OPENMP
		idx_thread = omp_get_thread_num();
#else
		idx_thread = 0;
#endif
		// compute EHH for the ancestral allele
		calc_ehh(data, nbr_chr, nbr_mrk, j, first_allele[j], lim_haplo, lim_ehh, phased, nhaplo[idx_thread],
				ehh[idx_thread]);
		
		// store nhaplo for ancestral allele at focal marker
		nhaplo_A[j] = nhaplo[idx_thread][j];
		
		// compute IHH for the ancestral allele
		ihhA[j] = integrate(map, ehh[idx_thread], nbr_mrk, j, lim_ehh, scale_gap, max_gap, discard_integration_at_border,
                      lower_ehh_y_bound);

		// compute EHH for the derived allele
		calc_ehh(data, nbr_chr, nbr_mrk, j, second_allele[j], lim_haplo, lim_ehh, phased, nhaplo[idx_thread],
				ehh[idx_thread]);
		
		// store nhaplo for derived allele at focal marker
		nhaplo_D[j] = nhaplo[idx_thread][j];
		
		// compute IHH for the derived allele
		ihhD[j] = integrate(map, ehh[idx_thread], nbr_mrk, j, lim_ehh, scale_gap, max_gap, discard_integration_at_border,
                      lower_ehh_y_bound);

		//compute EHHS for both Tang et al.'s (2007) and Sabeti et al.'s (2007) definitions
		calc_ehhs(data, nbr_chr, nbr_mrk, j, lim_haplo, lim_ehhs, phased, nhaplo[idx_thread],
				ehhs[idx_thread], nehhs[idx_thread]);
		//compute IES, using Tang et al.'s (2007) definition of EHHS
		ies[j] = integrate(map, ehhs[idx_thread], nbr_mrk, j, lim_ehhs, scale_gap, max_gap,
				discard_integration_at_border, lower_ehhs_y_bound);
		//compute IES, using Sabeti et al.'s (2007) definition of EHHS
		ines[j] = integrate(map, nehhs[idx_thread], nbr_mrk, j, lim_ehhs, scale_gap, max_gap,
				discard_integration_at_border, lower_ehhs_y_bound);
	}

	for (int i = 0; i < nbr_threads; i++) {
		free(ehh[i]);
		free(ehhs[i]);
		free(nehhs[i]);
		free(nhaplo[i]);
	}
	free(ehh);
	free(ehhs);
	free(nehhs);
	free(nhaplo);

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
