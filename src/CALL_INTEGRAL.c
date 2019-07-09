#include <R.h>
#include <Rinternals.h>
#include "integrate.h"

/**
 * Interface between R and C.
 * R objects are marked by a trailing underscore.
 */
SEXP CALL_INTEGRAL(SEXP map_, SEXP ehh_, SEXP mrk_, SEXP lim_ehh_, SEXP scale_gap_, SEXP max_gap_,
		SEXP discard_integration_at_border_, SEXP lower_y_bound_) {
	int nbr_mrk = length(ehh_);

	//get pointer to R data vectors
	double* map = REAL(map_);
	double* ehh = REAL(ehh_);

	//extract numbers from R objects
	int mrk = asInteger(mrk_) - 1;
	double lim_ehh = asReal(lim_ehh_);
	int scale_gap = asInteger(scale_gap_);
	int max_gap = asInteger(max_gap_);
	int discard_integration_at_border = asInteger(discard_integration_at_border_);
	double lower_y_bound = asReal(lower_y_bound_);

	double integral = integrate(map, ehh, nbr_mrk, mrk, lim_ehh, scale_gap, max_gap, discard_integration_at_border,
			lower_y_bound);

	return (ScalarReal(integral));
}
