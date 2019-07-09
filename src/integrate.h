#ifdef Conly
#define NA_REAL -1
#else
#include <Rinternals.h>
#endif

double integrate(const double *x_axis, const double *y_axis, const int n, const int mrk, const double threshold,
		const int scale_gap, const int max_gap, const int discard_integration_at_border, const double lower_y_bound);
