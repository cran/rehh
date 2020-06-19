#include "definitions.h"

double integrate(const double *x_axis, const double *y_axis, const int n, const int mrk, const double threshold,
		const int scale_gap, const int max_gap, const bool discard_integration_at_border, const double lower_y_bound,
		const bool interpolate);
