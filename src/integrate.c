#include "integrate.h"

double integrate(const double *x, const double *y, const int n, const int mrk, const double threshold,
		const int scale_gap, const int max_gap, const int discard_integration_at_border, const double lower_y_bound) {
	double area = 0.0;
	double dx;

	if (discard_integration_at_border && ((y[0] > threshold) || (y[n - 1] > threshold))) { // If the EHH or EHHS is larger than the minimum value at either end of the chromosome, ...
		return (NA_REAL);                                              // ... then do not compute the integral, and quit
	}

	if (y[mrk] <= threshold) {
		return (area);
	}

	for (int i = mrk; i > 0; i--) {
		dx = x[i] - x[i - 1];
		if (dx > max_gap) {   // If a gap larger than max_gap exists within the 'support' of the EHH or EHHS, ...
			if (discard_integration_at_border) {
				return (NA_REAL); // ... then do not compute the integral, and quit
			} else {
				break;            // ... stop integration
			}
		}
		if (dx > scale_gap) {
			dx = scale_gap;
		}
		if (y[i - 1] > threshold) {
			area += dx * ((y[i] + y[i - 1]) / 2 - lower_y_bound);
		} else {
			area += dx * (y[i] - lower_y_bound) * (y[i] - lower_y_bound) / (2 * y[i]);
			break;
		}
	}
	for (int i = mrk; i < n - 1; i++) {
		dx = x[i + 1] - x[i];
		if (dx > max_gap) {   // If a gap larger than max_gap exists within the 'support' of the EHH or EHHS, ...
			if (discard_integration_at_border) {
				return (NA_REAL); // ... then do not compute the integral, and quit
			} else {
				break;            // ... stop integration
			}
		}
		if (dx > scale_gap) {
			dx = scale_gap;
		}
		if (y[i + 1] > threshold) {
			area += dx * ((y[i + 1] + y[i]) / 2 - lower_y_bound);
		} else {
			area += dx * (y[i] - lower_y_bound) * (y[i] - lower_y_bound) / (2 * y[i]);
			break;
		}
	}
	return (area);
}
