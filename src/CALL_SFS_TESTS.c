#include <R.h>
#include "definitions.h"
#include "calc_sfs_tests.h"

/**
 * Interface between R and C.
 * R objects are marked by a trailing underscore.
 */
SEXP CALL_SFS_TESTS(SEXP data_, SEXP nbr_chr_, SEXP nbr_mrk_, SEXP map_,
                    SEXP polarized_, SEXP windows_, SEXP n_windows_, SEXP right_, 
                    SEXP min_n_mrk_, SEXP n_mrk_, SEXP results_) {
  
  //get pointer to R data vectors
  int* data = INTEGER(data_);
  double* windows = REAL(windows_);
  double* map = REAL(map_);
  int* n_mrk = INTEGER(n_mrk_);
  double *results = REAL(results_);
  
  //translate R vectors of size 1 to numbers
  int nbr_chr = asInteger(nbr_chr_);
  int nbr_mrk = asInteger(nbr_mrk_);
  bool polarized = asLogical(polarized_);
  bool right = asLogical(right_);
  int n_windows = asInteger(n_windows_);
  int min_n_mrk = asInteger(min_n_mrk_);
  
  calc_sfs_tests(data, nbr_chr, nbr_mrk, map, polarized,
                 windows, n_windows, right, min_n_mrk, 
                 n_mrk, results);
  
  return ScalarLogical(1);
}
