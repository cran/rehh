#include <R.h>
#include "definitions.h"
#include "calc_pairwise_haplen.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * Interface between R and C.
 * R objects are marked by a trailing underscore.
 */
SEXP CALL_SCAN_HH_FULL(SEXP data_, SEXP nbr_chr_, SEXP nbr_mrk_, SEXP first_allele_, SEXP second_allele_,
                   SEXP map_, SEXP max_gap_, SEXP max_extend_, SEXP phased_, SEXP discard_integration_at_border_, 
                   SEXP geometric_mean_, SEXP nbr_threads_) {
  
  //get pointer to R data vectors
  int* data = INTEGER(data_);
  double *map = REAL(map_);
  int* first_allele = INTEGER(first_allele_);
  int* second_allele = INTEGER(second_allele_);
  
  //translate R vectors of size 1 to numbers
  int nbr_chr = asInteger(nbr_chr_);
  int nbr_mrk = asInteger(nbr_mrk_);
  int max_gap = asInteger(max_gap_);
  int max_extend = asInteger(max_extend_);
  bool phased = asLogical(phased_);
  bool discard_integration_at_border = asLogical(discard_integration_at_border_);
  bool geometric_mean = asLogical(geometric_mean_);
  int nbr_threads = asInteger(nbr_threads_);
  
  //create R numerical vectors
  SEXP nhaplo_A_ = PROTECT(allocVector(INTSXP, nbr_mrk));
  SEXP nhaplo_D_ = PROTECT(allocVector(INTSXP, nbr_mrk));
  SEXP ihh_A_ = PROTECT(allocVector(REALSXP, nbr_mrk));
  SEXP ihh_D_ = PROTECT(allocVector(REALSXP, nbr_mrk));
  SEXP ies_ = PROTECT(allocVector(REALSXP, nbr_mrk));
  SEXP ines_ = PROTECT(allocVector(REALSXP, nbr_mrk));
  
  //get C references for R vectors
  int* nhaplo_A = INTEGER(nhaplo_A_);
  int* nhaplo_D = INTEGER(nhaplo_D_);
  double* ihh_A = REAL(ihh_A_);
  double* ihh_D = REAL(ihh_D_);
  double* ies = REAL(ies_);
  double* ines = REAL(ines_);
  
  //init output vectors by zeros
  for (int i = 0; i < nbr_mrk; i++) {
    nhaplo_A[i] = 0;
    nhaplo_D[i] = 0;
    ihh_A[i] = 0.0;
    ihh_D[i] = 0.0;
    ies[i] = 0.0;
    ines[i] = 0.0;
  }
  
  int j;
  
#ifdef _OPENMP
#pragma omp parallel num_threads(nbr_threads) private(j)
{
#endif
  
  double *pairwise_haplen = (double*) malloc(nbr_chr * nbr_chr * sizeof(double));
  
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 256)
#endif
  
  for (j = 0; j < nbr_mrk; j++) {
    for(int i = 0; i < nbr_chr * nbr_chr; i++){
      pairwise_haplen[i] = 0;
    }
    
    bool discard = calc_pairwise_haplen(data, nbr_chr, nbr_mrk, map, j, max_gap, max_extend, BOTH,
                                        phased, discard_integration_at_border, pairwise_haplen) == 1;
    
    if(discard){
      int nbr_first = 0;
      int nbr_second = 0;
      
      if(phased){
        for(int k = 0; k < nbr_chr; k++){
          int foc_allele1 = data[j * nbr_chr + k];
          if(foc_allele1 == first_allele[j]){
            nbr_first++;
          } else if(foc_allele1 == second_allele[j]){
            nbr_second++;
          }
        }
      }else{
        for(int k = 1; k < nbr_chr; k += 2){
          int foc_allele1 = data[j * nbr_chr + k];
          if(data[j * nbr_chr + k - 1] == foc_allele1){
            if(foc_allele1 == first_allele[j]){
              nbr_first += 2;
            }else if(foc_allele1 == second_allele[j]){
              nbr_second += 2;
            }
          }
        }
      }
      
      nhaplo_A[j] = nbr_first;
      nhaplo_D[j] = nbr_second;
      ihh_A[j] = NA_REAL;
      ihh_D[j] = NA_REAL;
      ies[j] = NA_REAL;
      ines[j] = NA_REAL;
      
    }else{
      
      if(geometric_mean){
        for(int k = 1; k < nbr_chr; k++){
          for(int l = 0; l < k; l++){
            if(pairwise_haplen[k * nbr_chr + l] > 0){
              pairwise_haplen[k * nbr_chr + l] = log(pairwise_haplen[k * nbr_chr + l]);
            }
          }
        }
      }
      int nbr_first = 0;
      int nbr_second = 0;
      int nbr_other = 0;
      double sum_first = 0.0;
      double sum_second = 0.0;
      double sum_other = 0.0;
      
      if(phased){
        //count alleles and sum mutual shared haplotype lengths
        for(int k = 0; k < nbr_chr; k++){
          int foc_allele1 = data[j * nbr_chr + k];
          if(foc_allele1 == first_allele[j]){
            nbr_first++;
            for(int l = 0; l < k; l++){
              if(data[j * nbr_chr + l] == foc_allele1){
                sum_first += pairwise_haplen[k * nbr_chr + l];
              }
            }
          } else if(foc_allele1 == second_allele[j]){
            nbr_second++;
            for(int l = 0; l < k; l++){
              if(data[j * nbr_chr + l] == foc_allele1){
                sum_second += pairwise_haplen[k * nbr_chr + l];
              }
            }
          } else{
            nbr_other++;
            for(int l = 0; l < k; l++){
              if(data[j * nbr_chr + l] == foc_allele1){
                sum_other += pairwise_haplen[k * nbr_chr + l];
              }
            }
          }
        }
      }else{
        for(int k = 1; k < nbr_chr; k += 2){
          int foc_allele1 = data[j * nbr_chr + k];
          if(data[j * nbr_chr + k - 1] == foc_allele1){
            double length = pairwise_haplen[k * nbr_chr + k - 1];
            if(foc_allele1 == first_allele[j]){
              nbr_first += 2;
              sum_first += length;
            }else if(foc_allele1 == second_allele[j]){
              nbr_second += 2;
              sum_second += length;
            } else {
              nbr_other += 2;
              sum_other += length;
            }
          }
        }
      }
      
      if(nbr_first > 0){
        nhaplo_A[j] = nbr_first;
        if(sum_first > 0){  // implies nhaplo_A > 1
          ihh_A[j] = sum_first * 2 / nbr_first;
          if(phased){
            ihh_A[j] /= nbr_first - 1;
          }
        }
      }
      if(nbr_second > 0){
        nhaplo_D[j] = nbr_second;
        if(sum_second > 0){ // implies nhaplo_D > 1
          ihh_D[j] = sum_second * 2 / nbr_second;
          if(phased){
            ihh_D[j] /= nbr_second - 1;
          }
        }
      }
      
      double total_sum = sum_first + sum_second + sum_other;
      
      int total_nbr_comparisons = phased ? nbr_first * (nbr_first-1) + 
        nbr_second * (nbr_second - 1) + nbr_other * (nbr_other - 1) : nbr_first + nbr_second + nbr_other;
      if(total_sum > 0){
        ines[j] = total_sum * 2 / total_nbr_comparisons;
      }
      
      if(phased){
        total_nbr_comparisons = (nbr_first + nbr_second + nbr_other) * 
          (nbr_first + nbr_second + nbr_other - 1);
        if(total_sum > 0){
          ies[j] = total_sum * 2 / total_nbr_comparisons;
        }
      }else{
        ies[j] = ines[j];
      }
      
      if(geometric_mean){
        if(ihh_A[j] != 0){
          ihh_A[j] = exp(ihh_A[j]);
        }
        if(ihh_D[j] != 0){
          ihh_D[j] = exp(ihh_D[j]);
        }
        ies[j] = NA_REAL;
        if(ines[j] != 0){
          ines[j] = exp(ines[j]);
        }
      }
    }
  }
  
  free(pairwise_haplen);
  
#ifdef _OPENMP
}
#endif

//create R list of length 6
SEXP list_ = PROTECT(allocVector(VECSXP, 6));

//add R vectors to list
SET_VECTOR_ELT(list_, 0, nhaplo_A_);
SET_VECTOR_ELT(list_, 1, nhaplo_D_);
SET_VECTOR_ELT(list_, 2, ihh_A_);
SET_VECTOR_ELT(list_, 3, ihh_D_);
SET_VECTOR_ELT(list_, 4, ies_);
SET_VECTOR_ELT(list_, 5, ines_);

//unprotect all created R objects
UNPROTECT(7);

//return R list
return list_;
}
