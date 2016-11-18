#include <R.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "hh_utils.h"

void CALL_SCAN_HH(int *Rdata,
                  int *nbr_snps,
                  int *nbr_chrs,
                  int *min_nbr_hapl,
                  double *min_ehh,
                  double *min_ehhs,
                  double *max_gap,
                  double *map,
                  double *ihh_vect,
                  double *ies_tang_vect,
                  double *ies_sabeti_vect,
                  int *nbr_threads)

{
  int i,j;
  int idx_thread;
  short int **data;
  int **tot_nbr_hapl;
  double **ehh,**ehhs_tang,**ehhs_sabeti;
  
  tot_nbr_hapl = (int **) malloc(*nbr_threads * sizeof(int *));                 // Allocate memory for thread-specific vectors
  ehh = (double **) malloc(*nbr_threads * sizeof(double *));
  ehhs_tang = (double **) malloc(*nbr_threads * sizeof(double *));
  ehhs_sabeti = (double **) malloc(*nbr_threads * sizeof(double *));
  for (i = 0; i < *nbr_threads; i++) {
    tot_nbr_hapl[i] = (int *) malloc(*nbr_snps * sizeof(int));
    ehh[i] = (double *) malloc(*nbr_snps * sizeof(double));
    ehhs_tang[i] = (double *) malloc(*nbr_snps * sizeof(double));
    ehhs_sabeti[i] = (double *) malloc(*nbr_snps * sizeof(double));
  }
  data = (short int **) malloc(*nbr_chrs * sizeof(short int *));                // Create a `data' array (nbr_chrs x nbr_snps) for C code
  for (i = 0; i < *nbr_chrs; i++) {
    data[i] = (short int *) malloc(*nbr_snps * sizeof(short int));
  }
  for (i = 0; i < *nbr_chrs; i++) {                                             // Copy information from the `Rdata' vector into the `data' array
    for (j = 0; j < *nbr_snps; j++) {
      switch (Rdata[(j * *nbr_chrs) + i]) {
          
        case 0:
          data[i][j] = MISSING;                                                 // missing data is encoded as `0' in the `Rdata' vector, and as MISSING (= `9') in the `data' array
          break;
          
        case 1:
          data[i][j] = ANCSTRL;                                                 // ancestral allele is encoded as `1' in the `Rdata' vector, and as ANCSTRL (= `0') in the `data' array
          break;
          
        case 2:
          data[i][j] = DERIVED;                                                 // derived allele is encoded as `2' in the `Rdata' vector, and as DERIVED (= `1') in the `data' array
          break;
          
      }
    }
  }
#ifdef _OPENMP
#pragma omp parallel for num_threads(*nbr_threads) private(i,idx_thread)
#endif
  for (j = 0; j < *nbr_snps; j++) {
#ifdef _OPENMP
    idx_thread = omp_get_thread_num();
#else
    idx_thread = 0;
#endif
    for (i = 0; i < *nbr_snps; i++) {
      ehh[idx_thread][i] = 0.0;                                                 // Initialize the `ehh' vector
      tot_nbr_hapl[idx_thread][i] = 0;                                          // Initialize the total number of haplotypes [N.B: might be unnecessary, tot_nbr_hapl[j] = 0; should be sufficient...]
    }
    compute_ehh(data,j,ANCSTRL,*nbr_snps,*nbr_chrs,tot_nbr_hapl[idx_thread],*min_nbr_hapl,*min_ehh,ehh[idx_thread]); // compute the EHH for the ancestral allele
    ihh_vect[j] = integrate(map,ehh[idx_thread],*nbr_snps,*min_ehh,*max_gap);   // compute the IHH for the ancestral allele
    for (i = 0; i < *nbr_snps; i++) {
      ehh[idx_thread][i] = 0.0;                                                 // Initialize the `ehh' vector
      tot_nbr_hapl[idx_thread][i] = 0;                                          // Initialize the total number of haplotypes [N.B: might be unnecessary, tot_nbr_hapl[j] = 0; should be sufficient...]
    }
    compute_ehh(data,j,DERIVED,*nbr_snps,*nbr_chrs,tot_nbr_hapl[idx_thread],*min_nbr_hapl,*min_ehh,ehh[idx_thread]); // compute the EHH for the derived allele
    ihh_vect[j + *nbr_snps] = integrate(map,ehh[idx_thread],*nbr_snps,*min_ehh,*max_gap); // compute the IHH for the derived allele
    for (i = 0; i < *nbr_snps; i++) {
      ehhs_tang[idx_thread][i] = 0.0;                                           // Initialize the `ehhs' vector
      ehhs_sabeti[idx_thread][i] = 0.0;                                         // Initialize the `ehhs' vector
      tot_nbr_hapl[idx_thread][i] = 0;                                          // Initialize the total number of haplotypes [N.B: might be unnecessary, tot_nbr_hapl[j] = 0; should be sufficient...]
    }
    compute_ehhs(data,j,*nbr_snps,*nbr_chrs,tot_nbr_hapl[idx_thread],*min_nbr_hapl,*min_ehhs,ehhs_tang[idx_thread],ehhs_sabeti[idx_thread]); // compute the EHHS for both Tang et al.'s (2007) and Sabeti et al.'s (2007) definitions
    ies_tang_vect[j] = integrate(map,ehhs_tang[idx_thread],*nbr_snps,*min_ehhs,*max_gap); // compute the IES, using Tang et al.'s (2007) definition
    ies_sabeti_vect[j] = integrate(map,ehhs_sabeti[idx_thread],*nbr_snps,*min_ehhs,*max_gap); // compute the IES, using Sabeti et al.'s (2007) definition
  }
  for (i = 0; i < *nbr_chrs; i++) {                                             // Free the memory for the `data' array
    free(data[i]);
  }
  free(data);
  for (i = 0; i < *nbr_threads; i++) {
    free(ehh[i]);
    free(ehhs_tang[i]);
    free(ehhs_sabeti[i]);
    free(tot_nbr_hapl[i]);
  }
  free(ehh);
  free(ehhs_tang);
  free(ehhs_sabeti);
  free(tot_nbr_hapl);
}
