#include <R.h>
#include "hh_utils.h"

void CALL_EHHS(int *Rdata,
               int *focl_snp,
               int *nbr_snps,
               int *nbr_chrs,
               int *tot_nbr_hapl,
               int *min_nbr_hapl,
               double *min_ehhs,
               double *max_gap,
               double *map,
               double *ehhs_tang,
               double *ies_tang,
               double *ehhs_sabeti,
               double *ies_sabeti)

{
  int i,j;
  short int **data;
  
  *focl_snp -= 1;                                                               // This is to account for C numbering of arrays
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
  for (j = 0; j < *nbr_snps; j++) {
    ehhs_tang[j] = 0.0;                                                         // Initialize the `ehhs_tang' vector, for Tang et al.'s (2007) definition
    ehhs_sabeti[j] = 0.0;                                                       // Initialize the `ehhs_sabeti' vector, for Sabeti et al.'s (2007) definition
    tot_nbr_hapl[j] = 0;                                                        // Initialize the total number of haplotypes [N.B: might be unnecessary, tot_nbr_hapl[focl_snp] = 0; should be sufficient...]
  }
  compute_ehhs(data,map,*focl_snp,*nbr_snps,*nbr_chrs,tot_nbr_hapl,*min_nbr_hapl,*min_ehhs,*max_gap,ehhs_tang,ehhs_sabeti); // compute the EHHS for both Tang et al.'s (2007) and Sabeti et al.'s (2007) definitions
  *ies_tang = integrate(map,ehhs_tang,*nbr_snps,*min_ehhs);                     // compute the IES, for Tang et al.'s (2007) definition
  *ies_sabeti = integrate(map,ehhs_sabeti,*nbr_snps,*min_ehhs);                 // compute the IES, for Sabeti et al.'s (2007) definition
  for (i = 0; i < *nbr_chrs; i++) {                                             // Free the memory for the `data' array
    free(data[i]);
  }
  free(data);
}
