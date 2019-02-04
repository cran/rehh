#include <R.h>
#include "hh_utils.h"

void CALL_EHH(int *Rdata,
              int *focl_snp,
              int *nbr_snps,
              int *nbr_chrs,
              int *tot_nbr_hapl,
              int *min_nbr_hapl,
              int *discard_integration_at_border,
              double *min_ehh,
              double *max_gap,
              double *map,
              double *ehh,
              double *ihh)

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
  for (j = 0; j < 2 * *nbr_snps; j++) {
    ehh[j] = 0.0;                                                               // Initialize the `ehh' vector
    tot_nbr_hapl[j] = 0;                                                        // Initialize the total number of haplotypes [N.B: might be unnecessary, tot_nbr_hapl[focl_snp] = 0; should be sufficient...]
  }
	compute_ehh(data,*focl_snp,ANCSTRL,*nbr_snps,*nbr_chrs,tot_nbr_hapl,*min_nbr_hapl,*min_ehh,ehh); // compute the EHH for the ancestral allele
  compute_ehh(data,*focl_snp,DERIVED,*nbr_snps,*nbr_chrs,tot_nbr_hapl + *nbr_snps,*min_nbr_hapl,*min_ehh,ehh + *nbr_snps); // compute the EHH for the derived allele
  ihh[ANCSTRL] = integrate(map,ehh,*nbr_snps,*min_ehh,*max_gap,*discard_integration_at_border);         // compute the IHH for the ancestral allele
  ihh[DERIVED] = integrate(map,ehh + *nbr_snps,*nbr_snps,*min_ehh,*max_gap,*discard_integration_at_border);    // compute the IHH for the derived allele
  for (i = 0; i < *nbr_chrs; i++) {                                             // Free the memory for the `data' array
		free(data[i]);
	}
	free(data);
}

