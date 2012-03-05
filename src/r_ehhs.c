#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "ehh_utils.h"

void r_ehhs(int *Rdata,
			int *number_SNPs,
			int *number_chromosomes,
			int *focal_SNP,
			double *map,
			int *number_haplotypes,
			double *EHHS,
			double *IES,
			int *min_number_haplotypes,
			double *min_EHH)

{
	int i,j;
	int n;
	int **data;
	double ehhs;
	
	*focal_SNP -= 1;	
	data = (int **) malloc(*number_chromosomes * sizeof(int *));
	for (i = 0; i < *number_chromosomes; i++) {
		data[i] = (int *) malloc(*number_SNPs * sizeof(int));
	}
	for (i = 0; i < *number_chromosomes; i++) {
		for (j = 0; j < *number_SNPs; j++) {	
			if (Rdata[(j * *number_chromosomes) + i] == 0) data[i][j] = 9;
			if (Rdata[(j * *number_chromosomes) + i] == 1) data[i][j] = 0;
			if (Rdata[(j * *number_chromosomes) + i] == 2) data[i][j] = 1;
		}
	}
	for (j = 0; j < *number_SNPs; j++) {	
		EHHS[j] = 0.0;
		number_haplotypes[j] = 0;
	}
	EHHS[*focal_SNP] = 1.0;
	number_haplotypes[*focal_SNP] = 0;
	for (i = 0; i < *number_chromosomes; i++) {
		if (data[i][*focal_SNP] == 0) number_haplotypes[*focal_SNP] += 1;
		if (data[i][*focal_SNP] == 1) number_haplotypes[*focal_SNP] += 1;
	}
	if (number_haplotypes[*focal_SNP] < *min_number_haplotypes) {
		EHHS[*focal_SNP] = 0.0;
	} else {
		for (j = (*focal_SNP - 1); j >= 0; j--) {
			compute_EHHS(data,*number_chromosomes,*number_SNPs,*focal_SNP,j,LEFT,&ehhs,&n);
					if (n > *min_number_haplotypes) {
						EHHS[j] = ehhs;
						number_haplotypes[j] = n;
							if(ehhs < *min_EHH){ 
									break;
									}
					} else {
						break;
					}


		}
		for (j = (*focal_SNP + 1); j < *number_SNPs; j++) {
			compute_EHHS(data,*number_chromosomes,*number_SNPs,*focal_SNP,j,RIGHT,&ehhs,&n);
					if (n > *min_number_haplotypes) {
						EHHS[j] = ehhs;
						number_haplotypes[j] = n;
							if(ehhs < *min_EHH){ 
									break;
									}
					} else {
						break;
					}

		}
	}
	*IES = integrate(map,EHHS,*number_SNPs,*min_EHH);
	for (i = 0; i < *number_chromosomes; i++) {
		free(data[i]);
	}
	free(data);	
}
