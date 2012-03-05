#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "ehh_utils.h"

void r_ehh(int *Rdata,
		   int *number_SNPs,
		   int *number_chromosomes,
		   int *focal_SNP,
		   double *map,
		   int *number_haplotypes,
		   double *EHH,
		   double *IHH,
		   int *min_number_haplotypes,
		   double *min_EHH)

{
	int i,j;
	int allele;
	int n;
	int **data;
	double ehh;
	double *tmp;

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
		for (allele = 0; allele < 2; allele++) {
			EHH[(allele * *number_SNPs) + j] = 0.0;
			number_haplotypes[(allele * *number_SNPs) + j] = 0;
		}
	}
	EHH[*focal_SNP] = EHH[*number_SNPs + *focal_SNP] = 1.0;
	number_haplotypes[*focal_SNP] = number_haplotypes[*number_SNPs + *focal_SNP] = 0;
	for (i = 0; i < *number_chromosomes; i++) {
		if (data[i][*focal_SNP] == 0) number_haplotypes[*focal_SNP] += 1;
		if (data[i][*focal_SNP] == 1) number_haplotypes[*number_SNPs + *focal_SNP] += 1;
	}
	for (allele = 0; allele < 2; allele++) {
		if (number_haplotypes[(allele * *number_SNPs) + *focal_SNP] < *min_number_haplotypes) {
			EHH[(allele * *number_SNPs) + *focal_SNP] = 0.0;
		} else {
			for (j = (*focal_SNP - 1); j >= 0; j--) {
				compute_EHH(data,*number_chromosomes,*number_SNPs,*focal_SNP,j,LEFT,allele,&ehh,&n);
				if (n > *min_number_haplotypes) { 
					EHH[(allele * *number_SNPs) + j] = ehh;
					number_haplotypes[(allele * *number_SNPs) + j] = n;
						if(ehh < *min_EHH){ 
								break;
								}
				} else {
					break;
				}


			}
			for (j = (*focal_SNP + 1); j < *number_SNPs; j++) {
				compute_EHH(data,*number_chromosomes,*number_SNPs,*focal_SNP,j,RIGHT,allele,&ehh,&n);
				if (n > *min_number_haplotypes) { 
					EHH[(allele * *number_SNPs) + j] = ehh;
					number_haplotypes[(allele * *number_SNPs) + j] = n;
						if(ehh < *min_EHH){ 
								break;
								}
				} else {
					break;
				}


			}
		}
	}
	tmp = (double *) malloc(*number_SNPs * sizeof(double));
	for (allele = 0; allele < 2; allele++) {
		for (j = 0; j < *number_SNPs; j++) {
			tmp[j] = EHH[(allele * *number_SNPs) + j];
		}
		IHH[allele] = integrate(map,tmp,*number_SNPs,*min_EHH);
	}	
	for (i = 0; i < *number_chromosomes; i++) {
		free(data[i]);
	}
	free(data);
	free(tmp);
}
