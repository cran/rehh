#include "ehh_utils.h"

void compute_EHH(int **data,
				 int number_chromosomes,
				 int number_SNPs,
				 int focal,
				 int end,
				 int direction,
				 int allele,
				 double *ehh,
				 int *n)

{
	int i,j,k;
	int down,up;
	int haplotypes;
	int found;
	int missing;
	int *count;
	int **haplotype;
	
	if (direction == LEFT) {
		down = end;
		up = focal;
	} else {
		down = focal;
		up = end;		
	}
	haplotype = (int **) malloc(number_chromosomes * sizeof(int *));
	count = (int *) malloc(number_chromosomes * sizeof(int));
	haplotypes = 0;
	for (i = 0; i < number_chromosomes; i++) {
		if (data[i][focal] == allele) {
			missing = 0;
			j = down;
			while ((j <= up) && (!missing)) {
				if(data[i][j] == 9) {
					missing = 1;
				}
				++j;
			}
			if (!missing) {
				if (haplotypes == 0) {
					count[haplotypes] = 1;
					haplotype[haplotypes] = data[i];
					++haplotypes;
				} else {
					j = found = 0;
					while ((j < haplotypes) && (!found)) {
						k = down;
						while((data[i][k] == haplotype[j][k]) && k <= up) {
							++k;
						}
						if (k > up) {
							++count[j];
							found = 1;
						}
						++j;
					}
					if (!found) {
						count[haplotypes] = 1;
						haplotype[haplotypes] = data[i];
						++haplotypes;						
					}
				}
			}
		}
	}	
	*n = 0;
	*ehh = 0.0;
	if (haplotypes > 0) {
		for (i = 0; i < haplotypes; i++) {
			*n += count[i];
			*ehh += count[i] * (count[i] - 1.0);
		}
		if (*n > 1) {
			*ehh /= (*n * (*n - 1.0));
		} else {
			*ehh = 0.0;
		}
	}
	free(haplotype);
	free(count);
}

void compute_EHHS(int **data,
				  int number_chromosomes,
				  int number_SNPs,
				  int focal,
				  int end,
				  int direction,
				  double *ehhs,
				  int *n)

{
	int i,j,k;
	int down,up;
	int haplotypes;
	int found;
	int missing;
	int **haplotype;
	int *count;
	int allele_count[2];
	double haplotype_homozygosity;
	double homozygosity;
	
	if (direction == LEFT) {
		down = end;
		up = focal;
	} else {
		down = focal;
		up = end;		
	}	
	haplotype = (int **) malloc(number_chromosomes * sizeof(int *));
	count = (int *) malloc(number_chromosomes * sizeof(int));
	haplotypes = 0;	
	for (i = 0; i < number_chromosomes; i++) {
		missing = 0;
		j = down;
		while ((j <= up) && (!missing)) {
			if(data[i][j] == 9) {
				missing = 1;
			}
			++j;
		}
		if (!missing) {
			if (haplotypes == 0) {
				count[haplotypes] = 1;
				haplotype[haplotypes] = data[i];
				++haplotypes;
			} else {
				j = found = 0;
				while ((j < haplotypes) && (!found)) {
					k = down;
					while((data[i][k] == haplotype[j][k]) && k <= up) {
						++k;
					}
					if (k > up) {
						++count[j];
						found = 1;
					}
					++j;
				}
				if (!found) {
					count[haplotypes] = 1;
					haplotype[haplotypes] = data[i];
					++haplotypes;						
				}
			}
		}
	}	
	*n = 0;
	*ehhs = 0.0;
	allele_count[0] = allele_count[1] = 0;
	haplotype_homozygosity = 0.0;
	if (haplotypes > 0) {
		for (i = 0; i < haplotypes; i++) {
			*n += count[i];
			haplotype_homozygosity += pow(count[i],2);
			if (haplotype[i][focal] == 0) allele_count[0] += count[i];
			if (haplotype[i][focal] == 1) allele_count[1] += count[i];
		}
		homozygosity = pow(allele_count[0],2) + pow(allele_count[1],2);
		haplotype_homozygosity /= pow(*n,2);
		homozygosity /= pow(*n,2);
		*ehhs = ((*n * haplotype_homozygosity) - 1) / ((*n * homozygosity) - 1);	
	}
	free(haplotype);
	free(count);
}

double integrate(double *x_axis,
				 double *y_axis,
				 int n,
				 double threshold)

{
	int i;
	double height,width;
	double area = 0.0;
	
	for (i = 0; i < (n - 1); i++) {
		if ((y_axis[i] > threshold) || (y_axis[i + 1] > threshold)) {
			if ((y_axis[i] > threshold) && (y_axis[i + 1] > threshold)) {
				height = y_axis[i] + y_axis[i + 1] - 2.0 * threshold;
				width = x_axis[i + 1] - x_axis[i];
				area += width * height / 2;				
			} else {
				if (y_axis[i] > threshold) {
					height = y_axis[i] - threshold;
					width = ((x_axis[i + 1] * y_axis[i] - x_axis[i] * y_axis[i + 1]) + threshold * (x_axis[i] - x_axis[i + 1])) / (y_axis[i] - y_axis[i + 1]) - x_axis[i];
					area += width * height / 2;				
				}
				if (y_axis[i + 1] > threshold) {
					height = y_axis[i + 1] - threshold;
					width = x_axis[i + 1] - ((x_axis[i + 1] * y_axis[i] - x_axis[i] * y_axis[i + 1]) + threshold * (x_axis[i] - x_axis[i + 1]))/ (y_axis[i] - y_axis[i + 1]);
					area += width * height / 2;						
				}
				
			}
		}
	}			
	return(area);
}

