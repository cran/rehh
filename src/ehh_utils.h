#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define LEFT 0
#define RIGHT 1

double integrate(double *x_axis,double *y_axis,int n,double threshold);
void compute_EHH(int **data,int number_chromosomes,int number_SNPs,int focal,int end,int direction,int allele,double *ehh,int *n);
void compute_EHHS(int **data,int number_chromosomes,int number_SNPs,int focal,int end,int direction,double *ehhs,int *n);

