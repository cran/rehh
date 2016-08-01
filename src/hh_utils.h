#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define DEBUG

#define ANCSTRL  0
#define DERIVED  1
#define MISSING  9
#define UNDEFND -1

void compute_ehh(short int **data,int focl_snp, int allele,int nbr_snps,int nbr_chrs,int *tot_nbr_hapl,int min_nbr_hapl,double min_ehh,double *ehh);
void compute_ehhs(short int **data,int focl_snp,int nbr_snps,int nbr_chrs,int *tot_nbr_hapl,int min_nbr_hapl,double min_ehhs,double *ehhs_tang,double *ehhs_sabeti);
void init_hapl(short int **data,int snp,int frst_allele,int scnd_allele,int nbr_chrs,int *haplotype,int *nbr_dstnct_hapl,int *refrnce_hapl,int *hapl_count);
void updt_hapl(short int **data,int snp,int nbr_chrs,int *tot_nbr_hapl,int *haplotype,int *nbr_dstnct_hapl,int *refrnce_hapl,int *hapl_count);
double hapl_homzgsty(int tot_nbr_hapl,int nbr_dstnct_hapl,int *hapl_count);
double site_hapl_homzgsty(int tot_nbr_hapl,int nbr_dstnct_hapl,int *hapl_count,int *hapl_snp);
//double site_hapl_homzgsty_sabeti(int tot_nbr_hapl,int nbr_dstnct_hapl,int *hapl_count,int *hapl_snp);
double integrate(double *x_axis, double *y_axis, int n, double threshold);
