#include "hh_utils.h"

void compute_ehh(short int **data,
                 int focl_snp,
                 int allele,
                 int nbr_snps,
                 int nbr_chrs,
                 int *tot_nbr_hapl,
                 int min_nbr_hapl,
                 double min_ehh,
                 double *ehh)

{
  int i,j;
  int nbr_dstnct_hapl;                                                          // Number of distinct haplotypes
  int *haplotype;
  int *refrnce_hapl;
  int *hapl_count;

  haplotype = (int *) malloc(nbr_chrs * sizeof(int));                           // `haplotype' is the list of all haplotypes
  refrnce_hapl = (int *) malloc(nbr_chrs * sizeof(int));                        // `refrnce_hapl' is the list of referenced haplotypes
  hapl_count = (int *) malloc(nbr_chrs * sizeof(int));                          // `hapl_count' gives the counts of distinct haplotypes
  for (i = 0; i < nbr_chrs; i++) {
    if (data[i][focl_snp] == allele) tot_nbr_hapl[focl_snp] += 1;               // Count the number of haplotypes at the focal SNP that contain the `allele'
  }
  if (tot_nbr_hapl[focl_snp] < min_nbr_hapl) {                                  // Except if the total number of haplotypes is less than the minimum allowed, ...
    ehh[focl_snp] = 0.0;                                                        // ... in which case the EHH at the focal SNP is nought, ...
  }
  else {
    ehh[focl_snp] = 1.0;                                                        //... the EHH at the focal SNP is 1.0, by definition
    init_hapl(data,focl_snp,allele,UNDEFND,nbr_chrs,haplotype,&nbr_dstnct_hapl,refrnce_hapl,hapl_count); // initialize the array of haplotypes (for the focal SNP) that contain the `allele'
#ifdef DEBUG
    {
      for (i = 0; i < nbr_chrs; i++) {
        printf("[pos. %d] haplotype[%d] = %d\n",focl_snp,i,haplotype[i]);
      }
      printf("\n\toverll nbr of haplotypes = %d\n",tot_nbr_hapl[focl_snp]);
      printf("\tnbr of dstnct haplotypes = %d\n",nbr_dstnct_hapl);
      for (i = 0; i < nbr_dstnct_hapl; i++) {
        printf("\trefrnce_hapl[%d] = %d (%d)\n",i,refrnce_hapl[i],hapl_count[i]);
      }
      printf("\n");
    }
#endif
    for (j = (focl_snp - 1); j >= 0; j--) {                                     // walk along the chromosome from the focal SNP to the left-hand side, ...
      updt_hapl(data,j,nbr_chrs,&tot_nbr_hapl[j],haplotype,&nbr_dstnct_hapl,refrnce_hapl,hapl_count); // ... and update the list of haplotypes
#ifdef DEBUG
      {
        for (i = 0; i < nbr_chrs; i++) {
          printf("[pos. %d] haplotype[%d] = %d\n",j,i,haplotype[i]);
        }
        printf("\n\toverll nbr of haplotypes = %d\n",tot_nbr_hapl[j]);
        printf("\tnbr of dstnct haplotypes = %d\n",nbr_dstnct_hapl);
        for (i = 0; i < nbr_dstnct_hapl; i++) {
          printf("\trefrnce_hapl[%d] = %d (%d)\n",i,refrnce_hapl[i],hapl_count[i]);
        }
        printf("\n");
      }
#endif
      ehh[j] = hapl_homzgsty(tot_nbr_hapl[j],nbr_dstnct_hapl,hapl_count);       // Compute the haplotype homozygosity at the jth SNP
      if (tot_nbr_hapl[j] >= min_nbr_hapl) {                                    // If the total number of haplotypes at the jth SNP is larger than the min. number allowed, ...
        if (ehh[j] <= min_ehh) {                                                 // ... but if the EHH is lower than the min. number allowed, ...
          ehh[j] = 0.0;                                                         // ... then the EHH is nought, and ...
          break;                                                                // ... quit the loop
        }
      }
      else {                                                                    // If the total number of haplotypes at the jth SNP is lower than the min. number allowed, ...
        ehh[j] = 0.0;                                                           // ... then the EHH is nought, and ...
        break;                                                                  // ... quit the loop
      }
    }
    init_hapl(data,focl_snp,allele,UNDEFND,nbr_chrs,haplotype,&nbr_dstnct_hapl,refrnce_hapl,hapl_count); // initialize the array of haplotypes
#ifdef DEBUG
    {
      for (i = 0; i < nbr_chrs; i++) {
        printf("[pos. %d] haplotype[%d] = %d\n",focl_snp,i,haplotype[i]);
      }
      printf("\n\toverll nbr of haplotypes = %d\n",tot_nbr_hapl[focl_snp]);
      printf("\tnbr of dstnct haplotypes = %d\n",nbr_dstnct_hapl);
      for (i = 0; i < nbr_dstnct_hapl; i++) {
        printf("\trefrnce_hapl[%d] = %d (%d)\n",i,refrnce_hapl[i],hapl_count[i]);
      }
      printf("\n");
    }
#endif
    for (j = (focl_snp + 1); j < nbr_snps; j++) {                               // walk along the chromosome from the focal SNP to the right-hand side
      updt_hapl(data,j,nbr_chrs,&tot_nbr_hapl[j],haplotype,&nbr_dstnct_hapl,refrnce_hapl,hapl_count); // ... and update the list of haplotypes
#ifdef DEBUG
      {
        for (i = 0; i < nbr_chrs; i++) {
          printf("[pos. %d] haplotype[%d] = %d\n",j,i,haplotype[i]);
        }
        printf("\n\toverll nbr of haplotypes = %d\n",tot_nbr_hapl[j]);
        printf("\tnbr of dstnct haplotypes = %d\n",nbr_dstnct_hapl);
        for (i = 0; i < nbr_dstnct_hapl; i++) {
          printf("\trefrnce_hapl[%d] = %d (%d)\n",i,refrnce_hapl[i],hapl_count[i]);
        }
        printf("\n");
      }
#endif
      ehh[j] = hapl_homzgsty(tot_nbr_hapl[j],nbr_dstnct_hapl,hapl_count);       // Compute the haplotype homozygosity at the jth SNP
      if (tot_nbr_hapl[j] >= min_nbr_hapl) {                                    // If the total number of haplotypes at the jth SNP is larger than the min. number allowed, ...
        if (ehh[j] <= min_ehh) {                                                 // ... but if the EHH is lower than the min. number allowed, ...
          ehh[j] = 0.0;                                                         // ... then the EHH is nought, and ...
          break;                                                                // ... quit the loop
        }
      }
      else {                                                                    // If the total number of haplotypes at the jth SNP is lower than the min. number allowed, ...
        ehh[j] = 0.0;                                                           // ... then the EHH is nought, and ...
        break;                                                                  // ... quit the loop
      }
    }
  }
  free(haplotype);
  free(refrnce_hapl);
  free(hapl_count);
}

void compute_ehhs(short int **data,
                  int focl_snp,
                  int nbr_snps,
                  int nbr_chrs,
                  int *tot_nbr_hapl,
                  int min_nbr_hapl,
                  double min_ehhs,
                  double *ehhs_tang,
                  double *ehhs_sabeti)

{
  int i,j;
  int nbr_dstnct_hapl;                                                          // Number of distinct haplotypes
  int *haplotype;
  int *refrnce_hapl;
  int *hapl_count;
  int *hapl_snp;

  haplotype = (int *) malloc(nbr_chrs * sizeof(int));                           // `haplotype' is the list of all haplotypes
  refrnce_hapl = (int *) malloc(nbr_chrs * sizeof(int));                        // `refrnce_hapl' is the list of referenced haplotypes
  hapl_count = (int *) malloc(nbr_chrs * sizeof(int));                          // `hapl_count' gives the counts of distinct haplotypes
  hapl_snp = (int *) malloc(nbr_chrs * sizeof(int));
  for (i = 0; i < nbr_chrs; i++) {
    if ((data[i][focl_snp] == ANCSTRL) || (data[i][focl_snp] == DERIVED)) tot_nbr_hapl[focl_snp] += 1; // Count the number of haplotypes at the focal SNP
  }
  if (tot_nbr_hapl[focl_snp] < min_nbr_hapl) {                                  // Except if the total number of haplotypes is less than the minimum allowed, ...
    ehhs_tang[focl_snp] = 0.0;                                                  // ... in which case the EHHS (Tang et al's 2007) at the focal SNP is nought, ...
    ehhs_sabeti[focl_snp] = 0.0;                                                // ... in which case the EHHS (Sabeti et al's 2007) at the focal SNP is nought, ...
  }
  else {
    init_hapl(data,focl_snp,ANCSTRL,DERIVED,nbr_chrs,haplotype,&nbr_dstnct_hapl,refrnce_hapl,hapl_count); // initialize the array of haplotypes (for the focal SNP)
    ehhs_tang[focl_snp] = 1.0;                                                  // 1.0 by definition
    ehhs_sabeti[focl_snp] = hapl_homzgsty(tot_nbr_hapl[focl_snp],nbr_dstnct_hapl,hapl_count); // Compute the standardized haplotype homozygosity at the focal SNP, following Sabeti et al's (2007)
    for (j = (focl_snp - 1); j >= 0; j--) {                                     // walk along the chromosome from the focal SNP to the left-hand side, ...
      updt_hapl(data,j,nbr_chrs,&tot_nbr_hapl[j],haplotype,&nbr_dstnct_hapl,refrnce_hapl,hapl_count); // ... and update the list of haplotypes
      for (i = 0; i < nbr_dstnct_hapl; i++) {
        hapl_snp[i] = data[ refrnce_hapl[i] ][focl_snp];
      }
      ehhs_tang[j] = site_hapl_homzgsty(tot_nbr_hapl[j],nbr_dstnct_hapl,hapl_count,hapl_snp); // Compute the haplotype homozygosity at the jth SNP, following Tang et al's (2007)
      ehhs_sabeti[j] = hapl_homzgsty(tot_nbr_hapl[j],nbr_dstnct_hapl,hapl_count); // Compute the standardized haplotype homozygosity at the jth SNP, following Sabeti et al's (2007)
      if (tot_nbr_hapl[j] >= min_nbr_hapl) {                                    // If the total number of haplotypes at the jth SNP is larger than the min. number allowed, ...
        if ((ehhs_tang[j] <= min_ehhs) && (ehhs_sabeti[j] <= min_ehhs)) {         // ... but if the EHHS is lower than the min. number allowed, ...
          ehhs_tang[j] = 0.0;                                                   // ... then the EHHS (Tang et al's 2007) is nought, and ...
          ehhs_sabeti[j] = 0.0;                                                 // ... then the EHHS (Sabeti et al's 2007) is nought, and ...
          break;                                                                // ... then quit the loop
        }
        else if (ehhs_tang[j] <= min_ehhs) {
          ehhs_tang[j] = 0.0;
        }
        else if (ehhs_sabeti[j] <= min_ehhs) {
          ehhs_sabeti[j] = 0.0;
        }
        else {
          continue;
        }
      }
      else {                                                                    // If the total number of haplotypes at the jth SNP is lower than the min. number allowed, ...
        ehhs_tang[j] = 0.0;                                                     // ... then the EHHS (Tang et al's 2007) is nought, and ...
        ehhs_sabeti[j] = 0.0;                                                   // ... then the EHHS (Sabeti et al's 2007) is nought, and ...
        break;                                                                  // ... quit the loop
      }
    }
    init_hapl(data,focl_snp,ANCSTRL,DERIVED,nbr_chrs,haplotype,&nbr_dstnct_hapl,refrnce_hapl,hapl_count); // initialize the array of haplotypes (for the focal SNP)
    for (j = (focl_snp + 1); j < nbr_snps; j++) {                               // walk along the chromosome from the focal SNP to the right-hand side
      updt_hapl(data,j,nbr_chrs,&tot_nbr_hapl[j],haplotype,&nbr_dstnct_hapl,refrnce_hapl,hapl_count); // ... and update the list of haplotypes
      for (i = 0; i < nbr_dstnct_hapl; i++) {
        hapl_snp[i] = data[ refrnce_hapl[i] ][focl_snp];
      }
      ehhs_tang[j] = site_hapl_homzgsty(tot_nbr_hapl[j],nbr_dstnct_hapl,hapl_count,hapl_snp); // Compute the standardized haplotype homozygosity at the jth SNP
      ehhs_sabeti[j] = hapl_homzgsty(tot_nbr_hapl[j],nbr_dstnct_hapl,hapl_count); // Compute the standardized haplotype homozygosity at the jth SNP, following Sabeti et al's (2007)
      if (tot_nbr_hapl[j] >= min_nbr_hapl) {                                    // If the total number of haplotypes at the jth SNP is larger than the min. number allowed, ...
        if ((ehhs_tang[j] <= min_ehhs) && (ehhs_sabeti[j] <= min_ehhs)) {         // ... but if the EHHS is lower than the min. number allowed, ...
          ehhs_tang[j] = 0.0;                                                   // ... then the EHHS (Tang et al's 2007) is nought, and ...
          ehhs_sabeti[j] = 0.0;                                                 // ... then the EHHS (Sabeti et al's 2007) is nought, and ...
          break;                                                                // ... then quit the loop
        }
        else if (ehhs_tang[j] <= min_ehhs) {
          ehhs_tang[j] = 0.0;
        }
        else if (ehhs_sabeti[j] <= min_ehhs) {
          ehhs_sabeti[j] = 0.0;
        }
        else {
          continue;
        }
      }
      else {                                                                    // If the total number of haplotypes at the jth SNP is lower than the min. number allowed, ...
        ehhs_tang[j] = 0.0;                                                     // ... then the EHHS (Tang et al's 2007) is nought, and ...
        ehhs_sabeti[j] = 0.0;                                                   // ... then the EHHS (Sabeti et al's 2007) is nought, and ...
        break;                                                                  // ... quit the loop
      }
    }
  }
  free(haplotype);
  free(refrnce_hapl);
  free(hapl_count);
  free(hapl_snp);
}

void init_hapl(short int **data,
               int snp,
               int frst_allele,
               int scnd_allele,
               int nbr_chrs,
               int *haplotype,
               int *nbr_dstnct_hapl,
               int *refrnce_hapl,
               int *hapl_count)


{
  int i;
  int frst_allele_hapl = UNDEFND;                                               // Initialize the haplotype that contains the the `frst_allele'
  int scnd_allele_hapl = UNDEFND;                                               // Initialize the haplotype that contains the the `scnd_allele'

  *nbr_dstnct_hapl = 0;                                                         // Initialize the number of distinct haplotypes
  hapl_count[ *nbr_dstnct_hapl ] = 0;                                           // Initialize the haplotype counts
  for (i = 0; i < nbr_chrs; i++) {                                              // Loop over all chromosomes in the sample
    if (data[i][snp] == frst_allele) {                                          // If the focal SNP from the ith chromosome contains the `frst_allele', ...
      if (frst_allele_hapl == UNDEFND) {                                        // ... and if the `frst_allele_hapl' has not been defined yet, ...
        refrnce_hapl[ *nbr_dstnct_hapl ] = i;                                   // ... then, `refrnce_hapl' refers to the ith chromosome, ...
        frst_allele_hapl = *nbr_dstnct_hapl;                                    // ... and the `frst_allele_hapl' is defined, ...
        hapl_count[ *nbr_dstnct_hapl ] += 1;                                    // ... and the haplotype count is increased, ...
        *nbr_dstnct_hapl += 1;                                                  // ... and the number of distinct haplotypes is incremented
        hapl_count[ *nbr_dstnct_hapl ] = 0;                                     // ... and the haplotype count for the next haplotype is initialized
      }
      else {
        hapl_count[ frst_allele_hapl ] += 1;                                    // ... otherwise (if the `frst_allele_hapl' has already been defined), then increment the haplotype count
      }
      haplotype[i] = frst_allele_hapl;                                          // The ith haplotype is defined
    }
    else if (data[i][snp] == scnd_allele) {                                     // If the focal SNP from the ith chromosome contains the `scnd_allele', ...
      if (scnd_allele_hapl == UNDEFND) {                                        // ... and if the `scnd_allele_hapl' has not been defined yet, ...
        refrnce_hapl[ *nbr_dstnct_hapl ] = i;                                   // ... then, `refrnce_hapl' refers to the ith chromosome, ...
        scnd_allele_hapl = *nbr_dstnct_hapl;                                    // ... and the `scnd_allele_hapl' is defined, ...
        hapl_count[ *nbr_dstnct_hapl ] += 1;                                    // ... and the haplotype count is increased, ...
        *nbr_dstnct_hapl += 1;                                                  // ... and the number of distinct haplotypes is incremented
        hapl_count[ *nbr_dstnct_hapl ] = 0;                                     // ... and the haplotype count for the next haplotype is initialized
      }
      else {
        hapl_count[ scnd_allele_hapl ] += 1;                                    // ... otherwise (if the `scnd_allele_hapl' has already been defined), then increment the haplotype count
      }
      haplotype[i] = scnd_allele_hapl;                                          // The ith haplotype is defined
    }
    else {
      haplotype[i] = UNDEFND;                                                   // If the focal SNP from the ith chromosome contains neither the `frst_allele', nor the `scnd_allele', then it is undefined
    }
  }
}

void updt_hapl(short int** data,
               int snp,
               int nbr_chrs,
               int *tot_nbr_hapl,
               int *haplotype,
               int *nbr_dstnct_hapl,
               int *refrnce_hapl,
               int *hapl_count)

{
  int i,k;
  int *derived_hapl = (int *) malloc(nbr_chrs * sizeof(int));
  int *ancstrl_hapl = (int *) malloc(nbr_chrs * sizeof(int));

  *tot_nbr_hapl = 0;                                                            // Initialize the total number of haplotypes
  for (i = 0; i < nbr_chrs; i++) {
    derived_hapl[i] = UNDEFND;                                                  // All derived haplotypes are set to UNDEFND
    ancstrl_hapl[i] = UNDEFND;
  }
  for (i = 0; i < nbr_chrs; i++) {
    if (haplotype[i] != UNDEFND) {                                              // The following is only for defined haplotypes (e.g., for haplotypes that carry the reference allele)
      if (data[i][snp] == MISSING) {                                            // If the data is missing at the SNP considered for the ith haplotype, ...
        if (refrnce_hapl[ haplotype[i] ] == i) {                                // ... and if the ith haplotype is a referenced copy of itself, ...
          for (k = 0; k < nbr_chrs; k++) {                                      // ... then loop over all chromosomes, ...
            if ((i != k) && (haplotype[k] == haplotype[i])) {                   // ... to look for another copy (k) of the referenced haplotype copy
              refrnce_hapl[ haplotype[i] ] = k;
              hapl_count[ haplotype[i] ] -= 1;                                  // ... and decrease the haplotype count for haplotype[i]
              break;
            }
          }
          if (k == nbr_chrs) {                                                  // If no copy of the referenced haplotype is found (which means that the referenced haploype no longer exists), ...
            for (k = 0; k < nbr_chrs; k++) {                                    // ... then loop over all chromosomes, ...
              if (haplotype[k] == (*nbr_dstnct_hapl - 1)) {                     // ... to replace the haplotype reference which no longer exists, with the last element of the list of reference haplotypes
                haplotype[k] = haplotype[i];
              }
            }

            if (ancstrl_hapl[ (*nbr_dstnct_hapl - 1) ] != UNDEFND) {            // If the haplotype reference that is being erased had a derived haplotype, ...
              derived_hapl[ ancstrl_hapl[ (*nbr_dstnct_hapl - 1) ] ] = haplotype[i]; // ... then re-assign the derived haplotype, ...
              ancstrl_hapl[ haplotype[i] ] = ancstrl_hapl[ (*nbr_dstnct_hapl - 1) ]; // ... and the ancestral haplotype
            }

            refrnce_hapl[ haplotype[i] ] = refrnce_hapl[ *nbr_dstnct_hapl - 1 ]; // Replace the (current) referenced haplotype with the last one in the list of reference haplotypes, ...
            hapl_count[ haplotype[i] ] = hapl_count[ *nbr_dstnct_hapl - 1 ];    // ... and take the corresponding haplotype count
            *nbr_dstnct_hapl -= 1;                                              // ... and decrease the number of distinct haplotypes by one
          }
          haplotype[i] = UNDEFND;                                               // ... then (if the data is missing at the SNP considered) the ith haplotype is undefined
        }
        else {                                                                  // ... otherwise (if the ith haplotype is NOT a referenced copy of itself), ...
          hapl_count[ haplotype[i] ] -= 1;                                      // ... then decrease the haplotype count for haplotype[i]
          haplotype[i] = UNDEFND;                                               // ... and set the ith haplotype to the undefined state
        }
      }
      else if (data[i][snp] != data[ refrnce_hapl[ haplotype[i] ] ][snp]) {     // If (at the position considered) the ith chromosome carries a new allele that differs from the referenced haplotype
        if (derived_hapl[ haplotype[i] ] != UNDEFND) {                          // ... and if a derived haplotype already exists, ...
          hapl_count[ haplotype[i] ] -= 1;                                      // ... then decrease the haplotype count for the `old' haplotype[i], ...
          haplotype[i] = derived_hapl[ haplotype[i] ];                          // ... and the ith haplotype is set to its derived haplotype, ...
          hapl_count[ haplotype[i] ] += 1;                                      // ... and increase the haplotype count for the `new' haplotype[i]
        }
        else {                                                                  // ... otherwise (if a derived haplotype does not exist yet), ...
          refrnce_hapl[ *nbr_dstnct_hapl ] = i;                                 // ... then `refrnce_hapl' refers to the ith chromosome, ...
          derived_hapl[ haplotype[i] ] = *nbr_dstnct_hapl;                      // ... define the derived haplotype, ...
          ancstrl_hapl[ *nbr_dstnct_hapl ] = haplotype[i];                      // ... and define the ancestral haplotype
          hapl_count[ haplotype[i] ] -= 1;                                      // ... then decrease the haplotype count for the `old' haplotype[i]
          haplotype[i] = derived_hapl[ haplotype[i] ];                          // ... the ith haplotype is set to its derived haplotype, ...
          hapl_count[ haplotype[i] ] = 1;                                       // ... and set the haplotype count for the `new' haplotype[i] to one, ...
          *nbr_dstnct_hapl += 1;                                                // ... and the number of distinct haplotypes is incremented
        }
        *tot_nbr_hapl += 1;                                                     // Increment the total number of haplotypes (for newly defined haplotypes)
      }
      else {
        *tot_nbr_hapl += 1;                                                     // Increment the total number of haplotypes (for haplotypes that already exist)
      }
    }
  }
	free(derived_hapl);
  free(ancstrl_hapl);
}

double hapl_homzgsty(int tot_nbr_hapl,
                     int nbr_dstnct_hapl,
                     int *hapl_count)

{
  int i;
  double homzgsty;

  homzgsty = 0.0;                                                               // Initialize the homozygosity
  for (i = 0; i < nbr_dstnct_hapl; i++) {                                       // Loop over all distinct haplotypes, ...
    homzgsty += (double) hapl_count[i] * (hapl_count[i] - 1.0);                 // ... to compute the homozygosity
  }
  if (tot_nbr_hapl > 1) {                                                       // If the total number of haplotypes is larger than 1, ...
    homzgsty /= (double) (tot_nbr_hapl * (tot_nbr_hapl - 1.0));                 // ... then finalize the computation of homozygosity, ...
  }
  else {
    homzgsty = 0.0;                                                             // ... otherwise, homozygosity is nought
  }
  return(homzgsty);
}

double site_hapl_homzgsty(int tot_nbr_hapl,
                          int nbr_dstnct_hapl,
                          int *hapl_count,
                          int *hapl_snp)

{
  int i;
  int allele_count[2];
  double hapl_homzgsty,homzgsty;

  hapl_homzgsty = 0.0;                                                          // Initialize the haplotype homozygosity
  allele_count[ANCSTRL] = allele_count[DERIVED] = 0;                            // Initialize the allele counts
  if (nbr_dstnct_hapl > 0) {
    for (i = 0; i < nbr_dstnct_hapl; i++) {                                     // Loop over all distinct haplotypes, ...
      hapl_homzgsty += pow(hapl_count[i],2.0);                                  // ... to compute the haplotype homozygosity, ...
      if (hapl_snp[i] == ANCSTRL) allele_count[ANCSTRL] += hapl_count[i];       // ... and the allele counts for each of the ANCESTRL and DERIVED SNPs
      if (hapl_snp[i] == DERIVED) allele_count[DERIVED] += hapl_count[i];
    }
    hapl_homzgsty /= pow(tot_nbr_hapl,2.0);                                     // Finalize the computation of haplotype homozygosity
    homzgsty = pow(allele_count[ANCSTRL],2.0) + pow(allele_count[DERIVED],2.0);
    homzgsty /= pow(tot_nbr_hapl,2.0);
    return ((((double) tot_nbr_hapl * hapl_homzgsty) - 1.0) / (((double) tot_nbr_hapl * homzgsty) - 1.0));
  }
  else {
    return (0.0);
  }
}

double integrate(double *x_axis,
                 double *y_axis,
                 int n,
                 double threshold,
                 double max_gap,
                 int discard_integration_at_border)

{
  int i;
  double height,width;
  double area = 0.0;

  if (discard_integration_at_border && ((y_axis[0] > threshold) || (y_axis[n - 1] > threshold))) {  // If the EHH or EHHS is larger than the minimum value at either end of the chromosome, ...
    return (UNDEFND);                                                           // ... then do not compute the integral, and quit
  }
  for (i = 0; i < (n - 1); i++) {
    if ((y_axis[i] > threshold) || (y_axis[i + 1] > threshold)) {
      if (fabs(x_axis[i + 1] - x_axis[i]) > max_gap) {                          // If a gap larger than max_gap exists within the 'support' of the EHH or EHHS, ...
        return (UNDEFND);                                                       // ... then do not compute the integral, and quit
      }
      if ((y_axis[i] > threshold) && (y_axis[i + 1] > threshold)) {
        height = y_axis[i] + y_axis[i + 1] - 2.0 * threshold;
        width = x_axis[i + 1] - x_axis[i];
        area += width * height / 2;
      }
      else {
        if (y_axis[i] > threshold) {
          height = y_axis[i] - threshold;
          width = ((x_axis[i + 1] * y_axis[i] - x_axis[i] * y_axis[i + 1]) + threshold * (x_axis[i] - x_axis[i + 1]))	/ (y_axis[i] - y_axis[i + 1]) - x_axis[i];
          area += width * height / 2;
        }
        if (y_axis[i + 1] > threshold) {
          height = y_axis[i + 1] - threshold;
          width = x_axis[i + 1] - ((x_axis[i + 1] * y_axis[i] - x_axis[i] * y_axis[i + 1]) + threshold * (x_axis[i] - x_axis[i + 1])) / (y_axis[i] - y_axis[i + 1]);
          area += width * height / 2;
        }
      }
    }
  }
  return (area);
}
