#include "calc_sfs_tests.h"
#include "sfs_moments.h"


void calc_sfs_tests(int* data, int nbr_chr, int nbr_mrk, double* map, 
                    bool polarized, double* windows, int n_windows, bool right, 
                    int min_n_mrk, int* n_mrk, double* results){
  int n_tests = 3;
  int fs_size;
  
  double* T = (double*) malloc(n_tests * sizeof(double));
  double* alpha = (double*) malloc(n_tests * sizeof(double));
  double* beta = (double*) malloc(n_tests * sizeof(double));
  double** secMom0;
  
  double*** Omega = (double***) malloc(n_tests * sizeof(double**));
  
  for(int test = 0; test < n_tests; test++){
    Omega[test] = (double**) malloc(nbr_chr * sizeof(double*));
    for(int n = 1; n < nbr_chr; n++){
      /* weights for sample size n */
      Omega[test][n] = (double*) malloc(n * sizeof(double));
    }
  }
  
  double** omega_S = (double**) malloc(nbr_chr * sizeof(double*));
  double** omega_PI = (double**) malloc(nbr_chr * sizeof(double*));
  double** omega_L = (double**) malloc(nbr_chr * sizeof(double*));
  double** fs0 = (double**) malloc(nbr_chr * sizeof(double*));
  
  /* get weights and null spectrum for each possible (sub-)sample size */
  if (polarized) {
    for(int n = 2; n <= nbr_chr; n++){
      omega_S[n - 1] = getOmega(n, 'S');
      omega_PI[n - 1] = getOmega(n, 'P');
      omega_L[n - 1] = getOmega(n, 'L');
      fs0[n - 1] = getXi0(n);
    }
    secMom0 = getSigma(nbr_chr);
  }else{
    for(int n = 2; n <= nbr_chr; n++){
      omega_S[n - 1] = getOmegaStar(n,'S');
      omega_PI[n - 1] = getOmegaStar(n, 'P');
      omega_L[n - 1] = getOmegaStar(n, 'L');
      fs0[n - 1] = getEta0(n);
    }
    secMom0 = getRho(nbr_chr);
  }
  
  for(int n = 1; n < nbr_chr; n++){
    for(int i = 0; i < n; i++){
      Omega[0][n][i] = omega_PI[n][i] - omega_S[n][i];
      Omega[1][n][i] = omega_PI[n][i] - omega_L[n][i];
      Omega[2][n][i] = omega_L[n][i] - omega_S[n][i];
    }
  }
  
  if (polarized) {
    fs_size = nbr_chr - 1;
  }else{
    fs_size = nbr_chr / 2;
  }
  
  /*
   * compute coefficients for theta-square estimation
   */
  double an = getWeightedFirstMoment(omega_S[nbr_chr - 1], fs0[nbr_chr - 1], fs_size);
  double bn = getWeightedSecondMoment(omega_S[nbr_chr - 1], fs0[nbr_chr - 1], secMom0, fs_size);
  
  /*
   * compute coefficients for the variance in the nominator
   */
  for (int test = 0; test < n_tests; test++) {
    alpha[test] = getWeightedFirstMoment(Omega[test][nbr_chr - 1], fs0[nbr_chr - 1], fs_size);
    beta[test] = getWeightedSecondMoment(Omega[test][nbr_chr - 1], fs0[nbr_chr - 1], secMom0, fs_size);
  }
  
  int start_mrk = 0;
  for(int window = 0; window < n_windows; window++){
    /* increase marker until it exceeds left window boundary */
    while(((right && map[start_mrk] <= windows[window])
             ||(!right && map[start_mrk] < windows[window]))
            && start_mrk < nbr_mrk){
            start_mrk++;
    }
    
    double theta_S = 0.0, theta_PI = 0.0, theta_L = 0.0; 
    
    int mrk = start_mrk;
    int n_polymorphic_mrk = 0;
    
    /* increase marker until it exceeds right window boundary */
    while(((right && map[mrk] <= windows[n_windows + window])
             || (!right && map[mrk] < windows[n_windows + window]))
            && mrk < nbr_mrk){
            
      int x = 0; 
      int n = 0;
      int first_non_missing_allele = -1;
      for(int chr = 0; chr < nbr_chr; chr++){
        if(data[mrk * nbr_chr + chr] != MISSING_VALUE){
          n++;
          if(polarized){
            /* any non-zero allele is taken as derived */ 
            if(data[mrk * nbr_chr + chr] > 0){
              x++;
            }
            /* count first non-missing allele */
          }else{
            if(first_non_missing_allele == -1){
              first_non_missing_allele = data[mrk * nbr_chr + chr];
              x++;
            }else if(data[mrk * nbr_chr + chr] == first_non_missing_allele){
              x++;
            }
          }
        }
      }
      /* fold spectrum */
      if(!polarized && x > (n + 1) / 2){
        x = n - x;
      }
      
      /* if polymorphic ... */
      if(x > 0 && x < n){
        n_polymorphic_mrk++;
        theta_S += (1. / fs0[n - 1][x - 1] ) * omega_S[n - 1][x - 1];
        theta_PI += (1. / fs0[n - 1][x - 1]) * omega_PI[n - 1][x - 1];
        theta_L += (1. / fs0[n - 1][x - 1]) * omega_L[n - 1][x - 1];
      }
      mrk++;
    }
    
    if(n_polymorphic_mrk > min_n_mrk){
      results[window] = theta_S;
      results[n_windows + window] = theta_PI;
      results[n_windows * 2 + window] = theta_L;
      
      double theta_S_squared = (theta_S * theta_S - an * theta_S) / (1. + bn);
      
      T[0] = theta_PI - theta_S;
      T[1] = theta_PI - theta_L;
      T[2] = theta_L - theta_S;
      
      for(int test = 0; test < n_tests; test++){
        if (T[test] != 0) {
          T[test] /= sqrt(alpha[test] * theta_S + beta[test] * theta_S_squared);
        }
        results[n_windows * (test + 3) + window] = T[test];
      }
    }
    
    n_mrk[window] = n_polymorphic_mrk;
  }
  
  for(int test = 0; test < n_tests; test++){
    for(int n = 1; n < nbr_chr; n++){
      free(Omega[test][n]);
    }
    free(Omega[test]);
  }
  free(Omega);
  for(int n = 1; n < nbr_chr; n++){
    free(omega_S[n]);
    free(omega_PI[n]);
    free(omega_L[n]);
    free(fs0[n]);
  }
  free(omega_S);
  free(omega_PI);
  free(omega_L);
  free(fs0);
  for(int i = 0; i < fs_size; i++){
    free(secMom0[i]);
  }
  free(secMom0);
  free(alpha);
  free(beta);
  free(T);
}
