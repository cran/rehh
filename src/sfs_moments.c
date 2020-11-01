#include "sfs_moments.h"

/**
 * Harmonic numbers: 1+1/2+1/3+...+1/i
 */
double* getHarmonicNumbers(int n) {
  double* HarmonicNumbers = (double*) malloc(n*sizeof(double));
  
  if (!HarmonicNumbers)
    return 0;
  
  HarmonicNumbers[0] = 0;
  
  for (int i = 1; i < n; i++) {
    HarmonicNumbers[i] = HarmonicNumbers[i - 1] + 1.0 / i;
  }
  
  return HarmonicNumbers;
}

/**
 * The expected Xi spectrum under constant population size
 * cf. equation (1) of Fu 1995.
 */
double* getXi0(int n) {
  double* xi0 = (double*) malloc((n-1)*sizeof(double));
  for (int i = 1; i < n; i++) {
    xi0[i - 1] = 1. / i;
  }
  return xi0;
}

/**
 * the expected Eta spectrum under constant population size
 * cf. equations (6) and (7) of Fu 1995.
 */
double* getEta0(int n) {
  double* eta0 = (double*) malloc(n / 2 * sizeof(double));
  for (int i = 1; i <= n / 2; i++) {
    eta0[i - 1] = (1. / i + 1. / (n - i)) / (1. + KRONECKER(i, n - i));
  }
  return eta0;
}

/*
 * cf. equation (6) of Fu 1995.
 */
double getBeta(int i, double *HarmonicNumbers, int n) {
  double ai = HarmonicNumbers[i - 1], an = HarmonicNumbers[n - 1];
  double beta = 0;
  
  beta = 2.0 * n * (an + (1.0 / n) - ai) / ((n - i + 1.0) * (n - i))
    - 2.0 / (n - i);
  
  return beta;
}

/*
 * cf. equation (2) of Fu 1995.
 */
double getSigma_ii(int i, double *HarmonicNumbers, int n) {
  double sigma_ii = 0;
  double ai = HarmonicNumbers[i - 1], an = HarmonicNumbers[n - 1];
  
  if (2 * i < n) {
    sigma_ii = getBeta(i + 1, HarmonicNumbers, n);
  } else {
    if (2 * i == n) {
      sigma_ii = 2.0 * (an - ai) / (n - i) - 1.0 / (i * i);
    } else {
      sigma_ii = getBeta(i, HarmonicNumbers, n) - 1.0 / (i * i);
    }
    
  }
  
  return sigma_ii;
}

/*
 * cf. equation (3) of Fu 1995.
 */
double getSigma_ij(int i, int j, double *HarmonicNumbers, int n) {
  double sigma_ij = 0;
  
  if (i == j) {
    return getSigma_ii(i, HarmonicNumbers, n);
  }
  
  if (i < j) {
    int tmp = i;
    i = j;
    j = tmp;
  }
  
  double ai = HarmonicNumbers[i - 1], aj = HarmonicNumbers[j - 1], an =
    HarmonicNumbers[n - 1];
  
  if (i + j < n) {
    sigma_ij = (getBeta(i + 1, HarmonicNumbers, n)
                  - getBeta(i, HarmonicNumbers, n)) / 2.0;
  } else {
    if (i + j == n) {
      sigma_ij = ((an - ai) / (n - i) + (an - aj) / (n - j))
      - ((getBeta(i, HarmonicNumbers, n)
            + getBeta(j + 1, HarmonicNumbers, n)) / 2.0)
            - (1.0 / (i * j));
    } else {
      sigma_ij = ((getBeta(j, HarmonicNumbers, n)
                     - getBeta(j + 1, HarmonicNumbers, n)) / 2.0)
      - (1.0 / (i * j));
    }
  }
  return sigma_ij;
}


/**
 * sigma matrix = (sigma)ij
 */
double** getSigma(int n) {
  int i, j;
  double **sigma = (double**) malloc((n - 1) * sizeof(double*));
  double *HarmonicNumbers = getHarmonicNumbers(n);
  
  for (i = 0; i < n - 1; i++) {
    sigma[i] = malloc((n - 1) * sizeof(double));
  }
  
  for (i = 0; i < n - 1; i++) {
    for (j = i; j < n - 1; j++) {
      if (i == j) {
        sigma[i][i] = getSigma_ii(i + 1, HarmonicNumbers, n);
      } else {
        sigma[j][i] = sigma[i][j] = getSigma_ij(i + 1, j + 1,
                                                HarmonicNumbers, n);
      }
    }
  }
  free(HarmonicNumbers);
  return sigma;
}


/*
 * cf. equation (9) of Fu 1995
 */
double getRho_ii(int i, double *HarmonicSums, int n) {
  double rho_ii;
  
  rho_ii = getSigma_ii(i, HarmonicSums, n)
    + getSigma_ii(n - i, HarmonicSums, n)
    + 2 * getSigma_ij(i, n - i, HarmonicSums, n);
    rho_ii /= (1.0 + KRONECKER(i, n - i)) * (1.0 + KRONECKER(i, n - i));
    
    return rho_ii;
}

/*
 * cf. equation (9) of Fu 1995
 */
double getRho_ij(int i, int j, double *HarmonicSums, int n) {
  double rho_ij;
  
  rho_ij = getSigma_ij(i, j, HarmonicSums, n)
    + getSigma_ij(i, n - j, HarmonicSums, n)
    + getSigma_ij(n - i, j, HarmonicSums, n)
    + getSigma_ij(n - i, n - j, HarmonicSums, n);
    rho_ij /= ((1.0 + KRONECKER(i, n - i)) * (1.0 + KRONECKER(j, n - j)));
    
    return rho_ij;
}

/**
 * rho matrix = (rho)ij
 */
double** getRho(int n) {
  int i, j;
  double **rho = (double**) malloc(n / 2 * sizeof(double*));
  double *HarmonicNumbers = getHarmonicNumbers(n);
  
  for (i = 0; i < n / 2; i++) {
    rho[i] = (double*) malloc(n / 2 * sizeof(double));
  }
  
  for (i = 0; i < n / 2; i++) {
    for (j = i; j < n / 2; j++) {
      if (i == j) {
        rho[i][i] = getRho_ii(i + 1, HarmonicNumbers, n);
      } else {
        rho[j][i] = rho[i][j] = getRho_ij(i + 1, j + 1, HarmonicNumbers,
                                          n);
      }
    }
  }
  return rho;
}

/*
 * first coefficient in nominator of test statistic
 */
double getWeightedFirstMoment(double* weight, double* fs0, int fs_size) {
  int i;
  double alpha = 0;
  
  for (i = 0; i < fs_size; i++) {
    alpha += weight[i] * weight[i] / fs0[i];
  }
  
  return alpha;
}

/*
 * second coefficient in nominator of test statistic
 */
double getWeightedSecondMoment(double* weight, double* fs0, double** secMom0,
                               int fs_size) {
  int i, j;
  double beta = 0;
  
  for (i = 0; i < fs_size; i++) {
    for (j = 0; j < fs_size; j++) {
      beta += (weight[i] / fs0[i]) * secMom0[i][j] * (weight[j] / fs0[j]);
    }
  }
  
  return beta;
}

/*
 * compute Omega for folded theta estimators
 * cf. Table 1 of Achaz (2009)
 */
double* getOmegaStar(int n, char type) {
  double* omega = (double*) malloc((n / 2) * sizeof(double));

  switch(type){
  case 'S':  
    for (int i = 1; i <= n / 2; i++) {
      omega[i - 1] = n / (double) (i * (n - i) * (1 + KRONECKER(i, n - i))) ;
    }
    break;
  case 'P': 
    for (int i = 1; i <= n / 2; i++) {
      omega[i - 1] = n / (double) (1 + KRONECKER(i, n - i));
    }
    break;
  case 'L':
    for (int i = 1; i <= n / 2; i++) {
      omega[i - 1] = 1. / (double) (1 + KRONECKER(i, n - i));
    }
    break;
  default: return NULL;
  }

  double sum = 0.;  
  for (int i = 0; i < n / 2; i++) {
    sum += omega[i];
  }
  
  for (int i = 0; i < n / 2; i++) {
    omega[i] /= sum;
  }
  
  return omega;
}

/*
 * compute omega for unfolded theta estimators
 * cf. Table 1 of Achaz (2009)
 */
double* getOmega(int n, char type) {
  double* omega = (double*) malloc((n - 1) * sizeof(double));

  switch(type){
    case 'S':  
      for (int i = 1; i < n; i++) {
        omega[i - 1] = 1.0 / i;
      }
      break;
    case 'P': 
      for (int i = 1; i < n; i++) {
        omega[i - 1] = (double) (n - i);
      }
      break;
    case 'L':
      for (int i = 1; i < n; i++) {
        omega[i - 1] = 1.0;
      }
      break;
    default: return NULL;
  }
  
  double sum = 0.;
  for (int i = 0; i < n - 1; i++) {
    sum += omega[i];
  }
  
  for (int i = 0; i < n - 1; i++) {
    omega[i] /= sum;
  }
  
  return omega;
}
