#include <R.h>
#include "definitions.h"

#define KRONECKER( A, B )  (((A)==(B))?1:0)

double* getXi0(int n);
double* getEta0(int n);
double** getSigma(int n);
double** getRho(int n);
double getWeightedFirstMoment(double* weight, double* fs0, int fs_size);
double getWeightedSecondMoment(double* weight, double* fs0, double** secMom0, int fs_size);
double* getOmega(int n, char test);
double* getOmegaStar(int n, char test);
