#include <stdlib.h>
#include <stdbool.h>

#ifdef DEBUG
#include <stdio.h>
#endif

#ifdef Conly
#ifndef MISSING_VALUE
#define MISSING_VALUE -1
#endif
#else
#include <Rinternals.h>
#ifndef MISSING_VALUE
#define MISSING_VALUE NA_INTEGER
#endif
#endif
