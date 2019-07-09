#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP CALL_EHH(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CALL_EHHS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CALL_FURCATION(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CALL_INTEGRAL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CALL_SCAN_HH(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CALL_ASNEWICK(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"CALL_EHH",       (DL_FUNC) &CALL_EHH,        8},
    {"CALL_EHHS",      (DL_FUNC) &CALL_EHHS,       7},
    {"CALL_FURCATION", (DL_FUNC) &CALL_FURCATION,  7},
    {"CALL_INTEGRAL",  (DL_FUNC) &CALL_INTEGRAL,   8},
    {"CALL_SCAN_HH",   (DL_FUNC) &CALL_SCAN_HH,   16},
    {"CALL_ASNEWICK",  (DL_FUNC) &CALL_ASNEWICK,   6},
    {NULL, NULL, 0}
};

void R_init_rehh(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
