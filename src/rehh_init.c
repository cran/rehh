#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void CALL_EHH(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void CALL_EHHS(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void CALL_SCAN_HH(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"CALL_EHH",     (DL_FUNC) &CALL_EHH,     12},
    {"CALL_EHHS",    (DL_FUNC) &CALL_EHHS,    14},
    {"CALL_SCAN_HH", (DL_FUNC) &CALL_SCAN_HH, 13},
    {NULL, NULL, 0}
};

void R_init_rehh(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
