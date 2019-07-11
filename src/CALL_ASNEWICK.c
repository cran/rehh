#include <R.h>
#include <Rinternals.h>
#include "asnewick.h"

/**
 * Interface between R and C.
 * R objects are marked by a trailing underscore.
 */
SEXP CALL_ASNEWICK(SEXP tmp_file_name_, SEXP node_parent_, SEXP label_parent_, SEXP node_pos_, SEXP xlim_, SEXP hap_name_) {
	const int nbr_nodes = length(node_parent_);
	const int nbr_labels = length(label_parent_);
	const double xlim = asReal(xlim_);

	//get pointer to R data vectors
	const int* node_parent = INTEGER(node_parent_);
	const int* label_parent = INTEGER(label_parent_);
	const double* node_pos = REAL(node_pos_);

	FILE *stream;
	stream = fopen(CHAR(asChar(tmp_file_name_)), "w");
	
	if (stream == NULL) {
	  return ScalarLogical(0);
	}
	
	//copy string vector
	char** hap_name = (char**) malloc(sizeof(char*) * nbr_labels);
	for (int i = 0; i < nbr_labels; i++) {
	  hap_name[i] = (char*) malloc(strlen(CHAR(STRING_ELT(hap_name_, i))) + 1);
	  strcpy(hap_name[i], CHAR(STRING_ELT(hap_name_, i)));
	}
	
	asnewick(stream, &nbr_labels, label_parent, &nbr_nodes, node_parent, node_pos, &xlim, hap_name);

	fclose(stream);

	for (int i = 0; i < nbr_labels; i++) {
		free(hap_name[i]);
	}
	free(hap_name);

	return ScalarLogical(1);
}
