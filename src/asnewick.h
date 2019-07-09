#include <stdio.h>
#include <stdlib.h>

void asnewick(FILE* stream, const int* nbr_labels, const int* label_parent, const int* nbr_nodes,
		const int* node_parent, const double* node_pos, const double* xlim, char** const hap_name);
