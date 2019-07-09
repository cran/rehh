#include <R.h>
#include "asnewick.h"
// #define DEBUG
#define precision 0 // precision (number digits after decimal point) of branch lengths

void asnewick_subtree(FILE* stream, const int* nbr_labels, const int* label_parent, const int* nbr_nodes,
		const int* node_parent, const double* node_pos, const double* xlim, char** const hap_name, const int node,
		const int is_last_sibling) {

	int last_sibling = 0;
	for (int i = node + 1; i < *nbr_nodes; i++) {
		if (node_parent[i] == node) {
			last_sibling = i;
		}
	}

	if (last_sibling > 0) {
		fprintf(stream, "(");
		for (int i = node + 1; i < *nbr_nodes; i++) {
			if (node_parent[i] == node) {
				asnewick_subtree(stream, nbr_labels, label_parent, nbr_nodes, node_parent, node_pos, xlim, hap_name, i,
						i == last_sibling);
			}
		}
	}
	
  // count and print up to two leave labels
	int node_size = 0;
	for (int i = 0; i < *nbr_labels; i++) {
		if (label_parent[i] == node) {
			if (node_size > 0) {
				if (node_size < 3) {
					fprintf(stream, "/");
				} else {
					fprintf(stream, "+");
				}
			}
			if (node_size < 3) {
				fprintf(stream, "%s", hap_name[i]);
			}
			node_size++;
		}
	}
	
#ifdef DEBUG
	// append external node number to leave label
	fprintf(stream, "_%i", node + 1);
#endif

	if (node_size == 1) {
	  fprintf(stream, ":%i", 0);
	} else if (node_size > 1) {
		double branch_length = *xlim - node_pos[node];
		if (branch_length < 0) {   // left tree
			branch_length *= -1;
		}
		fprintf(stream, ":%.*f", precision, branch_length);
	}

	if (!is_last_sibling) {
		fprintf(stream, ",");
	} else {
		fprintf(stream, ")");
		if (node != 0) {
			double branch_length = node_pos[node] - node_pos[node_parent[node]];
			if (branch_length < 0) {   // left tree
				branch_length *= -1;
			}
#ifdef DEBUG
			// include internal node number
			fprintf(stream, "%i", node_parent[node] + 1);
#endif
			fprintf(stream, ":%.*f", precision, branch_length);
		}
	}
	
}

void asnewick(FILE* stream, const int* nbr_labels, const int* label_parent, const int* nbr_nodes,
		const int* node_parent, const double* node_pos, const double* xlim, char** const hap_name) {
	fprintf(stream, "(");
	asnewick_subtree(stream, nbr_labels, label_parent, nbr_nodes, node_parent, node_pos, xlim, hap_name, 0, 1);
	fprintf(stream, ";");
}
