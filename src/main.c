#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fe_section.h"

int main(int argc, char** argv) {
	if (argc != 9) {
		printf("Error in input list; check entry syntax.\nSyntax: solver.out [A] [B] [d1] [d2] [function field file] [start] [end] [number of elements]\nNote that all inputs are required.\n\n");
		return 1;
	}

	// Take out the numbers from the input
	double a_input = atof(argv[1]);
	double b_input = atof(argv[2]);
	double d1_input = atof(argv[3]);
	double d2_input = atof(argv[4]);

	FILE *function_file = fopen(argv[5], "r");
	if (function_file == NULL) {
		fprintf(stderr, "ERROR: Function field file %s does not exist.\n", argv[5]);
		return 1;
	}

	// Take in the start and end points of the domain, and the number of elements
	double start_point = atof(argv[6]);
	double end_point = atof(argv[7]);
	int num_elements = atoi(argv[8]);

	double ss = (end_point - start_point)/(num_elements - 1);

	FILE *input_file = fopen("input_mesh.in", "w");
	fprintf(input_file, "%d\n", num_elements);

	for (int i = 0; i < num_elements; i++) {
		double point = i*ss + start_point;
		fprintf(input_file, "%f\n", point);

	}

	fclose(input_file);

	// Open up the mesh file to be used
	FILE *mesh_file = fopen("input_mesh.in", "r");

	// Define the mesh, function field and solution object here
	struct Mesh m;
	struct ODE_Solution s;
	struct Function_Field ff;

	// Populate the field here
	if (input_function_field(&ff, function_file)) {
		fprintf(stderr, "There was an error.\n");
		return 1;
	}

	if (parse_input_file(mesh_file, &m, LINEAR)) return 1;

	int status = solve_ode_constant(&m, &s, a_input, b_input, d1_input, d2_input, &ff, false);
	if (status) return 1;

	// Output the results from the solution struct
	status = output_solution_data(&m, &s);
	if (status) return 1;

	fclose(function_file);
	fclose(mesh_file);
	free_mesh_memory(&m);
	free_solution_memory(&s);
	free_function_field(&ff);

	return 0;

}
