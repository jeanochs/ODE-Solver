#include "fe_section.h"

#include "shape_functions.c"
#include "composition_functions.c"

// Function to free the memory allocated by the Mesh object
void free_mesh_memory(struct Mesh* input_mesh) {
	free(input_mesh->connectivity_grid);
	free(input_mesh->elements);

}

struct Element_Linear create_element_L2(double node1, double node2) {
	// Create element object
	struct Element_Linear e;
	e.kind = LINEAR;

	e.element.L2.node_coord[0] = node1;
	e.element.L2.node_coord[1] = node2;

	// Populate the shape functions and their derivatives
	e.element.L2.shape_func[0] = L2_N0;
	e.element.L2.shape_func[1] = L2_N1;
	e.element.L2.shape_derv[0] = L2_N0_D;
	e.element.L2.shape_derv[1] = L2_N1_D;

	// Populate the Jacobian equation and the isoparametric-physical relation
	// Jacobian assembly into GSL function
	struct L2_N_P p = {node1, node2};
	gsl_function j;
	j.function = &J_L2;
	j.params = &p;

	e.element.L2.jacobian = j;

	// Isoparametric assembly into GSL function
	struct Iso_Phy_Funcs p1 = {
		2,
		e.element.L2.node_coord,
		e.element.L2.shape_func
	};

	gsl_function iso;
	iso.function = &isoparametric_physical_conversion;
	iso.params = &p1;
	
	e.element.L2.iso_conversion = iso;

	return e;

}

struct Element_Linear create_element_L3(double node1, double node2, double node3) {
	// Create element object
	struct Element_Linear e;
	e.kind = LINEAR;

	e.element.L3.node_coord[0] = node1;
	e.element.L3.node_coord[1] = node2;
	e.element.L3.node_coord[2] = node3;

	// Populate the shape functions and their derivatives
	e.element.L3.shape_func[0] = L3_N0;
	e.element.L3.shape_func[1] = L3_N1;
	e.element.L3.shape_func[2] = L3_N2;
	e.element.L3.shape_derv[0] = L3_N0_D;
	e.element.L3.shape_derv[1] = L3_N1_D;
	e.element.L3.shape_derv[2] = L3_N2_D;

	// Populate the Jacobian equation and the isoparametric-physical relation
	// Jacobian assembly into GSL function
	struct L3_N_P p = {node1, node2, node3};
	gsl_function j;
	j.function = &J_L3;
	j.params = &p;

	e.element.L3.jacobian = j;

	// Isoparametric assembly into GSL function
	struct Iso_Phy_Funcs p1 = {
		3,
		e.element.L3.node_coord,
		e.element.L3.shape_func
	};

	gsl_function iso;
	iso.function = &isoparametric_physical_conversion;
	iso.params = &p1;
	
	e.element.L3.iso_conversion = iso;

	return e;

}

int parse_input_file(FILE* input_stream, struct Mesh* mesh_object, Element_2D_Type mesh_kind) {
	// Read each line from the input file.
	// Note that first line should be parsed as an unsinged integer; it gives node count
	char buffer[100];

	if (fgets(buffer, 100, input_stream) == NULL) {
		printf("Error opening file, or file is empty or invalid. Please check.\nAborting...");
		return 1;
	}

	// Convert to an integer
	int num_nodes = strtol(buffer, NULL, 10);
	mesh_object->num_nodes = num_nodes;

	// Allocate an array
	// TODO: Please rename
	double* node_coors = (double*) malloc(num_nodes*sizeof(double));
	if (node_coors == NULL) {
		printf("Error in allocating node coordinate array of length %d.\nAborting...", num_nodes);
		return 1;
	}

	// Parse the remaining as floats and produce the node array
	int counter = 0;
	double prev_node_coord; // Nodes must be organized in ascending sequential order. 
	while (fgets(buffer, 100, input_stream) != NULL) {
		double coord = strtof(buffer, NULL);
		node_coors[counter++] = coord;

		// Don't compare the first node coordiante to anything; simply set it as the first previoud node coordinate.
		if (counter == 1) {
			prev_node_coord = coord;
		}
		else {
			if (coord <= prev_node_coord) {
					printf("CRITICAL ERROR: Mesh file is malformed; coordinate of node %d (%f) is equal or less than the previous node %d (%f).\n", counter, coord, counter - 1, prev_node_coord);
					free(node_coors);
					return 1;
			}
			prev_node_coord = coord;
		}

	}

	if (counter != num_nodes) {
		printf("CRITICAL ERROR: Mesh file is malformed; number of nodes reported (%d) is not equal to the number of nodes scanned (%d).\n", num_nodes, counter);
		free(node_coors);
		return 1;
	}

	// Now, determine the number of elements and produce them
	switch (mesh_kind) {
		case LINEAR: {
			mesh_object->num_elements = num_nodes - 1;

			// Allocate the element array and element connectivity grid
			struct Element_Linear* element_array = (struct Element_Linear*) malloc(mesh_object->num_elements*sizeof(struct Element_Linear));
			struct Element_Conn* conn_array = (struct Element_Conn*) malloc(mesh_object->num_elements*sizeof(struct Element_Conn));

			if (element_array == NULL || conn_array == NULL) {
				printf("Error allocating array for elements of length %d. Please check.\nAborting...", num_nodes - 1);
				return 1;
			}

			for (int i = 1; i < num_nodes; i++) {
				struct Element_Linear ele = create_element_L2(node_coors[i - 1], node_coors[i]);

				// Add the ready element to the element array
				element_array[i - 1] = ele;

				// Populate the connectivity grid.
				struct Element_Conn conn_entry;
				conn_entry.kind = LINEAR;
				conn_entry.node_list.L2.node_id[0] = i - 1;
				conn_entry.node_list.L2.node_id[1] = i;
				conn_array[i - 1] = conn_entry;
			}

			// Finally, add the completed element array to the mesh object
			mesh_object->elements = element_array;
			mesh_object->connectivity_grid = conn_array;
			break;
		}

		case QUAD: {
			mesh_object->num_elements = (num_nodes - 1)/2;

			// Allocate the element array
			struct Element_Linear* element_array = (struct Element_Linear*) malloc(mesh_object->num_elements*sizeof(struct Element_Linear));
			struct Element_Conn* conn_array = (struct Element_Conn*) malloc(mesh_object->num_elements*sizeof(struct Element_Conn));

			if (element_array == NULL || conn_array == NULL) {
				printf("Error allocating array for elements of length %d. Please check.\nAborting...", num_nodes - 1);
				return 1;
			}

			for (int i = 2; i < num_nodes; i++) {
				struct Element_Linear ele = create_element_L3(node_coors[i - 2], node_coors[i - 1],  node_coors[i]);


				// Add the ready element to the element array
				element_array[i - 1] = ele;
				
				// Populate the connectivity grid.
				struct Element_Conn conn_entry;
				conn_entry.kind = QUAD;
				conn_entry.node_list.L3.node_id[0] = i - 2;
				conn_entry.node_list.L3.node_id[1] = i - 1;
				conn_entry.node_list.L3.node_id[2] = i;
				conn_array[i - 1] = conn_entry;
			}

			// Finally, add the completed element array to the mesh object
			mesh_object->elements = element_array;
			mesh_object->connectivity_grid = conn_array;
			break;
		}
	}

	return 0;

}

gsl_vector* output_constant_vector(struct Element_Linear* element, double (*driving_func) (double)) {
	// Check element type. This will determine iteration length
	switch (element->kind) {
		case LINEAR: { // Two nodes
					   // Create the vector
						 gsl_vector* v = gsl_vector_alloc(2);

						 for (int n = 0; n < 2; n++) {
							 gsl_function constant_vector;

							 struct Constant_Vector_Funcs p = {
								 driving_func,
								 element->element.L2.shape_func[n],
								 &element->element.L2.iso_conversion
							 };

							 constant_vector.function = &constant_vector_composition;
							 constant_vector.params = &p;

							 // The crux of this operation.
							 // Perform Gauss quadrature and add to the n-th entry
							 // Make the workspace first
							 gsl_integration_fixed_workspace* w = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre,
									 1,
									 -1,
									 1,
									 0, // Ignored
									 0 // Ignored
									 );
							 // Finally integrate.
							 double f_i1;
							 gsl_integration_fixed(&constant_vector, &f_i1, w);
							 gsl_vector_set(v, n, f_i1);
							 gsl_integration_fixed_free(w);
						 }

						 return v;
						 break;
					 }
		case QUAD: { // Three nodes
					 // Create the vector
					   gsl_vector* v = gsl_vector_alloc(3);

					   for (int n = 0; n < 3; n++) {
						   gsl_function constant_vector;

						   struct Constant_Vector_Funcs p = {
							   driving_func,
							   element->element.L3.shape_func[n],
							   &element->element.L3.iso_conversion
						   };

						   constant_vector.function = &constant_vector_composition;
						   constant_vector.params = &p;

						   // The crux of this operation.
						   // Perform Gauss quadrature and add to the n-th entry
						   // Make the workspace first
						   gsl_integration_fixed_workspace* w = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre,
								   1,
								   -1,
								   1,
								   0, // Ignored
								   0 // Ignored
								   );
						   // Finally integrate.
						   double f_i1;
						   gsl_integration_fixed(&constant_vector, &f_i1, w);
						   gsl_vector_set(v, n, f_i1);
						   gsl_integration_fixed_free(w);
					   }

					   return v;
					   break;
				   }
	}

}

gsl_matrix* output_coefficient_matrix(struct Element_Linear* element, double a, double b) {
	switch (element->kind) {
		case LINEAR: {
						 gsl_matrix* m = gsl_matrix_alloc(2, 2);

						 for (int i = 0; i < 2; i++) {
							 for (int j = 0; j < 2; j++) {
								 gsl_function coefficient_matrix;

								 struct Coefficient_Matrix_Funcs p = {
									 a, b,
									 element->element.L2.shape_func[i],
									 element->element.L2.shape_func[j],
									 element->element.L2.shape_derv[i],
									 element->element.L2.shape_derv[j],
									 &element->element.L2.jacobian
								 };

								 coefficient_matrix.function = &coefficient_matrix_composition;
								 coefficient_matrix.params = &p;

								 // The crux of this operation.
								 // Perform Gauss quadrature and add to the n-th entry
								 // Make the workspace first
								 gsl_integration_fixed_workspace* w = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre,
										 1,
										 -1,
										 1,
										 0, // Ignored
										 0 // Ignored
										 );
								 // Finally integrate.
								 double k_ij;
								 gsl_integration_fixed(&coefficient_matrix, &k_ij, w);
								 gsl_matrix_set(m, i, j, k_ij);
								 gsl_integration_fixed_free(w);
							 }
						 }

						 return m;
						 break;
					 }
		case QUAD: {
						 gsl_matrix* m = gsl_matrix_alloc(3, 3);

						 for (int i = 0; i < 3; i++) {
							 for (int j = 0; j < 3; j++) {
								 gsl_function coefficient_matrix;

								 struct Coefficient_Matrix_Funcs p = {
									 a, b,
									 element->element.L3.shape_func[i],
									 element->element.L3.shape_func[j],
									 element->element.L3.shape_derv[i],
									 element->element.L3.shape_derv[j],
									 &element->element.L3.jacobian
								 };

								 coefficient_matrix.function = &coefficient_matrix_composition;
								 coefficient_matrix.params = &p;

								 // The crux of this operation.
								 // Perform Gauss quadrature and add to the n-th entry
								 // Make the workspace first
								 gsl_integration_fixed_workspace* w = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre,
										 1,
										 -1,
										 1,
										 0, // Ignored
										 0 // Ignored
										 );
								 // Finally integrate.
								 double k_ij;
								 gsl_integration_fixed(&coefficient_matrix, &k_ij, w);
								 gsl_matrix_set(m, i, j, k_ij);
								 gsl_integration_fixed_free(w);
							 }
						 }

						 return m;
						 break;
					 }
	}

}

int solve_ode_constant(struct Mesh* input_mesh, struct ODE_Solution* solution, double a, double b, double (*func) (double)) {
	// First, check if the input mesh has valid node and element arrays
	if (input_mesh->connectivity_grid == NULL || input_mesh->elements == NULL) {
		printf("ERROR: Provided mesh is not properly loaded with element and node information.\nPlease ensure that the `parse_input_file` function has been called to populate the object, or check for other errors.\n");
		return 1;
	}

	// Initialize the constant vector and coefficient matrix with zeros
	// Coeffcient matrix size : (num_nodes, num_nodes)
	// Constant vector size : (num_nodes, 1)
	gsl_matrix* K_coeff = gsl_matrix_alloc(input_mesh->num_nodes, input_mesh->num_nodes);
	gsl_vector* F_const = gsl_vector_alloc(input_mesh->num_nodes);

	gsl_vector_set_zero(F_const);
	gsl_matrix_set_zero(K_coeff);

	// Now, iterate through the elements and solve for the local coefficient matrix and constant vectors 
	for (int e = 0; e < input_mesh->num_elements; e++) {
		// Calculate the constant vector
		gsl_vector* constant_local = output_constant_vector(&input_mesh->elements[e], func);

		// Calculate the coefficient matrix
		gsl_matrix* coefficient_local = output_coefficient_matrix(&input_mesh->elements[e], a, b);

		// Add to their respective entries in the global matrices
		int starting_point, size;
		switch (input_mesh->elements[e].kind) {
			case LINEAR: 
				starting_point = input_mesh->connectivity_grid[e].node_list.L2.node_id[0];
				size = 2;
				break;
			case QUAD: 
				starting_point = input_mesh->connectivity_grid[e].node_list.L2.node_id[0];
				size = 3;
				break;
		}

		gsl_matrix_view submatrix = gsl_matrix_submatrix(K_coeff, starting_point, starting_point, size, size);
		gsl_vector_view subvector = gsl_vector_subvector(F_const, starting_point, size);

		// Add the local matrices to the views
		gsl_matrix_add(&submatrix.matrix, coefficient_local);
		gsl_vector_add(&subvector.vector, constant_local);

	}

	// With populated matrix and vector, solve the linear equation [K][y] = [F]
	// Uisng the LU decomp solver
	gsl_permutation* perm = gsl_permutation_alloc(input_mesh->num_nodes);
	gsl_permutation_init(perm);
	int signum;
	gsl_linalg_LU_decomp(K_coeff, perm, &signum);

	// Initialize solution vector here
	gsl_vector* variable_vector = gsl_vector_alloc(input_mesh->num_nodes);

	// Solve.
	gsl_linalg_LU_solve(K_coeff, perm, F_const, variable_vector);

	// Pass the now solved variable vector to the Solution output.
	solution->solution_coeff = variable_vector;

	// Done.
	
	return 0;

}


