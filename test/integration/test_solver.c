#include <stdlib.h>
#include <check.h>
#include <string.h>

#include "fe_section.h"

#define TOL 1e-3

struct Function_Field *field = NULL;

// Directory where the input meshes are.
char input_mesh_dir[250];

double driving_func(double x) {
	return x*x + x + 3;

}

static void setup_function_field() {
	// This would normally be handled by a dedicated function, but this hasn't been written yet.
	// This setup will probably stay, though
	// Create the Function Field object
	struct Function_Field *f_field_temp = malloc(sizeof(struct Function_Field));

	create_function_field(f_field_temp, 0, 15, 2001, driving_func);

	field = f_field_temp;

}

static void teardown_function_field() {
	free_function_field(field);
	free(field);
	field = NULL;

}

START_TEST(L2_solver) {
	printf("Solving ODE with the linear mesh.\n");
	// Read in the linear mesh
	struct Mesh m;
	struct ODE_Solution sol;

	FILE* linear_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "linear_mesh.in");

	linear_mesh = fopen(dir, "r");
	if (linear_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	int status = parse_input_file(linear_mesh, &m, LINEAR);
	if (status) {
		printf("Error reading in the mesh file.\n");
		exit(1);
	}

	status = solve_ode_constant(&m, &sol, 4., 4., 0, 5, field, true);
	if (status) {
		printf("Error in solving mesh. Check.\n");
		exit(1);
	}

	char dir1[250], dir2[250], dir3[250];
	memcpy(dir1, input_mesh_dir, 250);
	memcpy(dir2, input_mesh_dir, 250);
	memcpy(dir3, input_mesh_dir, 250);

	strcat(dir1, "linear_matrix_global.output");
	strcat(dir2, "linear_vector_global.output");
	strcat(dir3, "linear_solution_vector.output");

	FILE* reference_matrix = fopen(dir1, "r");
	FILE* reference_vector = fopen(dir2, "r");
	FILE* reference_solution = fopen(dir3, "r");

	// Check the global matrices and vectors
	int size = m.num_nodes;
	gsl_matrix* K_global_ref = gsl_matrix_alloc(size, size);
	gsl_vector* F_global_ref = gsl_vector_alloc(size);
	gsl_vector* solution_ref = gsl_vector_alloc(size);

	status = gsl_matrix_fread(reference_matrix, K_global_ref);
	if (status) {
		printf("Error in reading matrix in. Check.\n");
		gsl_matrix_free(K_global_ref);
		gsl_vector_free(F_global_ref);
		gsl_vector_free(solution_ref);
		exit(1);
	}

	status = gsl_vector_fread(reference_vector, F_global_ref);
	if (status) {
		printf("Error in reading vector in. Check.\n");
		gsl_matrix_free(K_global_ref);
		gsl_vector_free(F_global_ref);
		gsl_vector_free(solution_ref);
		exit(1);
	}

	status = gsl_vector_fread(reference_solution, solution_ref);
	if (status) {
		printf("Error in reading vector in. Check.\n");
		gsl_matrix_free(K_global_ref);
		gsl_vector_free(F_global_ref);
		gsl_vector_free(solution_ref);
		exit(1);
	}

	fclose(reference_matrix);
	fclose(reference_vector);
	fclose(reference_solution);


	// Now, compare the matrices
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			printf("Matrix Index %d, %d\n\n", i, j);
			ck_assert_double_eq_tol(
				gsl_matrix_get(sol.coeff_matrix_global, i, j),
				gsl_matrix_get(K_global_ref, i, j),
				TOL
			);

		}
	}

	for (int i = 0; i < size; i++) {
		printf("Vector Index %d\n\n", i);
		ck_assert_double_eq_tol(
			gsl_vector_get(sol.const_vector_global, i),
			gsl_vector_get(F_global_ref, i),
			TOL
		);

	}

	// Check the solution
	for (int i = 0; i < size; i++) {
		printf("Solution Vector Index %d\n\n", i);
		ck_assert_double_eq_tol(
			gsl_vector_get(sol.solution_coeff, i),
			gsl_vector_get(solution_ref, i),
			TOL
		);

	}

	// Need to move to a setup/teardown fixture
	gsl_matrix_free(K_global_ref);
	gsl_vector_free(F_global_ref);
	gsl_vector_free(solution_ref);
	free_mesh_memory(&m);
	free_solution_memory(&sol);

}
END_TEST

START_TEST(L3_solver) {
	printf("Solving ODE with the quadratic mesh.\n");
	// Read in the linear mesh
	struct Mesh m;
	struct ODE_Solution sol;

	FILE* quad_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "quadratic_mesh.in");

	quad_mesh = fopen(dir, "r");
	if (quad_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	int status = parse_input_file(quad_mesh, &m, QUAD);
	if (status) {
		printf("Error reading in the mesh file.\n");
		exit(1);
	}

	status = solve_ode_constant(&m, &sol, 4., 4., 0, 5, field, true);
	if (status) {
		printf("Error in solving mesh. Check.\n");
		exit(1);
	}

	// Read in the reference files
	char dir1[250], dir2[250], dir3[250];
	memcpy(dir1, input_mesh_dir, 250);
	memcpy(dir2, input_mesh_dir, 250);
	memcpy(dir3, input_mesh_dir, 250);

	strcat(dir1, "quad_matrix_global.output");
	strcat(dir2, "quad_vector_global.output");
	strcat(dir3, "quad_solution_vector.output");

	FILE* reference_matrix = fopen(dir1, "r");
	FILE* reference_vector = fopen(dir2, "r");
	FILE* reference_solution = fopen(dir3, "r");

	// Check the global matrices and vectors
	int size = m.num_nodes;
	gsl_matrix* K_global_ref = gsl_matrix_alloc(size, size);
	gsl_vector* F_global_ref = gsl_vector_alloc(size);
	gsl_vector* solution_ref = gsl_vector_alloc(size);

	status = gsl_matrix_fread(reference_matrix, K_global_ref);
	if (status) {
		printf("Error in reading matrix in. Check.\n");
		gsl_matrix_free(K_global_ref);
		gsl_vector_free(F_global_ref);
		gsl_vector_free(solution_ref);
		exit(1);
	}

	status = gsl_vector_fread(reference_vector, F_global_ref);
	if (status) {
		printf("Error in reading vector in. Check.\n");
		gsl_matrix_free(K_global_ref);
		gsl_vector_free(F_global_ref);
		gsl_vector_free(solution_ref);
		exit(1);
	}

	status = gsl_vector_fread(reference_solution, solution_ref);
	if (status) {
		printf("Error in reading vector in. Check.\n");
		gsl_matrix_free(K_global_ref);
		gsl_vector_free(F_global_ref);
		gsl_vector_free(solution_ref);
		exit(1);
	}

	fclose(reference_matrix);
	fclose(reference_vector);
	fclose(reference_solution);


	// Now, compare the matrices
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			printf("Matrix Index %d, %d\n\n", i, j);
			ck_assert_double_eq_tol(
				gsl_matrix_get(sol.coeff_matrix_global, i, j),
				gsl_matrix_get(K_global_ref, i, j),
				TOL
			);

		}
	}

	for (int i = 0; i < size; i++) {
		printf("Vector Index %d\n\n", i);
		ck_assert_double_eq_tol(
			gsl_vector_get(sol.const_vector_global, i),
			gsl_vector_get(F_global_ref, i),
			TOL
		);

	}

	// Check the solution
	for (int i = 0; i < size; i++) {
		printf("Solution Vector Index %d\n\n", i);
		ck_assert_double_eq_tol(
			gsl_vector_get(sol.solution_coeff, i),
			gsl_vector_get(solution_ref, i),
			TOL
		);

	}

	// Need to move to a setup/teardown fixture
	gsl_matrix_free(K_global_ref);
	gsl_vector_free(F_global_ref);
	gsl_vector_free(solution_ref);
	free_mesh_memory(&m);
	free_solution_memory(&sol);

}
END_TEST

Suite* solver_suite() {
	Suite *s;
	TCase *tc_linear, *tc_quad;

	s = suite_create("Solver Integration Tests");

	// LInear Mesh Core Test Case
	tc_linear = tcase_create("Linear Mesh Case");
	tcase_add_checked_fixture(
		tc_linear,
		setup_function_field,
		teardown_function_field
	);

	tcase_add_test(tc_linear, L2_solver);
	suite_add_tcase(s, tc_linear);

	tc_quad = tcase_create("Quadratic Mesh Case");
	tcase_add_checked_fixture(
		tc_quad,
		setup_function_field,
		teardown_function_field
	);
	tcase_add_test(tc_quad, L3_solver);
	suite_add_tcase(s, tc_quad);

	return s;

}

int main() {
	// Assign directory here
	const char *dir_name = "./integration/test_reference_array/";
	strcpy(input_mesh_dir, dir_name);

	int number_failed;
	Suite *s_solver;
	SRunner *sr_solver;

	// Set the suites here.
	s_solver = solver_suite();

	// Setting the runners here
	sr_solver = srunner_create(s_solver);

	// Run the suites here.
	srunner_set_fork_status(sr_solver, CK_NOFORK);
	srunner_run_all(sr_solver, CK_NORMAL);

	number_failed = srunner_ntests_failed(sr_solver);

	srunner_free(sr_solver);

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
