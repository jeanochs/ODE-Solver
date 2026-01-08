#include <stdlib.h>
#include <check.h>
#include <string.h>

#include "fe_section.h"

#define TOL 1e-3

// Directory where the input meshes are.
char input_mesh_dir[250];

START_TEST(linear_test_parser_1) {
	printf("Now testing with file `linear_mesh_1.in`.\n");
	// Passing
	struct Mesh m;

	FILE* linear_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "linear_mesh_1.in");

	linear_mesh = fopen(dir, "r");
	if (linear_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	// Parse the file, populate the mesh object
	int status = parse_input_file(linear_mesh, &m, LINEAR);	

	// Checks
	ck_assert_uint_eq(m.num_nodes, 4);
	ck_assert_uint_eq(m.num_elements, 3);

	// Element 1 (index 0)
	// First, check if the kind is linear
	ck_assert_int_eq(m.elements[0].kind, LINEAR);

	int n_id1 = m.connectivity_grid[0].node_list.L2.node_id[0];
	int n_id2 = m.connectivity_grid[0].node_list.L2.node_id[1];

	double n1 = m.elements[0].element.L2.node_coord[0];
	double n2 = m.elements[0].element.L2.node_coord[1];

	double N_0 = m.elements[0].element.L2.shape_func[0](0);
	double N_1 = m.elements[0].element.L2.shape_func[1](0);
	double N_0_D = m.elements[0].element.L2.shape_derv[0](0);
	double N_1_D = m.elements[0].element.L2.shape_derv[1](0);
	double J = GSL_FN_EVAL(&m.elements[0].element.L2.jacobian, 0);
	double x_phys = GSL_FN_EVAL(&m.elements[0].element.L2.iso_conversion, 0);

	ck_assert_int_eq(n_id1, 0);
	ck_assert_int_eq(n_id2, 1);
	ck_assert_double_eq_tol(n1, 0., TOL);
	ck_assert_double_eq_tol(n2, 1., TOL);
	ck_assert_double_eq_tol(N_0, 0.5, TOL);
	ck_assert_double_eq_tol(N_1, 0.5, TOL);
	ck_assert_double_eq_tol(N_0_D, -0.5, TOL);
	ck_assert_double_eq_tol(N_1_D, 0.5, TOL);
	ck_assert_double_eq_tol(J, 0.5, TOL);
	ck_assert_double_eq_tol(x_phys, 0.5, TOL);

	// Element 3 (index 2)
	// First, check if the kind is linear
	ck_assert_int_eq(m.elements[2].kind, LINEAR);

	n_id1 = m.connectivity_grid[2].node_list.L2.node_id[0];
	n_id2 = m.connectivity_grid[2].node_list.L2.node_id[1];

	n1 = m.elements[2].element.L2.node_coord[0];
	n2 = m.elements[2].element.L2.node_coord[1];

	N_0 = m.elements[2].element.L2.shape_func[0](0.0);
	N_1 = m.elements[2].element.L2.shape_func[1](0.0);
	N_0_D = m.elements[2].element.L2.shape_derv[0](0.0);
	N_1_D = m.elements[2].element.L2.shape_derv[1](0.0);
	J = GSL_FN_EVAL(&m.elements[2].element.L2.jacobian, 0);
	x_phys = GSL_FN_EVAL(&m.elements[2].element.L2.iso_conversion, 0);

	ck_assert_int_eq(n_id1, 2);
	ck_assert_int_eq(n_id2, 3);
	ck_assert_double_eq_tol(n1, 2., TOL);
	ck_assert_double_eq_tol(n2, 3., TOL);
	ck_assert_double_eq_tol(N_0, 0.5, TOL);
	ck_assert_double_eq_tol(N_1, 0.5, TOL);
	ck_assert_double_eq_tol(N_0_D, -0.5, TOL);
	ck_assert_double_eq_tol(N_1_D, 0.5, TOL);
	ck_assert_double_eq_tol(J, 0.5, TOL);
	ck_assert_double_eq_tol(x_phys, 2.5, TOL);

	// Parser should pass a code of sucess
	ck_assert_int_eq(status, EXIT_SUCCESS);

	if (!status) free_mesh_memory(&m);

}
END_TEST

START_TEST(linear_test_parser_2) {
	printf("Now testing with file `linear_mesh_2.in`.\n");
	// Passing
	struct Mesh m;

	FILE* linear_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "linear_mesh_2.in");

	linear_mesh = fopen(dir, "r");
	if (linear_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	// Parse the file, populate the mesh object
	int status = parse_input_file(linear_mesh, &m, LINEAR);	

	// Should've failed
	ck_assert_int_eq(status, EXIT_FAILURE);
	if (!status) free_mesh_memory(&m);
}
END_TEST

START_TEST(linear_test_parser_3) {
	printf("Now testing with file `linear_mesh_3.in`.\n");
	// Passing
	struct Mesh m;

	FILE* linear_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "linear_mesh_3.in");

	linear_mesh = fopen(dir, "r");
	if (linear_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	// Parse the file, populate the mesh object
	int status = parse_input_file(linear_mesh, &m, LINEAR);	

	// Should've failed
	ck_assert_int_eq(status, EXIT_FAILURE);
	if (!status) free_mesh_memory(&m);

}
END_TEST

START_TEST(linear_test_parser_4) {
	printf("Now testing with file `linear_mesh_4.in`.\n");
	// Passing
	struct Mesh m;

	FILE* linear_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "linear_mesh_4.in");

	linear_mesh = fopen(dir, "r");
	if (linear_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	// Parse the file, populate the mesh object
	int status = parse_input_file(linear_mesh, &m, LINEAR);	

	// Should've failed
	ck_assert_int_eq(status, EXIT_FAILURE);
	if (!status) free_mesh_memory(&m);

}
END_TEST

START_TEST(linear_test_parser_5) {
	printf("Now testing with file `linear_mesh_5.in`.\n");
	// Passing
	struct Mesh m;

	FILE* linear_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "linear_mesh_5.in");

	linear_mesh = fopen(dir, "r");
	if (linear_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}
	
	// Parse the file, populate the mesh object
	int status = parse_input_file(linear_mesh, &m, LINEAR);	

	// Should've failed
	ck_assert_int_eq(status, EXIT_FAILURE);
	if (!status) free_mesh_memory(&m);

}
END_TEST

START_TEST(quadratic_test_parser_1) {
	printf("Now testing with file `quadratic_mesh_1.in`.\n");
	// Passing
	struct Mesh m;

	FILE* quad_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "quadratic_mesh_1.in");

	quad_mesh = fopen(dir, "r");
	if (quad_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	// Parse the file, populate the mesh object
	int status = parse_input_file(quad_mesh, &m, QUAD);	

	// Checks
	ck_assert_uint_eq(m.num_nodes, 9);
	ck_assert_uint_eq(m.num_elements, 4);

	// Element 1 (index 0)
	// First, check if the kind is quadratic
	ck_assert_int_eq(m.elements[0].kind, QUAD);

	int n_id1 = m.connectivity_grid[0].node_list.L3.node_id[0];
	int n_id2 = m.connectivity_grid[0].node_list.L3.node_id[1];
	int n_id3 = m.connectivity_grid[0].node_list.L3.node_id[2];

	double n1 = m.elements[0].element.L3.node_coord[0];
	double n2 = m.elements[0].element.L3.node_coord[1];
	double n3 = m.elements[0].element.L3.node_coord[2];

	double N_0 = m.elements[0].element.L3.shape_func[0](0.1);
	double N_1 = m.elements[0].element.L3.shape_func[1](-0.3);
	double N_2 = m.elements[0].element.L3.shape_func[2](0.5);
	double N_0_D = m.elements[0].element.L3.shape_derv[0](0.73);
	double N_1_D = m.elements[0].element.L3.shape_derv[1](-0.34);
	double N_2_D = m.elements[0].element.L3.shape_derv[2](0.01);
	double J = GSL_FN_EVAL(&m.elements[0].element.L3.jacobian, -0.5);
	double x_phys = GSL_FN_EVAL(&m.elements[0].element.L3.iso_conversion, 0.4);

	ck_assert_int_eq(n_id1, 0);
	ck_assert_int_eq(n_id2, 1);
	ck_assert_int_eq(n_id3, 2);
	ck_assert_double_eq_tol(n1, 0, TOL);
	ck_assert_double_eq_tol(n2, 1.25, TOL);
	ck_assert_double_eq_tol(n3, 2.5, TOL);
	ck_assert_double_eq_tol(N_0, -0.045, TOL);
	ck_assert_double_eq_tol(N_1, 0.91, TOL);
	ck_assert_double_eq_tol(N_2, 0.375, TOL);
	ck_assert_double_eq_tol(N_0_D, 0.23, TOL);
	ck_assert_double_eq_tol(N_1_D, 0.68, TOL);
	ck_assert_double_eq_tol(N_2_D, 0.51, TOL);
	ck_assert_double_eq_tol(J, 1.25, TOL);
	ck_assert_double_eq_tol(x_phys, 1.75, TOL);

	// Element 3 (index 2)
	// First, check if the kind is quadratic
	ck_assert_int_eq(m.elements[2].kind, QUAD);

	n_id1 = m.connectivity_grid[2].node_list.L3.node_id[0];
	n_id2 = m.connectivity_grid[2].node_list.L3.node_id[1];
	n_id3 = m.connectivity_grid[2].node_list.L3.node_id[2];

	n1 = m.elements[2].element.L3.node_coord[0];
	n2 = m.elements[2].element.L3.node_coord[1];
	n3 = m.elements[2].element.L3.node_coord[2];

	N_0 = m.elements[2].element.L3.shape_func[0](-0.43);
	N_1 = m.elements[2].element.L3.shape_func[1](-0.99);
	N_2 = m.elements[2].element.L3.shape_func[2](0.48);
	N_0_D = m.elements[2].element.L3.shape_derv[0](0.13);
	N_1_D = m.elements[2].element.L3.shape_derv[1](-0.13);
	N_2_D = m.elements[2].element.L3.shape_derv[2](0.57);
	J = GSL_FN_EVAL(&m.elements[2].element.L3.jacobian, -0.07);
	x_phys = GSL_FN_EVAL(&m.elements[2].element.L3.iso_conversion, 0.1);

	ck_assert_int_eq(n_id1, 4);
	ck_assert_int_eq(n_id2, 5);
	ck_assert_int_eq(n_id3, 6);
	ck_assert_double_eq_tol(n1, 6., TOL);
	ck_assert_double_eq_tol(n2, 6.8, TOL);
	ck_assert_double_eq_tol(n3, 7., TOL);
	ck_assert_double_eq_tol(N_0, 0.30745, TOL);
	ck_assert_double_eq_tol(N_1, 0.02, TOL);
	ck_assert_double_eq_tol(N_2, 0.3552, TOL);
	ck_assert_double_eq_tol(N_0_D, -0.37, TOL);
	ck_assert_double_eq_tol(N_1_D, 0.26, TOL);
	ck_assert_double_eq_tol(N_2_D, 1.07, TOL);
	ck_assert_double_eq_tol(J, 0.542, TOL);
	ck_assert_double_eq_tol(x_phys, 6.847, TOL);

	// Parser should pass a code of success
	ck_assert_int_eq(status, EXIT_SUCCESS);

	if (!status) free_mesh_memory(&m);

}
END_TEST

START_TEST(quadratic_test_parser_2) {
	printf("Now testing with file `quadratic_mesh_2.in`.\n");
	// Passing
	struct Mesh m;

	FILE* quadratic_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "quadratic_mesh_2.in");

	quadratic_mesh = fopen(dir, "r");
	if (quadratic_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	// Parse the file, populate the mesh object
	int status = parse_input_file(quadratic_mesh, &m, LINEAR);	

	// Should've failed
	ck_assert_int_eq(status, EXIT_FAILURE);
	if (!status) free_mesh_memory(&m);

}
END_TEST

START_TEST(quadratic_test_parser_3) {
	printf("Now testing with file `quadratic_mesh_3.in`.\n");
	// Passing
	struct Mesh m;

	FILE* quadratic_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "quadratic_mesh_3.in");

	quadratic_mesh = fopen(dir, "r");
	if (quadratic_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	// Parse the file, populate the mesh object
	int status = parse_input_file(quadratic_mesh, &m, LINEAR);	

	// Should've failed
	ck_assert_int_eq(status, EXIT_FAILURE);
	if (!status) free_mesh_memory(&m);

}
END_TEST

START_TEST(quadratic_test_parser_4) {
	printf("Now testing with file `quadratic_mesh_4.in`.\n");
	// Passing
	struct Mesh m;

	FILE* quadratic_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "quadratic_mesh_4.in");

	quadratic_mesh = fopen(dir, "r");
	if (quadratic_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}

	// Parse the file, populate the mesh object
	int status = parse_input_file(quadratic_mesh, &m, LINEAR);	

	// Should've failed
	ck_assert_int_eq(status, EXIT_FAILURE);
	if (!status) free_mesh_memory(&m);

}
END_TEST

START_TEST(quadratic_test_parser_5) {
	printf("Now testing with file `quadratic_mesh_5.in`.\n");
	// Passing
	struct Mesh m;

	FILE* quadratic_mesh;
	char dir[250];
	memcpy(dir, input_mesh_dir, 250);

	strcat(dir, "quadratic_mesh_5.in");

	quadratic_mesh = fopen(dir, "r");
	if (quadratic_mesh == NULL) {
		printf("The file has not been found, or other error opening.\n");
		exit(1);
	}
	
	// Parse the file, populate the mesh object
	int status = parse_input_file(quadratic_mesh, &m, LINEAR);	

	// Should've failed
	ck_assert_int_eq(status, EXIT_FAILURE);
	if (!status) free_mesh_memory(&m);

}
END_TEST

// Setup and Execution
Suite* linear_parser_suite() {
	Suite *s;
	TCase *tc_core;

	s = suite_create("Parser Unit Tests");

	/* Core test case */
	tc_core = tcase_create("Core");

	tcase_add_test(tc_core, linear_test_parser_1);
	tcase_add_test(tc_core, linear_test_parser_2);
	tcase_add_test(tc_core, linear_test_parser_3);
	tcase_add_test(tc_core, linear_test_parser_4);
	tcase_add_test(tc_core, linear_test_parser_5);
	suite_add_tcase(s, tc_core);

	return s;
}

Suite* quadratic_parser_suite() {
	Suite* s;
	TCase* tc_core;

	s = suite_create("Quadratic Mesh Parser Tests");

	tc_core = tcase_create("Core");

	tcase_add_test(tc_core, quadratic_test_parser_1);
	tcase_add_test(tc_core, quadratic_test_parser_2);
	tcase_add_test(tc_core, quadratic_test_parser_3);
	tcase_add_test(tc_core, quadratic_test_parser_4);
	tcase_add_test(tc_core, quadratic_test_parser_5);
	suite_add_tcase(s, tc_core);

	return s;

}


int main() {
	// Assign directory here
	const char* dir_name = "./integration/test_meshes/";
	strcpy(input_mesh_dir, dir_name);

	int number_failed;
	Suite *s_parser_lin, *s_parser_quad;
	SRunner *sr_parser;

	// Set the suites here.
	s_parser_lin = linear_parser_suite();
	s_parser_quad = quadratic_parser_suite();

	// Setting the runners here
	sr_parser = srunner_create(s_parser_lin);
	srunner_add_suite(sr_parser, s_parser_quad);

	// Run the suites here.
	srunner_run_all(sr_parser, CK_NOFORK);

	number_failed = srunner_ntests_failed(sr_parser);

	srunner_free(sr_parser);

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}

	


