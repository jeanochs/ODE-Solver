#include <stdlib.h>
#include <check.h>
#include <string.h>

#include "fe_section.h"

#define TOL 1e-5

// Directory where the input meshes are.
char input_mesh_dir[250];

START_TEST(test_parser_1) {
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

	int n1 = m.elements[0].element.L2.node_id[0];
	int n2 = m.elements[0].element.L2.node_id[1];

	double N_0 = m.elements[0].element.L2.shape_func[0](0);
	double N_1 = m.elements[0].element.L2.shape_func[1](0);
	double N_0_D = m.elements[0].element.L2.shape_derv[0](0);
	double N_1_D = m.elements[0].element.L2.shape_derv[1](0);

	ck_assert_int_eq(n1, 0);
	ck_assert_int_eq(n2, 1);
	ck_assert_double_eq_tol(N_0, 0.5, TOL);
	ck_assert_double_eq_tol(N_1, 0.5, TOL);
	ck_assert_double_eq_tol(N_0_D, -0.5, TOL);
	ck_assert_double_eq_tol(N_1_D, 0.5, TOL);

	// Element 3 (index 2)
	// First, check if the kind is linear
	ck_assert_int_eq(m.elements[2].kind, LINEAR);

	n1 = m.elements[2].element.L2.node_id[0];
	n2 = m.elements[2].element.L2.node_id[1];

	N_0 = m.elements[2].element.L2.shape_func[0](0.0);
	N_1 = m.elements[2].element.L2.shape_func[1](0.0);
	N_0_D = m.elements[2].element.L2.shape_derv[0](0.0);
	N_1_D = m.elements[2].element.L2.shape_derv[1](0.0);

	ck_assert_int_eq(n1, 2);
	ck_assert_int_eq(n2, 3);
	ck_assert_double_eq_tol(N_0, 0.5, TOL);
	ck_assert_double_eq_tol(N_1, 0.5, TOL);
	ck_assert_double_eq_tol(N_0_D, -0.5, TOL);
	ck_assert_double_eq_tol(N_1_D, 0.5, TOL);

	// Parser should pass a code of sucess
	ck_assert_int_eq(status, EXIT_SUCCESS);

	if (!status) free_mesh_memory(&m);

}
END_TEST

START_TEST(test_parser_2) {
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

START_TEST(test_parser_3) {
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

START_TEST(test_parser_4) {
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

START_TEST(test_parser_5) {
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


// Setup and Execution
Suite *parser_suite() {
	Suite *s;
	TCase *tc_core;

	s = suite_create("Parser Unit Tests");

	/* Core test case */
	tc_core = tcase_create("Core");

	tcase_add_test(tc_core, test_parser_1);
	tcase_add_test(tc_core, test_parser_2);
	tcase_add_test(tc_core, test_parser_3);
	tcase_add_test(tc_core, test_parser_4);
	tcase_add_test(tc_core, test_parser_5);
	suite_add_tcase(s, tc_core);

	return s;
}


int main() {
	// Assign directory here
	const char* dir_name = "./unit/test_meshes/";
	strcpy(input_mesh_dir, dir_name);

	int number_failed;
	Suite *s_parser;
	SRunner *sr_parser;

	// Set the suites here.
	s_parser = parser_suite();

	// Setting the runners here
	sr_parser = srunner_create(s_parser);

	// Run the suites here.
	srunner_run_all(sr_parser, CK_NORMAL);

	number_failed = srunner_ntests_failed(sr_parser);

	srunner_free(sr_parser);

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}

	


