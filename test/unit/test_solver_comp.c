#include <stdlib.h>
#include <check.h>
#include <string.h>

#include "fe_section.h"

#define TOL 1e-3

double driving_func(double x) {
	return (x*x) + 2*x  - 1;

}

// Constants
double a = 4.0f;
double b = 4.0f;

// Predefined matrices
// L2 Element Matrices and Vectors
double L2_K[2][2] = {
	{7/6, 17/6},
	{-7/6, 31/6}
};

double L2_F[2] = {23/3, 11};

// L3 Element Matrices and Vectors
double L3_K[3][3] = {
	{0.2099, 2.2006, -0.9438},
	{-3.133, 8.2577, 2.8750},
	{0.3895, -2.4583, 4.6021}
};

double L3_F[3] = {0.4423, 13.7783, 8.2794};


START_TEST(L2_element_comp) {
	// Create the linear element
	struct Element_Linear ele = create_element_L2(1, 3);

	gsl_vector* v = output_constant_vector(&ele, driving_func);
	gsl_matrix* m = output_coefficient_matrix(&ele, a, b);

	// Loading the matrices and vectors
	gsl_vector* v_check = gsl_vector_alloc(2);
	gsl_matrix* m_check = gsl_matrix_alloc(2, 2);
	for (int i = 0; i < 2; i++) {
		gsl_vector_set(v_check, i, L2_F[i]);
	}

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			gsl_matrix_set(m_check, i, j, L2_K[i][j]);
		}
	}

	ck_assert_mem_eq(v->data, v_check->data, 2);
	ck_assert_mem_eq(m->data, m_check->data, 4);

} 
END_TEST

START_TEST(L3_element_comp) {
	// Create the linear element
	struct Element_Linear ele = create_element_L3(0, 1.3, 3);

	gsl_vector* v = output_constant_vector(&ele, driving_func);
	gsl_matrix* m = output_coefficient_matrix(&ele, a, b);

	// Loading the matrices and vectors
	gsl_vector* v_check = gsl_vector_alloc(3);
	gsl_matrix* m_check = gsl_matrix_alloc(3, 3);
	for (int i = 0; i < 3; i++) {
		gsl_vector_set(v_check, i, L2_F[i]);
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			gsl_matrix_set(m_check, i, j, L2_K[i][j]);
		}
	}

	ck_assert_mem_eq(v->data, v_check->data, 3);
	ck_assert_mem_eq(m->data, m_check->data, 9);

} 
END_TEST


// Setup and Execution
Suite* composition_suite() {
	Suite* s;
	TCase* tc_core;

	s = suite_create("Matrix and Vector Composition Tests");

	/* Core test case */
	tc_core = tcase_create("Core");

	tcase_add_test(tc_core, L2_element_comp);
	tcase_add_test(tc_core, L3_element_comp);
	suite_add_tcase(s, tc_core);

	return s;
}


int main() {
	// Assign directory here
	int number_failed;
	Suite *s_parser;
	SRunner *sr_parser;

	// Set the suites here.
	s_parser = composition_suite();

	// Setting the runners here
	sr_parser = srunner_create(s_parser);

	// Run the suites here.
	srunner_run_all(sr_parser, CK_NORMAL);

	number_failed = srunner_ntests_failed(sr_parser);

	srunner_free(sr_parser);

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
