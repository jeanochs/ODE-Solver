#include <stdlib.h>
#include <check.h>
#include <string.h>

#include "fe_section.h"

#define TOL 1e-3

double driving_func(double x) {
	return (x*x) + x + 3;

}

// Constants
double a = 4.0f;
double b = 4.0f;

// Predefined matrices
// L2 Element Matrices and Vectors
double L2_K[2][2] = {
	{7./6, 17./6},
	{-7./6, 31./6}
};

double L2_F[2] = {23./3, 11.};

// L3 Element Matrices and Vectors
double L3_K[3][3] = {
	{0.2099, 2.2006, -0.9438},
	{-3.133, 8.2577, 2.8750},
	{0.3895, -2.4583, 4.6021}
};

double L3_F[3] = {0.4423, 13.7783, 8.2794};

// Answer fields for composition functions
// Inputs
double L2_constant_comp_input[2] = {-0.35, 0.2};
double L2_coefficient_comp_input[2][2] = {
	{-0.75, 0.1},
	{0.75, -0.2}
};
// Answers
double L2_constant_comp[2] = {4.97644, 6.024};
double L2_coefficient_comp[2][2] = {
	{25./16, 41./25},
	{-25./16, 169./100}
};

// Inputs
double L3_constant_comp_input[3] = {-0.3, 0.35, 0.55};
double L3_coefficient_comp_input[3][3] = {
	{-0.1, -0.5, -0.75},
	{0.3, 0.65, -0.1},
	{0, 0.9, -0.53}
};
// Answers
double L3_constant_comp[3] = {1.2436, 11.9016, 7.3036};
double L3_coefficient_comp[3][3] = {
	{0.1322, 2.1933, -0.6911},
	{-1.273, 0.3051, 1.3786},
	{-0.1667, -6.3022, 0.0956}
};

// Tests to ensure proper construction of matrix and vector integrands
START_TEST(L2_element_comp) {
	// Create element
	struct Element_Linear ele;
	create_element_L2(&ele, 1, 3);

	// Check the composition function for the constant vector
	for (int i = 0; i < 2; i++) {
		struct Constant_Vector_Funcs p = {
			driving_func,
			ele.element.L2.shape_func[i],
			&ele.element.L2.iso_conversion,
			&ele.element.L2.jacobian
		};

		double comp_value = constant_vector_composition(
				L2_constant_comp_input[i],
				&p
		);

		// Compare the values
		printf("Constant Vector Composition Function, L2 %d\n", i);
		ck_assert_double_eq_tol(comp_value, L2_constant_comp[i], TOL);

	}

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			struct Coefficient_Matrix_Funcs p = {
				a, b, 
				ele.element.L2.shape_func[i],
				ele.element.L2.shape_func[j],
				ele.element.L2.shape_derv[i],
				ele.element.L2.shape_derv[j],
				&ele.element.L2.jacobian
			};

			double comp_value = coefficient_matrix_composition(L2_coefficient_comp_input[i][j], &p);

			// Compare the values
			printf("Coefficient Matrix Composition Function L2 %d,%d\n", i, j);
			ck_assert_double_eq_tol(comp_value, L2_coefficient_comp[i][j], TOL);

		}
	}

	free_element_memory(&ele);

}
END_TEST

START_TEST(L3_element_comp) {
	// Create element
	struct Element_Linear ele;
	create_element_L3(&ele, 0, 1.3, 3);

	// Check the composition function for the constant vector
	for (int i = 0; i < 3; i++) {
		struct Constant_Vector_Funcs p = {
			driving_func,
			ele.element.L3.shape_func[i],
			&ele.element.L3.iso_conversion,
			&ele.element.L3.jacobian

		};

		double comp_value = constant_vector_composition(
				L3_constant_comp_input[i],
				&p
		);

		// Compare the values
		printf("Constant Vector Composition Function, L3 %d\n", i);
		ck_assert_double_eq_tol(comp_value, L3_constant_comp[i], TOL);

	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			struct Coefficient_Matrix_Funcs p = {
				a, b, 
				ele.element.L3.shape_func[i],
				ele.element.L3.shape_func[j],
				ele.element.L3.shape_derv[i],
				ele.element.L3.shape_derv[j],
				&ele.element.L3.jacobian
			};

			double comp_value = coefficient_matrix_composition(L3_coefficient_comp_input[i][j], &p);

			// Compare the values
			printf("Coefficient Matrix Composition Function L3 %d,%d\n", i, j);
			ck_assert_double_eq_tol(comp_value, L3_coefficient_comp[i][j], TOL);

		}
	}

	free_element_memory(&ele);

}
END_TEST

// Tests for element-wise matrix and vector assembly
START_TEST(L2_element_assembly) {
	// Create element
	struct Element_Linear ele;
	create_element_L2(&ele, 1, 3);

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

	// Now check everything by iteration
	printf("Checking data for the vector assembly of L2\n");
	for (int i = 0; i < 2; i++) {
		ck_assert_double_eq_tol(
			gsl_vector_get(v, i),
			gsl_vector_get(v_check, i),
			TOL);
	}

	printf("Checking data for the matrix assembly of L2\n");
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			ck_assert_double_eq_tol(
				gsl_matrix_get(m, i, j),
				gsl_matrix_get(m_check, i, j),
				TOL);
		}
	}

	free_element_memory(&ele);
	gsl_vector_free(v); gsl_vector_free(v_check);
	gsl_matrix_free(m); gsl_matrix_free(m_check);

} 
END_TEST

START_TEST(L3_element_assembly) {
	// Create the linear element
	struct Element_Linear ele;
	create_element_L3(&ele, 0, 1.3, 3);

	gsl_vector* v = output_constant_vector(&ele, driving_func);
	gsl_matrix* m = output_coefficient_matrix(&ele, a, b);

	// Loading the matrices and vectors
	gsl_vector* v_check = gsl_vector_alloc(3);
	gsl_matrix* m_check = gsl_matrix_alloc(3, 3);
	for (int i = 0; i < 3; i++) {
		gsl_vector_set(v_check, i, L3_F[i]);
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			gsl_matrix_set(m_check, i, j, L3_K[i][j]);
		}
	}

	// Now check everything by iteration
	printf("Checking data for the vector assembly of L3\n");
	for (int i = 0; i < 3; i++) {
		ck_assert_double_eq_tol(
			gsl_vector_get(v, i),
			gsl_vector_get(v_check, i),
			TOL);
	}

	printf("Checking data for the matrix assembly of L3\n");
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			ck_assert_double_eq_tol(
				gsl_matrix_get(m, i, j),
				gsl_matrix_get(m_check, i, j),
				TOL);
		}
	}

	free_element_memory(&ele);
	gsl_vector_free(v); gsl_vector_free(v_check);
	gsl_matrix_free(m); gsl_matrix_free(m_check);

} 
END_TEST


// Setup and Execution
Suite* composition_suite() {
	Suite* s;
	TCase* tc_core;

	s = suite_create("Integrand Composition Tests");

	/* Core test case */
	tc_core = tcase_create("Core");

	tcase_add_test(tc_core, L2_element_comp);
	tcase_add_test(tc_core, L3_element_comp);
	suite_add_tcase(s, tc_core);

	return s;
}

Suite* local_assembly_suite() {
	Suite* s;
	TCase* tc_core;

	s = suite_create("Matrix and Vector Composition Tests");

	/* Core test case */
	tc_core = tcase_create("Core");

	tcase_add_test(tc_core, L2_element_assembly);
	tcase_add_test(tc_core, L3_element_assembly);
	suite_add_tcase(s, tc_core);

	return s;
}


int main() {
	// Assign directory here
	int number_failed;
	Suite *s_composition, *s_assembly;
	SRunner *sr_composition, *sr_assembly;

	// Set the suites here.
	s_composition = composition_suite();
	s_assembly = local_assembly_suite();

	// Setting the runners here
	sr_composition = srunner_create(s_composition);
	sr_assembly = srunner_create(s_assembly);

	// Run the suites here.
	srunner_run_all(sr_composition, CK_NORMAL);
	srunner_run_all(sr_assembly, CK_NORMAL);

	number_failed = srunner_ntests_failed(sr_composition);
	number_failed += srunner_ntests_failed(sr_assembly);

	srunner_free(sr_composition);
	srunner_free(sr_assembly);

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
