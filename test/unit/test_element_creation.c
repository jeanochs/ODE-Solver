#include <stdlib.h>
#include <check.h>
#include <string.h>

#include "fe_section.h"

#define TOL 1e-3

START_TEST(L2_element_creation) {
	// Create the linear element
	struct Element_Linear ele;
	create_element_L2(&ele, 1, 3);

	// Calculate the shape functions and their derivatives
	double N0 = ele.element.L2.shape_func[0](-0.3);
	double N1 = ele.element.L2.shape_func[1](0.43);
	double N0_D = ele.element.L2.shape_derv[0](-0.1);
	double N1_D = ele.element.L2.shape_derv[1](-0.5);

	// Calculate the Jacobian
	double J = GSL_FN_EVAL(&ele.element.L2.jacobian, 0);

	// Calculate the physical coordiante from the isoparamtric one
	double phy_x = GSL_FN_EVAL(&ele.element.L2.iso_conversion, 0.1);

	// Check the results
	ck_assert_double_eq_tol(N0, 0.65, TOL);
	ck_assert_double_eq_tol(N1, 0.715, TOL);
	ck_assert_double_eq_tol(N0_D, -0.5, TOL);
	ck_assert_double_eq_tol(N1_D, 0.5, TOL);
	ck_assert_double_eq_tol(J, 1.0f, TOL);
	ck_assert_double_eq_tol(phy_x, 2.1, TOL);

	free_element_memory(&ele);

} 
END_TEST

START_TEST(L3_element_creation) {
	// Create the linear element
	struct Element_Linear ele;
	create_element_L3(&ele, 0, 1.3, 3);

	// Calculate the shape functions and their derivatives
	double N0 = ele.element.L3.shape_func[0](-0.1);
	double N1 = ele.element.L3.shape_func[1](-0.2);
	double N2 = ele.element.L3.shape_func[2](0.33);
	double N0_D = ele.element.L3.shape_derv[0](0.13);
	double N1_D = ele.element.L3.shape_derv[1](0.5);
	double N2_D = ele.element.L3.shape_derv[2](0.6);

	// Calculate the Jacobian
	double J = GSL_FN_EVAL(&ele.element.L3.jacobian, 0.15);

	// Calculate the physical coordiante from the isoparamtric one
	double phy_x = GSL_FN_EVAL(&ele.element.L3.iso_conversion, 0.2);

	ck_assert_double_eq_tol(N0, 0.055, TOL);
	ck_assert_double_eq_tol(N1, 0.96, TOL);
	ck_assert_double_eq_tol(N2, 0.2194, TOL);
	ck_assert_double_eq_tol(N0_D, -0.37, TOL);
	ck_assert_double_eq_tol(N1_D, -1, TOL);
	ck_assert_double_eq_tol(N2_D, 1.1, TOL);
	ck_assert_double_eq_tol(J, 1.56, TOL);
	ck_assert_double_eq_tol(phy_x, 1.608, TOL);

	free_element_memory(&ele);

} 
END_TEST

// Setup and Execution
Suite* element_creation_suite() {
	Suite* s;
	TCase* tc_core;

	s = suite_create("Element Creation Tests");

	/* Core test case */
	tc_core = tcase_create("Core");

	tcase_add_test(tc_core, L2_element_creation);
	tcase_add_test(tc_core, L3_element_creation);
	suite_add_tcase(s, tc_core);

	return s;
}


int main() {
	// Assign directory here
	int number_failed;
	Suite *s_parser;
	SRunner *sr_parser;

	// Set the suites here.
	s_parser = element_creation_suite();

	// Setting the runners here
	sr_parser = srunner_create(s_parser);

	// Run the suites here.
	srunner_run_all(sr_parser, CK_NORMAL);

	number_failed = srunner_ntests_failed(sr_parser);

	srunner_free(sr_parser);

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
