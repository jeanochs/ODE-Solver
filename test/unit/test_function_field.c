#include <stdlib.h>
#include <check.h>
#include <string.h>

#include "fe_section.h"

#define TOL 1e-3

double test_function(double x) {
	return 3*x*x - 5*x + 3;

}

START_TEST(Function_Field_Check) {
	// This would normally be handled by a dedicated function, but this hasn't been written yet.
	// This setup will probably stay, though
	// Create the (temporary) Function Field object
	struct Function_Field f_field_temp;

	// Create the driving function points
	create_function_field(&f_field_temp, 0, 10, 1001, test_function);

	// Calculate some points
	double ans1, ans2, ans3;
	f_eval(&f_field_temp, 2, &ans1);
	f_eval(&f_field_temp, 5.5, &ans2);
	f_eval(&f_field_temp, 4.67, &ans3);

	// Compare
	ck_assert_double_eq_tol(ans1, 5, TOL);
	ck_assert_double_eq_tol(ans2, 66.25, TOL);
	ck_assert_double_eq_tol(ans3, 45.0767, TOL);

	// Free up the thing
	free_function_field(&f_field_temp);

}
END_TEST

// Setup and Execution
Suite* function_suite() {
	Suite* s;
	TCase* tc_core;

	s = suite_create("Function Field Test");

	/* Core test case */
	tc_core = tcase_create("Core");

	tcase_add_test(tc_core, Function_Field_Check);
	suite_add_tcase(s, tc_core);

	return s;
}



int main() {
	// Assign directory here
	int number_failed;
	Suite *s_function;
	SRunner *sr_function;

	// Set the suites here.
	s_function = function_suite();

	// Setting the runners here
	sr_function = srunner_create(s_function);

	// Run the suites here.
	srunner_run_all(sr_function, CK_NOFORK);

	number_failed = srunner_ntests_failed(sr_function);

	srunner_free(sr_function);

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
