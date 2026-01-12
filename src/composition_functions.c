/* Ancillary source file for the composition functions that will be used in integration
 *
 */

/* ODE Solver */
/* Internal Function Parameter Data Structures
 * These are forward defined since they are temporary objects.
 * There structure will have to be inferred from the below composition functions that use them.
 */

double constant_vector_composition(double zeta, void* func_params) {
	struct Constant_Vector_Funcs* p = (struct Constant_Vector_Funcs*) func_params;
	struct Function_Field *data_field = p->function_field;
	double (*shape) (double) = p->shape;

	double jacobian = GSL_FN_EVAL(p->jacobian, zeta);
	double x_value = GSL_FN_EVAL(p->converter, zeta); // Convert from ispoarametric to physical
	double f_result;
	f_eval(data_field, x_value, &f_result);

	return f_result*shape(zeta)*jacobian;

}

double coefficient_matrix_composition(double zeta, void* func_params) {
	struct Coefficient_Matrix_Funcs* p = (struct Coefficient_Matrix_Funcs*) func_params;

	// The first -1 below has caused me a lot of problems forgetting it. 
	double term1 = -1*p->shape_derv_i(zeta)*p->shape_derv_j(zeta)*(1/GSL_FN_EVAL(p->jacobian, zeta));
	double term2 = p->a*p->shape_i(zeta)*p->shape_derv_j(zeta);
	double term3 = p->b*p->shape_i(zeta)*p->shape_j(zeta)*GSL_FN_EVAL(p->jacobian, zeta);

	return term1 + term2 + term3;

}





