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
	double (*f) (double) = p->nonhomo;
	double (*shape) (double) = p->shape;

	// Since the function is given in terms of x, teh conversion from zeta to the crresponding x-value is done by the below equation.
	double x_value = GSL_FN_EVAL(p->converter, zeta);

	return f(x_value)*shape(zeta);

}

double coefficient_matrix_composition(double zeta, void* func_params) {
	struct Coefficient_Matrix_Funcs* p = (struct Coefficient_Matrix_Funcs*) func_params;

	double term1 = p->shape_derv_i(zeta)*p->shape_derv_j(zeta)*(1/GSL_FN_EVAL(p->jacobian, zeta));
	double term2 = p->a*p->shape_i(zeta)*p->shape_derv_j(zeta);
	double term3 = p->b*p->shape_i(zeta)*p->shape_j(zeta)*GSL_FN_EVAL(p->jacobian, zeta);

	return term1 + term2 + term3;

}





