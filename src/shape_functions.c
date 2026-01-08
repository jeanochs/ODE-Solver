// Source file purely for shape functions

/* 1D Elements
 *
 */

/* L2: Linear 1D Element (2 Nodes) */

/* Shape Functions */
double L2_N0(double zeta) {
	return (1 - zeta)/2;
}

double L2_N1(double zeta) {
	return (1 + zeta)/2;
}

/* Derivative of Shape Functions */
double L2_N0_D(double zeta) {
	return -1.0f/2;
}

double L2_N1_D(double zeta) {
	return 1.0f/2;
}

/* Jacobian Structure */
struct L2_N_P {
	double x1, x2;
};

double J_L2(double zeta, void* node_params) {
	struct L2_N_P* p = (struct L2_N_P*) node_params;
	double x1 = p->x1;
	double x2 = p->x2;

	return (x2 - x1)/2;

}

/* L3: Quadratic 1D Element (3 Nodes) */

/* Shape Functions */
double L3_N0(double zeta) {
	return zeta*(zeta - 1)/2;
}

double L3_N1(double zeta) {
	return 1 - zeta*zeta;
}

double L3_N2(double zeta) {
	return zeta*(zeta + 1)/2;
}

/* Derivative of Shape Functions */
double L3_N0_D(double zeta) {
	return zeta - 0.5;
}

double L3_N1_D(double zeta) {
	return -2*zeta;
}

double L3_N2_D(double zeta) {
	return zeta + 0.5;
}

/* Jacobian Structure */
struct L3_N_P {
	double x1, x2, x3;
};

double J_L3(double zeta, void* node_params) {
	struct L3_N_P* p = (struct L3_N_P*) node_params;
	double x1 = p->x1;
	double x2 = p->x2;
	double x3 = p->x3;

	return (zeta - 0.5)*x1 - 2*zeta*x2 + (zeta + 0.5)*x3;

}

// Helper Function
// Data structure used to assemble the phsucal-to-isoparametric functions.
// Since the conversion is just the addition of the products of the individual shape functions and nodes, the structure has arrays for iteration.
struct Iso_Phy_Funcs {
	uint8_t num_points; // This is a small type since I don't expect an element to have more than 20 points.
	double* x_coords;
	double (**shape_funcs) (double);
};

double isoparametric_physical_conversion(double zeta, void* params) {
	struct Iso_Phy_Funcs* p = (struct Iso_Phy_Funcs*) params;

	double sum = 0;
	for (int i = 0; i < p->num_points; i++) {
		sum += (p->x_coords[i]*p->shape_funcs[i](zeta));
	}

	return sum;

}
