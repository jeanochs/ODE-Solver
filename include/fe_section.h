#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>


typedef enum {
	LINEAR,
	QUAD
} Element_2D_Type;

typedef union {
	// Used with the LINEAR tag
	struct {
		double node_coord[2];
		double (*shape_func[2]) (double); // Shape functions
		double (*shape_derv[2]) (double); // Shape function derivatives
		gsl_function jacobian;
		gsl_function iso_conversion; // Function that converts physical coordiantes to isoparametric relationships
	} L2;

	// Used with the QUAD (quadratic) tag
	struct {
		double node_coord[3];
		double (*shape_func[3]) (double); // Shape functions
		double (*shape_derv[3]) (double); // Shape function derivatives
		gsl_function jacobian;
		gsl_function iso_conversion; // Function that converts physical coordiantes to isoparametric relationships
	} L3;

} Element_2D;

struct Element_Linear {
	Element_2D_Type kind;
	Element_2D element;
};

typedef union {
	struct {int node_id[2]; } L2;
	struct {int node_id[3]; } L3;

} Node_Conn_2D;

struct Element_Conn {
	Element_2D_Type kind;
	Node_Conn_2D node_list;

};

struct Mesh {
	struct Element_Conn* connectivity_grid;
	struct Element_Linear* elements;
	uint16_t num_nodes;
	uint16_t num_elements;

};

struct ODE_Solution {
	gsl_vector* solution_coeff;
	// These fields can be NULL; it is an optional output
	gsl_matrix* coeff_matrix_global;
	gsl_vector* const_vector_global;

};

/* Function Prototypes */

// Main Functions
int parse_input_file(FILE* input_stream, struct Mesh* mesh_object, Element_2D_Type mesh_kind);
int solve_ode_constant(struct Mesh* input_mesh, struct ODE_Solution* solution, double a, double b, double d1, double d2, double (*func) (double), bool output_global_arrays);

// Creation Functions
gsl_vector* output_constant_vector(struct Element_Linear* element, double (*driving_func) (double));
gsl_matrix* output_coefficient_matrix(struct Element_Linear* element, double a, double b);
void create_element_L2(struct Element_Linear* e, double node1, double node2);
void create_element_L3(struct Element_Linear* e, double node1, double node2, double node3);
void free_element_memory(struct Element_Linear* ele);
void free_mesh_memory(struct Mesh* input_mesh);
void free_solution_memory(struct ODE_Solution* solution);

// Composition functions

struct Constant_Vector_Funcs {
	double (*nonhomo) (double);
	double (*shape) (double);
	gsl_function* converter;
	gsl_function* jacobian;

};

struct Coefficient_Matrix_Funcs {
	double a, b;
	double (*shape_i) (double);
	double (*shape_j) (double);
	double (*shape_derv_i) (double);
	double (*shape_derv_j) (double);
	gsl_function* jacobian;

};

double constant_vector_composition(double zeta, void* func_params);
double coefficient_matrix_composition(double zeta, void* func_params);


