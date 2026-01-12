// Header file for the Function Field struct and its associated routines
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

struct Function_Field {
	double *f_values;
	double *x_values;
	size_t number_of_points;
	double step_size;

};

int create_function_field(struct Function_Field *field, double start, double end, double number_of_points, double (*generating_func) (double));
int output_function_field(struct Function_Field *field, char *filename);
int f_eval(struct Function_Field *field, double x, double *f);
void free_function_field(struct Function_Field *field);

