#include "function_field.h"

int create_function_field(struct Function_Field *field, double start, double end, double number_of_points, double (*generating_func) (double)) {
	if (end <= start) {
		printf("The end value %f is less than or equal to the start value %f.\n", end, start);
		return 1;

	}

	double step_size = (end - start)/number_of_points;

	// Allocate the arrays
	field->x_values = malloc(number_of_points*sizeof(double));
	field->f_values = malloc(number_of_points*sizeof(double));

	for (int i = 0; i < number_of_points; i++) {
		double x_value = i*step_size;
		double f_value = generating_func(x_value);

		field->x_values[i] = x_value;
		field->f_values[i] = f_value;

	}

	field->number_of_points = number_of_points;
	field->step_size = step_size;

	return 0;

}

int input_function_field(struct Function_Field *field, FILE *file_stream) {
	// Parse the first line of the file to get the step size and number of points
	char buffer[300];

	if (fgets(buffer, 300, file_stream) == NULL) {
		fprintf(stderr, "There was an issue with opening the file. Please try again.\n");
		return 1;
	}

	// Split up
	size_t num_points = atoi(strtok(buffer, "\t"));
	double step_size = atof(strtok(NULL, "\t"));

	// Allocate the arrays
	double *x_point = malloc(num_points*sizeof(double));
	double *f_point = malloc(num_points*sizeof(double));

	// Parse through each line, get the numbers, and add to the arrays
	int counter = 0;
	while (fgets(buffer, 300, file_stream) != NULL) {
		// Split up the numbers
		double x_p = atof(strtok(buffer, "\t"));
		double f_p = atof(strtok(NULL, "\t"));

		x_point[counter] = x_p;
		f_point[counter] = f_p;
		counter++;

	}

	// Setup the Function Field struct
	field->f_values = f_point;
	field->x_values = x_point;
	field->number_of_points = num_points;
	field->step_size = step_size;

	return 0;

}


int output_function_field(struct Function_Field *field, char *filename) {
	FILE *dat_file = fopen(filename, "w");

	if (dat_file == NULL) {
		printf("Could not open the provided file %s' please check.\n", filename);
		return 1;

	}

	// Print the header
	double np = field->number_of_points;
	double ss = field->step_size;

	fprintf(dat_file, "%f\t%f\n", np, ss);

	for (int i = 0; i < field->number_of_points; i++) {
		double x = field->x_values[i];
		double f = field->f_values[i];
		fprintf(dat_file, "%f\t%f\n", x, f);

	}

	// Close the file once done
	fclose(dat_file);

	return 0;

}



int f_eval(struct Function_Field *field, double x, double *f) {
	// Check that the main value is within the bounds of the x-values
	if (x < field->x_values[0] || x > field->x_values[field->number_of_points - 1]) {
		printf("%f is not within the range given by the field.", x);
		return 1;
	}

	double starting_value = field->x_values[0];

	// First, divide by the step size and check if the result is an integer
	double indexing_number = (x - starting_value)/field->step_size;
	if (floor(indexing_number) == indexing_number) {
		*f = field->f_values[(int)indexing_number];

		return 0;

	}
	else {
		int lower_index = floor(indexing_number);
		int upper_index = floor(indexing_number) + 1;
		double slope = (field->f_values[upper_index] - field->f_values[lower_index])/(field->x_values[upper_index] - field->x_values[lower_index]);
		double interpolated_f_value = field->f_values[lower_index] + (x - field->x_values[lower_index])*slope;

		*f = interpolated_f_value;

		return 0;

	}

}

void free_function_field(struct Function_Field *field) {
	free(field->f_values);
	free(field->x_values);

}



	
