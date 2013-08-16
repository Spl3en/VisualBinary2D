#pragma once

// Includes

// Defines

// Structure
typedef struct
{
	float *values;
	float increment;
	float start;
	float range;

} LookUpTable;


// Prototypes
LookUpTable *
lookup_table_new (double start, int range, double increment, double (*function)(double));
double lookup_table_get (LookUpTable *table, double value);
double lookup_table_get_safe (LookUpTable *table, double value);
