#include "LookUpTable.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

LookUpTable *lookup_table_new (double start, int range, double increment, double (*function)(double))
{
	LookUpTable *this = malloc(sizeof(LookUpTable));
	if (!this)
		return NULL;
	
	this->start     = start;
	this->range     = range;
	this->increment = increment;
	this->values    = malloc(sizeof(double) * range * (1.0 / increment));

	if (!this->values)
	{
		printf("%s : allocation of %d elements failed\n", __FUNCTION__, range);
		return NULL;
	}

	double i; int pos;
	for (i = start, pos = 0;
		 i < start + range;
		 i += increment, pos++)
	{
		this->values[pos] = function(i);
	}

	return this;
}

double lookup_table_get (LookUpTable *this, double value)
{
	return this->values[(int)(value * (1.0 / this->increment) - this->start) + 1];
}

double lookup_table_get_safe (LookUpTable *this, double value)
{
	if (value >= this->start
	&&  value <= this->start + this->range)
	{
		return this->values[(int)(value * (1.0 / this->increment) - this->start) + 1];
	}
	else
	{
		printf("%s: Out of lookuptable range (%f out of <%f-%f>)\n", __FUNCTION__, value, this->start, this->start + this->range);
		return 0.0;
	}
}

