#include "thresfiltermpi.h"
#include <stdio.h>
#include <stdlib.h>
#include "pixel.h"

#define uint unsigned int

uint averageMPI(const unsigned size, pixel *src)
{
	unsigned sum = 0;
	for (unsigned i = 0; i < size; i++)
	{
		pixel *px = (src + i);
		sum += (uint)px->r + (uint)px->g + (uint)px->b;
		// printf("Sum at iteration %i is %i\n", i, sum);
	}
	return sum / size;
}

void thresfilterMPI(const unsigned size, const unsigned average, pixel *src)
{
	for (int i = 0; i < size; i++)
	{
		pixel *px = (src + i);
		uint v = (uint)px->r + (uint)px->g + (uint)px->b;
		if (v < average)
		{
			px->r = px->g = px->b = 0;
		}
		else
		{
			px->r = px->g = px->b = 255;
		}
	}
}
