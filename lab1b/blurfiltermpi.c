/*
  File: blurfilter.c
  Implementation of blurfilter function.
 */

#include <stdio.h>
#include <stdlib.h>
#include "blurfiltermpi.h"
#include "pixel.h"
#include "mpi.h"

void gather(int index, pixel *src, pixel *dst, int radius, const double *w, int jump)
{
	double r = 0, g = 0, b = 0, n = 0, wc = 0;
	for (int wi = -radius; wi <= radius; wi++)
	{
		// cache optimized if jump is 1
		int offset = index + (wi * jump); // Jump is 1 if x or xsize if y

		wc = w[abs(wi)];
		pixel *pixel = (src + offset);
		r += wc * pixel->r;
		g += wc * pixel->g;
		b += wc * pixel->b;
		n += wc;
	}
	(dst + index)->r = r / n;
	(dst + index)->g = g / n;
	(dst + index)->b = b / n;
}

void blurfilterMPI(const int xsize, const int ysize, pixel *src, const int radius, const double *w)
{
	pixel *dst = (pixel *)malloc(sizeof(pixel) * xsize * ysize);

	// transfer borders
	for (int y = 0; y < ysize; y++)
	{
		for (int x = 0; x < xsize; x++)
		{
			if (x < radius || x > xsize - radius || y < radius || y > ysize - radius)
			{
				int i = x + y * xsize;
				dst[i] = src[i];
			}
		}
	}

	// Horizontal
	for (int y = radius; y < ysize - radius; y++)
	{
		for (int x = radius; x < xsize - radius; x++)
		{
			int i = x + y * xsize;
			gather(i, src, dst, radius, w, 1);
		}
	}

	// Vertical
	for (int y = radius; y < ysize - radius; y++)
	{
		for (int x = radius; x < xsize - radius; x++)
		{
			int i = x + y * xsize;
			gather(i, dst, src, radius, w, xsize);
		}
	}

	free(dst);
}