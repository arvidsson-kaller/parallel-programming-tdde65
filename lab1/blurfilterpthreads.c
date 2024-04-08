/*
  File: blurfilter.c
  Implementation of blurfilter function.
 */

#include <stdio.h>
#include <stdlib.h>
#include "blurfilterpthreads.h"
#include "ppmio.h"
#include <pthread.h>
#include "pixel.h"

struct thread_data
{
	int threadId;

	pixel *src;
	pixel *dst;
	int size;
	int radius;
	double const *w;
	int jump;
};

#define NUM_THREADS 16
pthread_t threads[NUM_THREADS];
struct thread_data thread_data_array[NUM_THREADS];

void gather(int index, pixel *src, pixel *dst, int size, int radius, const double *w, int jump)
{
	double r = 0, g = 0, b = 0, n = 0, wc = 0;
	// printf("Pixel %d: \n", index);
	for (int wi = -radius; wi <= radius; wi++)
	{
		// cache optimized if jump is 1
		int offset = index + (wi * jump); // Jump is 1 if x or xsize if y
		if (offset >= 0 && offset < size)
		{
			wc = w[wi];
			pixel *pixel = (src + offset);
			// printf("(%d, %d, %d), ", index, pixel->r, pixel->g, pixel->b);
			r += wc * pixel->r;
			g += wc * pixel->g;
			b += wc * pixel->b;
			n += wc;
		}
	}
	// printf("\nfinished: (%f, %f, %f) with n: %f\n", r,g,b, n);
	(dst + index)->r = r / n;
	(dst + index)->g = g / n;
	(dst + index)->b = b / n;
}

void *blurFilterThread(void *tParam)
{
	struct thread_data *t;
	t = (struct thread_data *)tParam;
	int size = t->size / NUM_THREADS;
	int start = t->threadId * size;
	int end = (t->threadId - 1) == NUM_THREADS ? t->size : start + size;
	// If last thread, do extra work to make sure all is covered, if some was lost to rounding.

	printf("Thread %d, start %d to end %d\n", t->threadId, start, end);
	for (int i = start; i < end; i++)
	{
		gather(i, t->src, t->dst, t->size, t->radius, t->w, t->jump);
	}
}

void blurfilterPthreads(const int xsize, const int ysize, pixel *src, const int radius, const double *w)
{
	printf("Image %d x %d = %d\n", xsize, ysize, xsize * ysize);
	pixel *dst = (pixel *)malloc(sizeof(pixel) * MAX_PIXELS);

	printf("=== BEGIN HORIZONTAL ===\n");
	for (size_t i = 0; i < NUM_THREADS; i++)
	{
		thread_data_array[i].threadId = i;
		thread_data_array[i].size = xsize * ysize;
		thread_data_array[i].radius = radius;
		thread_data_array[i].w = w;

		thread_data_array[i].jump = 1;
		thread_data_array[i].src = src;
		thread_data_array[i].dst = dst;

		pthread_create(&threads[i], NULL, blurFilterThread, (void *)&thread_data_array[i]);
	}

	for (int i = 0; i < NUM_THREADS; i++)
	{
		pthread_join(threads[i], NULL);
	}

	printf("=== BEGIN VERTICAL ===\n");
	for (size_t i = 0; i < NUM_THREADS; i++)
	{
		thread_data_array[i].src = dst;
		thread_data_array[i].dst = src;
		thread_data_array[i].jump = xsize;

		pthread_create(&threads[i], NULL, blurFilterThread, (void *)&thread_data_array[i]);
	}

	for (int i = 0; i < NUM_THREADS; i++)
	{
		pthread_join(threads[i], NULL);
	}
	printf("=== END COMPUTATION ===\n");

	free(dst);
}
