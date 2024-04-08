#include "thresfilterpthreads.h"
#include <stdio.h>
#include <stdlib.h>
#include "ppmio.h"
#include <pthread.h>
#include "pixel.h"

#define uint unsigned int

struct thread_data
{
	int threadId;

	pixel *src;
	int size;
};

#define NUM_THREADS 16
pthread_t threads[NUM_THREADS];
struct thread_data thread_data_array[NUM_THREADS];

uint averageThreadResult[NUM_THREADS];
uint averageResult;

void *averageThread(void *tParam)
{
	struct thread_data *t;
	t = (struct thread_data *)tParam;
	int size = t->size / NUM_THREADS;
	int start = t->threadId * size;
	int end = (t->threadId - 1) == NUM_THREADS ? t->size : start + size;
	// If last thread, do extra work to make sure all is covered, if some was lost to rounding.

	uint sum = 0;
	printf("Thread %d, start %d to end %d\n", t->threadId, start, end);
	for (int i = start; i < end; i++)
	{
		pixel *px = (t->src + i);
		sum += (uint)px->r + (uint)px->g + (uint)px->b;
	}
	averageThreadResult[t->threadId] = sum / size;
}

void *thresholdThread(void *tParam)
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
		pixel *px = (t->src + i);
		uint v = (uint)px->r + (uint)px->g + (uint)px->b;
		if (v < averageResult)
		{
			px->r = px->g = px->b = 0;
		}
		else
		{
			px->r = px->g = px->b = 255;
		}
	}
}

void thresfilterPthreads(const int xsize, const int ysize, pixel *src)
{
	printf("Image %d x %d = %d\n", xsize, ysize, xsize * ysize);

	printf("=== BEGIN AVERAGE ===\n");
	for (size_t i = 0; i < NUM_THREADS; i++)
	{
		thread_data_array[i].threadId = i;
		thread_data_array[i].size = xsize * ysize;
		thread_data_array[i].src = src;

		pthread_create(&threads[i], NULL, averageThread, (void *)&thread_data_array[i]);
	}

	averageResult = 0;

	for (int i = 0; i < NUM_THREADS; i++)
	{
		pthread_join(threads[i], NULL);
		averageResult += averageThreadResult[i] / NUM_THREADS;
	}

	printf("=== BEGIN VERTICAL ===\n");
	for (size_t i = 0; i < NUM_THREADS; i++)
	{
		pthread_create(&threads[i], NULL, thresholdThread, (void *)&thread_data_array[i]);
	}

	for (int i = 0; i < NUM_THREADS; i++)
	{
		pthread_join(threads[i], NULL);
	}
	printf("=== END COMPUTATION ===\n");
}
