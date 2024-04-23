
//-----------------------------------------------------------------------
// Serial program for solving the heat conduction problem
// on a square using the Jacobi method.
// Written by August Ernstsson 2015-2019
//-----------------------------------------------------------------------

#define _POSIX_C_SOURCE 199309L
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>

double timediff(struct timespec *begin, struct timespec *end)
{
	double sec = 0.0, nsec = 0.0;
	if ((end->tv_nsec - begin->tv_nsec) < 0)
	{
		sec = (double)(end->tv_sec - begin->tv_sec - 1);
		nsec = (double)(end->tv_nsec - begin->tv_nsec + 1000000000);
	}
	else
	{
		sec = (double)(end->tv_sec - begin->tv_sec);
		nsec = (double)(end->tv_nsec - begin->tv_nsec);
	}
	return sec + nsec / 1E9;
}

void printm(int n, double *M)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%f\t", *(M + n * i + j));
		}
		printf("\n");
	}
	printf("\n");
}

void arrcpy(double *dst, double *src, int len)
{
	for (int it = 0; it < len; it++)
		dst[it] = src[it];
}

void laplsolv(int n, int maxiter, double tol)
{
	double T[n + 2][n + 2];
	double row_above[n], tmp2[n], row_after_end[n];
	int k;

	struct timespec starttime, endtime;

	// Set boundary conditions and initial values for the unknowns
	for (int i = 0; i <= n + 1; ++i)
	{
		for (int j = 0; j <= n + 1; ++j)
		{
			if (i == n + 1)
				T[i][j] = 2;
			else if (j == 0 || j == n + 1)
				T[i][j] = 1;
			else
				T[i][j] = 0;
		}
	}

	clock_gettime(CLOCK_MONOTONIC, &starttime);

	// Solve the linear system of equations using the Jacobi method
	for (k = 0; k < maxiter; ++k)
	{
		double error = -INFINITY;

		// Copy to temp buffers

#pragma omp parallel private(row_above, tmp2, row_after_end)
		{
			int me = omp_get_thread_num();
			int p = omp_get_num_threads();
			int rows_per = n / p;
			double my_error = -INFINITY;

			int start = me * rows_per;
			int end = start + rows_per;
			double below_value;
			if (me == p - 1)
			{
				end = n;
			}

			arrcpy(row_above, &T[start][1], n);
			arrcpy(row_after_end, &T[end + 1][1], n);

#pragma omp barrier

			for (int i = start + 1; i <= end; ++i)
			{

				// Apply the Jacobi algorithm to each element in this row
				for (int j = 1; j <= n; ++j)
				{
					if (i == end)
					{
						below_value = row_after_end[j - 1];
					}
					else
					{
						below_value = T[i + 1][j];
					}
					tmp2[j - 1] = (T[i][j - 1] + T[i][j + 1] + below_value + row_above[j - 1]) / 4.0;
					my_error = fmax(my_error, fabs(T[i][j] - tmp2[j - 1]));
				}
				arrcpy(row_above, &T[i][1], n);
				arrcpy(&T[i][1], tmp2, n);
			}
#pragma omp critical
			{
				error = fmax(error, my_error);
			}
		}
		if (error < tol)
			break;
	}

	clock_gettime(CLOCK_MONOTONIC, &endtime);

	printf("Time: %f\n", timediff(&starttime, &endtime));
	printf("Number of iterations: %d\n", k);
	printf("Temperature of element T(1,1): %.17f\n", T[1][1]);
}

int main(int argc, char *argv[])
{
	if (argc != 4)
	{
		printf("Usage: %s [size] [maxiter] [tolerance] \n", argv[0]);
		exit(1);
	}

	int size = atoi(argv[1]);
	int maxiter = atoi(argv[2]);
	double tol = atof(argv[3]);

	printf("Size %d, max iter %d and tolerance %f.\n", size, maxiter, tol);
	laplsolv(size, maxiter, tol);
	return 0;
}