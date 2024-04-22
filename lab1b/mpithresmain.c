#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "pixel.h"
#include "thresfiltermpi.h"
#include "mpi.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

#define PPM "im1.ppm"
#define IN "./data/"
#define OUT "./out/"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Datatype PIXEL_MPI;
	{
		pixel item;
		MPI_Datatype block_types[] = {MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR};
		int block_lenghts[] = {1, 1, 1};
		MPI_Aint start, displacement[3];

		MPI_Get_address(&item, &start);
		MPI_Get_address(&item.r, &displacement[0]);
		displacement[0] -= start;
		MPI_Get_address(&item.g, &displacement[1]);
		displacement[1] -= start;
		MPI_Get_address(&item.b, &displacement[2]);
		displacement[2] -= start;

		MPI_Type_create_struct(3, block_lenghts, displacement, block_types, &PIXEL_MPI);
		MPI_Type_commit(&PIXEL_MPI);
	}

	int p;	// number of started MPI processes
	int me; // my rank
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);

	int xsize, ysize, colmax;
	unsigned size;
	pixel *src;
	double stime, etime;
	char file[100];

	if (me == 0)
	{
		printf("Reading image\n");

		src = (pixel *)malloc(sizeof(pixel) * MAX_PIXELS);
		sprintf(file, "%s%s", IN, PPM);

		/* Read file */
		if (read_ppm(file, &xsize, &ysize, &colmax, (char *)src) != 0)
			exit(1);
		if (colmax > 255)
		{
			fprintf(stderr, "Too large maximum color-component value\n");
			exit(1);
		}
		printf("Has read the image with size %dx%d\n", xsize, ysize);
		size = xsize * ysize;
		stime = MPI_Wtime();
	}

	MPI_Bcast(&size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	// Use overestimate of chunksize to handle unevenly divisible
	unsigned max_chunk_size = size / p + 1;
	pixel *src_chunk = (pixel *)malloc(sizeof(pixel) * max_chunk_size);

	MPI_Scatter(src, max_chunk_size, PIXEL_MPI, src_chunk, max_chunk_size, PIXEL_MPI, 0, MPI_COMM_WORLD);
	if (me == 0)
	{
		free(src);
	}
	// chunks are oversized, make sure to not go out of bounds
	unsigned chunk_size = MIN(max_chunk_size * (me + 1), size) - max_chunk_size * me;

	unsigned chunk_sum = sumMPI(chunk_size, src_chunk);
	unsigned sum = 0;
	MPI_Allreduce(&chunk_sum, &sum, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

	unsigned average = sum / size;
	
	thresfilterMPI(chunk_size, average, src_chunk);

	pixel *dst;
	if (me == 0)
	{
		dst = (pixel *)malloc(sizeof(pixel) * max_chunk_size * p);
	}
	MPI_Gather(src_chunk, max_chunk_size, PIXEL_MPI, dst, max_chunk_size, PIXEL_MPI, 0, MPI_COMM_WORLD);
	free(src_chunk);

	if (me == 0)
	{
		etime = MPI_Wtime();
		printf("Filtering took: %g secs\n", (etime - stime));

		/* Write result */
		printf("Writing output file\n");

		sprintf(file, "%sthreshold-final-%s", OUT, PPM);
		if (write_ppm(file, xsize, ysize, (char *)dst) != 0)
			exit(1);
		free(dst);
	}
	MPI_Finalize();
	return 0;
}
