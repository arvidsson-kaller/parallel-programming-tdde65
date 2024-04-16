#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "blurfiltermpi.h"
#include "gaussw.h"
#include "pixel.h"
#include "mpi.h"
#include <math.h>

#define RADIUS 100
#define PPM "im1.ppm"
#define IN "./data/"
#define OUT "./out/"

#define TAG_SIZE_MPI 1
#define TAG_DATA_MPI 2

int is_power_of_two(unsigned int x)
{
	// A power of two has only one bit set. So, if x is a power of two,
	// only one bit will be set, and x - 1 will have all the bits set
	// to the right of the original set bit. For example:
	// 8 (1000) -> 7 (0111)
	// 16 (10000) -> 15 (01111)
	// If x is 0, it will return false.
	return (x && !(x & (x - 1)));
}

void calculate_dimensions(int p, unsigned *cols, unsigned *rows)
{
	if (is_power_of_two(p) && p != 1)
	{
		unsigned diff = -1;
		for (size_t i = 1; i < p; i *= 2)
		{
			unsigned tmp = abs(i - p / i);
			if (tmp < diff)
			{
				diff = tmp;
				*cols = i;
				*rows = p / i;
			}
		}
	}
	else
	{
		double root = sqrt(p);
		unsigned flo = floor(root);
		unsigned cei = ceil(root);
		if (cei * flo == p)
		{
			*cols = cei;
			*rows = flo;
		}
		else
		{
			*cols = p;
			*rows = 1;
		}
	}
	if (*rows > *cols)
	{
		unsigned tmp = *cols;
		*cols = *rows;
		*rows = tmp;
	}
}

int divide_image(pixel *src, pixel *cells_source, unsigned xsize, unsigned ysize, unsigned cols, unsigned rows, unsigned max_cell_size, unsigned *cell_sizes)
{
	unsigned xsize_cell = xsize / cols;
	unsigned ysize_cell = ysize / rows;
	for (size_t row = 0; row < rows; row++)
	{
		for (size_t col = 0; col < cols; col++)
		{
			unsigned rank = (col + row * cols);
			pixel *cell_src = cells_source + (max_cell_size)*rank;

			int y_start = row * ysize_cell - RADIUS;
			int y_end = (row + 1) * ysize_cell + RADIUS;
			if (row == rows - 1)
			{
				// Handle extra pixels in the last col
				y_end = ysize + RADIUS;
			}
			int x_start = col * xsize_cell - RADIUS;
			int x_end = (col + 1) * xsize_cell + RADIUS;
			if (col == cols - 1)
			{
				// Handle extra pixels in the last row
				x_end = xsize + RADIUS;
			}

			for (int y = y_start; y < y_end; y++)
			{

				for (int x = x_start; x < x_end; x++)
				{
					pixel ptmp;
					ptmp.r = ptmp.g = ptmp.b = 255;

					if (x >= 0 && y >= 0 && x < xsize && y < ysize)
					{
						ptmp = src[x + y * xsize];
					}

					cell_src[(x - x_start) + (y - y_start) * (x_end - x_start)] = ptmp;
				}
			}
			cell_sizes[rank * 2] = x_end - x_start;
			cell_sizes[rank * 2 + 1] = y_end - y_start;
		}
	}
}

int main(int argc, char **argv)
{
	int xsize, ysize, colmax;

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

	double w[RADIUS+1];

	unsigned cols = 0;
	unsigned rows = 0;

	calculate_dimensions(p, &cols, &rows);

	char me_have_work = cols * rows > me;

	unsigned max_cell_size = 0;
	unsigned cell_sizes[cols * rows * 2];
	char file[100];
	struct timespec stime, etime;

	pixel *cells_source;
	if (me == 0)
	{
		printf("Splitting image into %dx%d grid\n", cols, rows);

		pixel *src = (pixel *)malloc(sizeof(pixel) * MAX_PIXELS);
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

		max_cell_size = (xsize / cols + cols + RADIUS * 2) * (ysize / rows + rows + RADIUS * 2);
		cells_source = (pixel *)malloc(sizeof(pixel) * (max_cell_size)*cols * rows);

		divide_image(src, cells_source, xsize, ysize, cols, rows, max_cell_size, cell_sizes);
		free(src);
	}

	unsigned cell_size[2];

	if (me_have_work)
	{
		get_gauss_weights(RADIUS, w);

		clock_gettime(CLOCK_REALTIME, &stime);
		MPI_Bcast(&max_cell_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		MPI_Scatter(cell_sizes, 2, MPI_UNSIGNED, &cell_size, 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		pixel *src = (pixel *)malloc(sizeof(pixel) * max_cell_size);

		MPI_Scatter(cells_source, max_cell_size, PIXEL_MPI, src, max_cell_size, PIXEL_MPI, 0, MPI_COMM_WORLD);
		if (me == 0)
		{
			free(cells_source);
		}

		unsigned xsize_cell = cell_size[0];
		unsigned ysize_cell = cell_size[1];

		blurfilterMPI(xsize_cell, ysize_cell, src, RADIUS, w);
		
		pixel *dst = (pixel *)malloc(sizeof(pixel) * max_cell_size);
		unsigned dst_i = 0;
		for (int y = RADIUS; y < ysize_cell - RADIUS; y++)
		{
			for (int x = RADIUS; x < xsize_cell - RADIUS; x++)
			{
				unsigned i = x + y * xsize_cell;
				dst[dst_i++] = src[i];
			}
		}
		free(src);


		if (me == 0)
		{
			src = (pixel *)malloc(sizeof(pixel) * max_cell_size * cols * rows);
		}
		MPI_Gather(dst, max_cell_size, PIXEL_MPI, src, max_cell_size, PIXEL_MPI, 0, MPI_COMM_WORLD);
		free(dst);

		if (me == 0)
		{
			dst = (pixel *)malloc(sizeof(pixel) * xsize * ysize);
			for (size_t rank = 0; rank < cols * rows; rank++)
			{
				pixel *rank_src = src + max_cell_size * rank;
				unsigned xsize_cell = cell_sizes[rank * 2] - RADIUS * 2;
				unsigned ysize_cell = cell_sizes[rank * 2 + 1] - RADIUS * 2;

				unsigned rank_col = rank % cols;
				unsigned rank_row = rank / cols;

				for (size_t cell_y = 0; cell_y < ysize_cell; cell_y++)
				{
					for (size_t cell_x = 0; cell_x < xsize_cell; cell_x++)
					{
						unsigned x = (xsize / cols) * rank_col + cell_x;
						unsigned y = (ysize / rows) * rank_row + cell_y;
						// Transfer into global from rank
						dst[x + y * xsize] = rank_src[cell_x + cell_y * xsize_cell];
					}
				}
			}

			clock_gettime(CLOCK_REALTIME, &etime);
			printf("Filtering took: %g secs\n", (etime.tv_sec - stime.tv_sec) +
													1e-9 * (etime.tv_nsec - stime.tv_nsec));

			sprintf(file, "%s%d-final-%s", OUT, me, PPM);
			if (write_ppm(file, xsize, ysize, (char *)dst) != 0)
				exit(1);
			free(dst);
		}
	}

	MPI_Finalize();
	return 0;
}
