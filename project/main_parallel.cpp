// #include <stdlib.h>
#include <time.h>
#include <stdio.h>
// #include <math.h>
#include <stdbool.h>
#include <iostream>

#include "coordinate.h"
#include "definitions.h"
#include "physics.h"
#include "mpi.h"
#include "divide.h"
#include <vector>

// Define a macro for enabling/disabling debug prints
class DebugStream
{
public:
	template <typename T>
	DebugStream &operator<<(const T &value)
	{
		int my_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
		if (my_rank == 0)
			std::cout << value;
		return *this;
	}
};

struct particle_info
{
	pcord_t coords;
	bool collision;
};

DebugStream cout0;

int dims[2];
MPI_Comm GRID_COMM_MPI;
MPI_Datatype PCORD_MPI;

// Feel free to change this program to facilitate parallelization.

int grid_rank(int x, int y)
{
	int coords[2] = {x, y};
	int rank = -1;

	if (x >= 0 && x < dims[0] && y >= 0 && y < dims[1])
	{
		MPI_Cart_rank(GRID_COMM_MPI, coords, &rank);
	}
	return rank;
}

float rand1()
{
	return (float)(rand() / (float)RAND_MAX);
}

void wait_and_push(int rank, MPI_Request &req, pcord_t *recv_buf, std::vector<particle_info> &particles)
{
	if (rank != -1)
	{
		MPI_Status status;
		MPI_Wait(&req, &status);
		int count;
		MPI_Get_count(&status, PCORD_MPI, &count);
		for (size_t i = 0; i < count; i++)
		{
			particles.push_back({recv_buf[i], false});
		}
	}
}

int main(int argc, char **argv)
{

	MPI_Init(&argc, &argv);
	{
		pcord_t item;
		MPI_Datatype block_types[] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
		int block_lenghts[] = {1, 1, 1, 1};
		MPI_Aint start, displacement[4];

		MPI_Get_address(&item, &start);
		MPI_Get_address(&item.x, &displacement[0]);
		displacement[0] -= start;
		MPI_Get_address(&item.y, &displacement[1]);
		displacement[1] -= start;
		MPI_Get_address(&item.vx, &displacement[2]);
		displacement[2] -= start;
		MPI_Get_address(&item.vy, &displacement[3]);
		displacement[3] -= start;

		MPI_Type_create_struct(4, block_lenghts, displacement, block_types, &PCORD_MPI);
		MPI_Type_commit(&PCORD_MPI);
	}

	int nproc; // number of started MPI processes
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	calculate_dimensions(nproc, &dims[0], &dims[1]);
	MPI_Dims_create(nproc, 2, dims);
	cout0 << "dimensions: " << dims[0] << "x" << dims[1] << " with " << INIT_NO_PARTICLES << " particles\n";

	int periods[2];
	periods[0] = 0;
	periods[1] = 0;
	int reorder = 1;

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &GRID_COMM_MPI);

	int my_coords[2];
	int my_rank;

	MPI_Cart_get(GRID_COMM_MPI, 2, dims, periods, my_coords);
	MPI_Comm_rank(GRID_COMM_MPI, &my_rank);

	int left_rank = grid_rank(my_coords[0] - 1, my_coords[1]);
	int right_rank = grid_rank(my_coords[0] + 1, my_coords[1]);
	int up_rank = grid_rank(my_coords[0], my_coords[1] + 1);
	int down_rank = grid_rank(my_coords[0], my_coords[1] - 1);

	unsigned int time_stamp = 0, time_max;
	float pressure = 0;

	// parse arguments
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s simulation_time\n", argv[0]);
		fprintf(stderr, "For example: %s 10\n", argv[0]);
		exit(1);
	}

	time_max = atoi(argv[1]);

	struct timespec stime, etime;

	/* Initialize */
	// 1. set the walls
	cord_t wall;
	wall.y0 = wall.x0 = 0;
	wall.x1 = BOX_HORIZ_SIZE;
	wall.y1 = BOX_VERT_SIZE;

	// 2. allocate particle buffer and initialize the particles

	std::vector<particle_info> particles{};

	srand(time(NULL) + 1234);

	float r, a;

	int box_size[2] = {(int)BOX_HORIZ_SIZE / dims[0], (int)BOX_VERT_SIZE / dims[1]};

	cord_t box_border;
	box_border.x0 = box_size[0] * my_coords[0];
	box_border.y0 = box_size[1] * my_coords[1];
	box_border.x1 = box_size[0] * (my_coords[0] + 1);
	box_border.y1 = box_size[1] * (my_coords[1] + 1);

	for (int i = 0; i < INIT_NO_PARTICLES / nproc; i++) // We can lose some particles to rounding for now
	{
		// initialize random position
		pcord_t particle{};
		particle.x = wall.x0 + box_border.x0 + rand1() * box_size[0];
		particle.y = wall.y0 + box_border.y0 + rand1() * box_size[1];

		// initialize random velocity
		r = rand1() * MAX_INITIAL_VELOCITY;
		a = rand1() * 2 * PI;
		particle.vx = r * cos(a);
		particle.vy = r * sin(a);
		particles.push_back({particle, false});
	}

	pcord_t *particles_send_left = (pcord_t *)malloc(INIT_NO_PARTICLES * sizeof(pcord_t));
	pcord_t *particles_send_right = (pcord_t *)malloc(INIT_NO_PARTICLES * sizeof(pcord_t));
	pcord_t *particles_send_up = (pcord_t *)malloc(INIT_NO_PARTICLES * sizeof(pcord_t));
	pcord_t *particles_send_down = (pcord_t *)malloc(INIT_NO_PARTICLES * sizeof(pcord_t));

	pcord_t *particles_recv_left = (pcord_t *)malloc(INIT_NO_PARTICLES * sizeof(pcord_t));
	pcord_t *particles_recv_right = (pcord_t *)malloc(INIT_NO_PARTICLES * sizeof(pcord_t));
	pcord_t *particles_recv_up = (pcord_t *)malloc(INIT_NO_PARTICLES * sizeof(pcord_t));
	pcord_t *particles_recv_down = (pcord_t *)malloc(INIT_NO_PARTICLES * sizeof(pcord_t));
	if (my_rank == 0)
	{
		clock_gettime(CLOCK_REALTIME, &stime);
	}

	unsigned sent_over_stat = 0;

	/* Main loop */
	for (time_stamp = 0; time_stamp < time_max; time_stamp++)
	{ // for each time stamp

		for (auto p = particles.begin(); p != particles.end(); ++p)
		{
			{ // for all particles
				if (p->collision)
					continue;

				/* check for collisions */
				for (auto pp = std::next(p); pp != particles.end(); ++pp)
				{
					{
						if (pp->collision)
							continue;
						float t = collide(&p->coords, &pp->coords);
						if (t != -1)
						{ // collision
							p->collision = pp->collision = true;
							interact(&p->coords, &pp->coords, t);
							break; // only check collision of two particles
						}
					}
				}
			}
		}

		// move particles that has not collided with another
		for (particle_info &p : particles)
		{

			if (!p.collision)
			{
				feuler(&p.coords, 1);

				/* check for wall interaction and add the momentum */
				pressure += wall_collide(&p.coords, wall);
			}
			else
			{
				p.collision = false; // reset collision
			}
		}

		unsigned send_up_size = 0, send_left_size = 0, send_right_size = 0, send_down_size = 0;
		MPI_Request req_up, req_left, req_right, req_down;
		MPI_Request req_recv_up, req_recv_left, req_recv_right, req_recv_down;
		for (auto p = particles.begin(); p != particles.end();)
		{
			if (p->coords.x < box_border.x0 && left_rank != -1)
			{
				particles_send_left[send_left_size++] = p->coords;
				p = particles.erase(p);
			}
			else if (p->coords.x > box_border.x1 && right_rank != -1)
			{
				particles_send_right[send_right_size++] = p->coords;
				p = particles.erase(p);
			}
			else if (p->coords.y < box_border.y0 && down_rank != -1)
			{
				particles_send_down[send_down_size++] = p->coords;
				p = particles.erase(p);
			}
			else if (p->coords.y > box_border.y1 && up_rank != -1)
			{
				particles_send_up[send_up_size++] = p->coords;
				p = particles.erase(p);
			}
			else
			{
				++p;
			}
		}
		sent_over_stat += send_up_size + send_down_size + send_left_size + send_right_size;

		if (up_rank != -1)
		{
			MPI_Isend(particles_send_up, send_up_size, PCORD_MPI, up_rank, 0, GRID_COMM_MPI, &req_up);
			MPI_Irecv(particles_recv_up, INIT_NO_PARTICLES, PCORD_MPI, up_rank, 0, GRID_COMM_MPI, &req_recv_up);
		}
		if (down_rank != -1)
		{
			MPI_Isend(particles_send_down, send_down_size, PCORD_MPI, down_rank, 0, GRID_COMM_MPI, &req_down);
			MPI_Irecv(particles_recv_down, INIT_NO_PARTICLES, PCORD_MPI, down_rank, 0, GRID_COMM_MPI, &req_recv_down);
		}
		if (left_rank != -1)
		{
			MPI_Isend(particles_send_left, send_left_size, PCORD_MPI, left_rank, 0, GRID_COMM_MPI, &req_left);
			MPI_Irecv(particles_recv_left, INIT_NO_PARTICLES, PCORD_MPI, left_rank, 0, GRID_COMM_MPI, &req_recv_left);
		}
		if (right_rank != -1)
		{
			MPI_Isend(particles_send_right, send_right_size, PCORD_MPI, right_rank, 0, GRID_COMM_MPI, &req_right);
			MPI_Irecv(particles_recv_right, INIT_NO_PARTICLES, PCORD_MPI, right_rank, 0, GRID_COMM_MPI, &req_recv_right);
		}

		wait_and_push(up_rank, req_recv_up, particles_recv_up, particles);
		wait_and_push(down_rank, req_recv_down, particles_recv_down, particles);
		wait_and_push(left_rank, req_recv_left, particles_recv_left, particles);
		wait_and_push(right_rank, req_recv_right, particles_recv_right, particles);

		MPI_Status tmp_status;

		if (up_rank != -1)
		{
			MPI_Wait(&req_up, &tmp_status);
		}
		if (down_rank != -1)
		{
			MPI_Wait(&req_down, &tmp_status);
		}
		if (left_rank != -1)
		{
			MPI_Wait(&req_left, &tmp_status);
		}
		if (right_rank != -1)
		{
			MPI_Wait(&req_right, &tmp_status);
		}
	}

	float global_pressure = 0;
	MPI_Reduce(&pressure, &global_pressure, 1, MPI_FLOAT, MPI_SUM, 0, GRID_COMM_MPI);

	if (my_rank == 0)
	{
		clock_gettime(CLOCK_REALTIME, &etime);

		double time = (etime.tv_sec - stime.tv_sec) + 1e-9 * (etime.tv_nsec - stime.tv_nsec);
		printf("Time elapsed: %g s\n", time);
		printf("Time per timestep: %g s\n", time / time_max);
		printf("Average pressure = %f\n", global_pressure / (WALL_LENGTH * time_max));
	}

	// std::cout << my_rank << ", sent average: " << sent_over_stat * 1.0f / time_max << "\n";

	MPI_Finalize();
	return 0;
}
