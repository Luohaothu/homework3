#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <xmmintrin.h>
#include <mpi.h>
#include <omp.h>
#define stepnum 1000
#define x_size (2 << 7)
#define y_size (2 << 7)
#define z_size (2 << 7)
#define block_y_size (2 << 3)
#define block_z_size (2 << 7)
#define FINISHED_SIGNAL 1
#define min(x,y) (x < y? x : y)
#define GRID0 0
#define GRID1 1
void compute(double *grid0, double *grid1, int rank, int size);
double* creat_grid(int x_s, int y_s, int z_s, int size, int rank, int identifyer);
inline long index(long y_s, long z_s, long i, long j, long k);
void block(double* grid0, double* grid1, long x_s, long y_s, long z_s, int b_y_s, int b_z_s);
inline long index(long y_s, long z_s, long i, long j, long k);


//main program
int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int size, rank;
	long t1, t2;
	static int ranks[1] = { 0 };
	MPI_Request request1, request2, request3, request4;
	MPI_Status status, status1, status2, status3, status4;
	MPI_Group MPI_GROUP_WORLD, grprem;
	MPI_Comm commslave, newcomm;
	MPI_Comm commsp;
	MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
	MPI_Group_excl(MPI_GROUP_WORLD, 1, ranks, &grprem);
	MPI_Comm_create(MPI_COMM_WORLD, grprem, &commslave);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("Node %d in %d is ready\n", rank, size);
	//Initialize
	double *grid0 = creat_grid(x_size, y_size, z_size, size, rank, GRID0);
	double *grid1 = creat_grid(x_size, y_size, z_size, size, rank, GRID1);
	MPI_Barrier(MPI_COMM_WORLD);
	if (size != 1)
	{
		if (rank == 0)
		{
			for (int i = 1; i < size; i++)
			{
				int len = (i == size - 1) ? (x_size / (size - 1)) : (x_size / (size - 1) + x_size % (size - 1));
				MPI_Ssend(grid0 + (i - 1)*x_size / (size - 1), len*y_size*z_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
			}
		}
		else
		{
			for (int i = 1; i < size; i++)
			{
				if (rank == i)
				{
					int len = (i == size - 1) ? (x_size / (size - 1)) : (x_size / (size - 1) + x_size % (size - 1));
					MPI_Recv(grid0 + y_size*z_size, len*y_size*z_size, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &status);
				}
			}

		}
	}
	//Compute
	if (rank == 0) printf("Start computing...\n");
	if (rank != 0 && size > 1)
	{
		for (int t = 0; t < stepnum; t++)
		{
			compute(grid0, grid1, rank, size);
			//send right slice of data to next node, then receieve right slice of data form next node
			if (rank < size - 1)
			{
			MPI_Isend(grid1 + (1 + x_size / (size - 1))*y_size*z_size,
			y_size*z_size, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD, &request1);
			MPI_Irecv(grid1 + (1 + x_size / (size - 1))*y_size*z_size,
			y_size*z_size, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &request2);
			MPI_Wait(&request1, &status1);
			MPI_Wait(&request2, &status2);
			}
			//receieve left slice of data from perior node, then send left slice of data to perior node
			if (rank > 1)
			{
			MPI_Irecv(grid1, y_size*z_size, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &request3);
			MPI_Isend(grid1, y_size*z_size, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD, &request4);
			MPI_Wait(&request3, &status3);
			MPI_Wait(&request4, &status4);
			}
			double *temp;
			temp = grid0;
			grid0 = grid1;
			grid1 = temp;
			MPI_Barrier(commslave);
		}
	}
	else if (size == 1)
	{
		for (int t = 0; t < stepnum; t++)
		{
			compute(grid0, grid1, rank, size);
			double *temp;
			temp = grid0;
			grid0 = grid1;
			grid1 = temp;
		}
	}
	else{ ; }
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Rank %d finished computing!\n", rank);
	//Gather data form nodes to host
	if (size != 1)
	{
		if (stepnum % 2)
		{
			double *temp;
			temp = grid0;
			grid0 = grid1;
			grid1 = temp;
		}
		for (int i = 1; i < size; i++)
		{
			if (rank == i)
			{
				int len = (i == size - 1) ? (x_size / (size - 1)) : (x_size / (size - 1) + x_size % (size - 1));
				MPI_Isend(grid0 + y_size*z_size, len*y_size*z_size, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &request2);
				MPI_Wait(&request2, &status2);
			}
		}
		if (rank == 0)
		{
			for (int i = 1; i < size; i++)
			{
				int len = (i == size - 1) ? (x_size / (size - 1)) : (x_size / (size - 1) + x_size % (size - 1));
				MPI_Irecv(grid1, len*y_size*z_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request1);
				MPI_Wait(&request1, &status1);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0) printf("All work complete\n");
	MPI_Finalize();
	return 0;
}

double* creat_grid(int x_s, int y_s, int z_s, int size, int rank, int identifyer)
{
	if (size == 1)
	{
		double *grid0 = (double*)_mm_malloc(x_s*y_s*z_s*sizeof(double), 64);
		double *grid1 = (double*)_mm_malloc(x_s*y_s*z_s*sizeof(double), 64);
		memset(grid0, 0, x_s*y_s*z_s*sizeof(double));
		memset(grid1, 0, x_s*y_s*z_s*sizeof(double));
		printf("Grid%d is created for node %d\n", identifyer, rank);
		return (identifyer) ? (grid1) : (grid0);
	}
	else if (rank > 0)
	{
		int len = (rank == size - 1) ? (x_s / (size - 1)) : (x_s / (size - 1) + x_s % (size - 1));
		double *grid = (double*)_mm_malloc((len + 2)*y_s*z_s*sizeof(double), 64);
		memset(grid, 0, (len + 2)*y_s*z_s*sizeof(double));
		printf("Grid%d is created for node %d\n", identifyer, rank);
		return grid;
	}
	else
	{
		double *grid = (double*)_mm_malloc(x_s*y_s*z_s*sizeof(double), 64);
		memset(grid, 0, x_s*y_s*z_s*sizeof(double));
		printf("Grid%d is created for node %d\n", identifyer, rank);
		return grid;
	}
}

void compute(double *grid0, double *grid1, int rank, int size)
{
	long x_s, y_s, z_s;
	int b_y_s, b_z_s;

	if (size == 1)
	{
		x_s = x_size;
		y_s = y_size;
		z_s = z_size;
	}
	else if (rank > 0)
	{
		int len = (rank == size - 1) ? (x_size / (size - 1)) : (x_size / (size - 1) + x_size % (size - 1));
		x_s = len + 2;
		y_s = y_size;
		z_s = z_size;
	}
	else
	{
		x_s = x_size;
		y_s = y_size;
		z_s = z_size;
	}
	b_y_s = y_s / block_y_size;
	b_z_s = z_s / block_z_size;
	block(grid0, grid1, x_s, y_s, z_s, b_y_s, b_z_s);

}

void block(double* grid0, double* grid1, long x_s, long y_s, long z_s, int b_y_s, int b_z_s)
{
	register double one_six = 1 / 6;
#pragma vector aligned
#pragma omp parallel for num_threads(omp_get_num_procs()) schedule (guided)
	for (long b_j = 1; b_j < y_s - 1; b_j += b_y_s)
	{
		for (long b_k = 1; b_k < z_s - 1; b_k += b_z_s)
		{
			for (long i = 1; i < x_s - 1; i++)
			{
				for (long j = b_j; j < min(b_j + b_y_s, y_s - 1); j++)
				{
					for (long k = b_k; k < min(b_k + b_z_s, z_s - 1); k++)
					{
						grid1[index(y_s, z_s, i, j, k)] = (
							grid0[index(y_s, z_s, i, j, k - 1)] + grid0[index(y_s, z_s, i, j, k + 1)] +
							grid0[index(y_s, z_s, i, j - 1, k)] + grid0[index(y_s, z_s, i, j + 1, k)] +
							grid0[index(y_s, z_s, i - 1, j, k)] + grid0[index(y_s, z_s, i + 1, j, k)]) * one_six;
					}
				}
			}
		}
	}
}

inline long index(long y_s, long z_s, long i, long j, long k)
{
	return i * y_s * z_s + j * z_s + k;
}
