#ifndef PARALLEL_INFO_H
#define PARALLEL_INFO_H


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>

class ParallelInfo
{
public:

	int nthreads; // number of OpenMP threads
	int n_nodes; // number of MPI Nodes
	int tid; // OpenMP Thread ID
	int myid_mpi; // MPI id
	int device_type;



};

#endif /* PARALLEL_INFO_H */
