#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include <omp.h>

#include "../ParallelInfo.h"
#include "../ExplicitPIC.h"
#include "../ImplicitPIC.h"
#include "../PlasmaData.h"
#include "../Field_Solvers/ConstantSolver.h"
#include "../Problem_Initializers/ExB_Initializer.h"




int main(int argc,char* argv[])
{
	int rc,myid,num_nodes;

	rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS)
	{
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_size(MPI_COMM_WORLD,&num_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	int ndevices = 0;
	int nprocs = omp_get_num_procs();

	nprocs = 1;

	int nptcls_node,nptcls_total;
	int nthreads = ndevices + nprocs;




	// Setup plasma data
	PlasmaData pdata(argc,argv);
	pdata.nVelocity = 3;
	pdata.ndimensions = 2;
	pdata.iEM = 1;

	pdata.nz = 2;

	pdata.setup();

	nptcls_node = ndevices*pdata.nptcls_gpu + nprocs * pdata.nptcls_cpu;

	MPI_Allreduce(&nptcls_node,&nptcls_total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	pdata.nptcls_total = nptcls_total;

	ParallelInfo** pinfo = (ParallelInfo**)malloc((nthreads+1)*sizeof(ParallelInfo*));

	pinfo[0] = new ParallelInfo;
	pinfo[0] -> nthreads = nthreads;
	pinfo[0] -> n_nodes = num_nodes;
	pinfo[0] -> myid_mpi = myid;
	pinfo[0] -> device_type = 0;

	ConstantSolver* LOsolver = new ConstantSolver;
	ExB_Initializer* initializer = new ExB_Initializer(&pdata);


	ImplicitPIC* simulation = new ImplicitPIC(&pdata,LOsolver,initializer,pinfo,ndevices);

	simulation->simulate();







	MPI_Finalize();
}
