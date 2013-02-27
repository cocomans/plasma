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
#include "../Field_Solvers/Coupled1D_ES.h"
#include "../Field_Solvers/AmpereSolver.h"
#include "../Problem_Initializers/TwoStream_Initializer.h"




int main(int argc,char** argv)
{
	int rc,myid,num_nodes;




	rc = MPI_Init(&argc,&argv);

	setvbuf( stdout, NULL, _IONBF, 0 );

	if (rc != MPI_SUCCESS)
	{
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}



	MPI_Comm_size(MPI_COMM_WORLD,&num_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	MPI_Barrier(MPI_COMM_WORLD);



	int ndevices = 0;
	int nprocs = 1;

	//nprocs = 1;

	int nptcls_node,nptcls_total;
	int nthreads = ndevices + nprocs;

	srand(156*(32+5*myid));



	// Setup plasma data
	PlasmaData pdata(argc,argv);

	pdata.ny = 2;
	pdata.nz = 2;

	// This is a 1D, Electrostatic Problem
	pdata.ndimensions = 1;
	pdata.nVelocity = 1;
	pdata.iEM = 0;
	pdata.Bmag_avg = 0.0;

	pdata.setup();

	//nptcls_node = ndevices*pdata.nptcls_gpu + nprocs * pdata.nptcls_cpu;

	//MPI_Allreduce(&nptcls_node,&nptcls_total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	pdata.nptcls_total = pdata.nptcls*num_nodes;

	ParallelInfo** pinfo = (ParallelInfo**)malloc((nthreads+1)*sizeof(ParallelInfo*));

	pinfo[0] = new ParallelInfo;
	pinfo[0] -> nthreads = nthreads;
	pinfo[0] -> n_nodes = num_nodes;
	pinfo[0] -> myid_mpi = myid;
	pinfo[0] -> device_type = 0;

	Coupled1D_ES* LOsolver = new Coupled1D_ES;
//	AmpereSolver* LOsolver = new AmpereSolver;
	TwoStream_Initializer* initializer = new TwoStream_Initializer(&pdata);



	ImplicitPIC* simulation = new ImplicitPIC(&pdata,LOsolver,initializer,pinfo,ndevices);

	simulation->simulate();







	MPI_Finalize();
}
