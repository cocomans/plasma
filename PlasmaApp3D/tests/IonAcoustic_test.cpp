#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include <omp.h>

#include "../ParallelInfo.h"
#include "../ImplicitPIC.h"
#include "../ExplicitPIC.h"
#include "../PlasmaData.h"
#include "../Field_Solvers/AmpereSolver.h"
#include "../Field_Solvers/Coupled1D_ES.h"
#include "../Problem_Initializers/IonAcoustic_Initializer.h"
#include "../Problem_Initializers/TwoStream_Initializer.h"



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
	int nprocs = 1;

	nprocs = 1;

	int nptcls_node,nptcls_total;
	int nthreads = ndevices + nprocs;

	srand(1209876*(5+myid));

	system("echo $LD_LIBRARY_PATH");

	// Setup plasma data
	PlasmaData pdata(argc,argv);

	pdata.ny = 2;
	pdata.nz = 2;

	pdata.setup();

	pdata.ndimensions = 1;
	pdata.nVelocity = 1;
	pdata.Bmag_avg = 0;

	pdata.SimName = "IonAcousticShock";

	nptcls_node = ndevices*pdata.nptcls_gpu + nprocs * pdata.nptcls_cpu;

	MPI_Allreduce(&nptcls_node,&nptcls_total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	pdata.nptcls_total = nptcls_total;

	ParallelInfo** pinfo = (ParallelInfo**)malloc((nthreads+1)*sizeof(ParallelInfo*));

	pinfo[0] = new ParallelInfo;
	pinfo[0] -> nthreads = nthreads;
	pinfo[0] -> n_nodes = num_nodes;
	pinfo[0] -> myid_mpi = myid;
	pinfo[0] -> device_type = 0;

	Coupled1D_ES* LOsolver = new Coupled1D_ES;
//	AmpereSolver* LOsolver = new AmpereSolver;
	IonAcoustic_Initializer* initializer = new IonAcoustic_Initializer(&pdata);

	ImplicitPIC* simulation = new ImplicitPIC(&pdata,LOsolver,initializer,pinfo,ndevices);

	simulation->simulate();







	MPI_Finalize();
}
