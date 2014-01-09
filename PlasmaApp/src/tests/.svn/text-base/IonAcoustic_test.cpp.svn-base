#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "../cuda_defines.h"
#include <gnuplot_i.h>
#include <mpi.h>
#include <omp.h>

#include "../ParallelInfo.h"
#include "../ImplicitPIC.h"
#include "../PlasmaData.h"
#include "../RunData.h"
#include "../Field_Solvers/AmpereSolver.h"
//#include "../Field_Solvers/Coupled1D_ES.h"
#include "../Problem_Initializers/IonAcoustic_Initializer.h"
#include "../Problem_Initializers/TwoStream_Initializer.h"
#include "../Field_Solvers/HoLoInterface/HoLoInterface.h"




int main(int argc,char* argv[])
{
	int rc,myid,num_nodes;

	system("numactl --show");

	int provided;

	rc = MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);




	if (rc != MPI_SUCCESS)
	{
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_size(MPI_COMM_WORLD,&num_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	system("numactl --show");


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
	pdata.iEM = 0;

	pdata.rdata->SimName = "IonAcousticShock";


//	Coupled1D_ES* LOsolver = new Coupled1D_ES;
//	AmpereSolver* LOsolver = new AmpereSolver;
	HoLoInterface* LOsolver = new HoLoInterface(&pdata);
	IonAcoustic_Initializer* initializer = new IonAcoustic_Initializer(&pdata);

	ImplicitPIC* simulation = new ImplicitPIC(&pdata,LOsolver,initializer);

	simulation->simulate();







	MPI_Finalize();
}
