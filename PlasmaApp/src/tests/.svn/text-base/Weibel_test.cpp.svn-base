#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "../cuda_defines.h"
#include <gnuplot_i.h>
//#include <mpi.h>
#include <omp.h>

#include "../ParallelInfo.h"
#include "../ImplicitPIC.h"
#include "../PlasmaData.h"
//#include "../Field_Solvers/Coupled1D_ES.h"
//#include "../Field_Solvers/Coupled1D_ES_JFNK.h"
#include "../Field_Solvers/AmpereSolver.h"
#include "../Problem_Initializers/Weibel_Initializer.h"
#include "../Field_Solvers/HoLoInterface/HoLoInterface.h"
#include "../Field_Solvers/HoLoInterface2D/HoLoInterface2D.h"

#include "../Normalizations/NormElectronEM.h"



int main(int argc,char** argv)
{
	int rc,myid,num_nodes;

	system("numactl --show");

	rc = MPI_Init(&argc,&argv);
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);


	setvbuf( stdout, NULL, _IONBF, 0 );

	if (rc != MPI_SUCCESS)
	{
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	system("numactl --show");

	printf("made it past system\n");

	MPI_Comm_size(MPI_COMM_WORLD,&num_nodes);
	printf("made it past com size\n");

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	printf("made it past com rank\n");

	MPI_Barrier(MPI_COMM_WORLD);







	int ndevices = 0;
	int nprocs = 1;

	//nprocs = 1;

	int nptcls_node,nptcls_total;
	int nthreads = ndevices + nprocs;

	srand(156*(32+5*myid));

	Normalization* mynorms = new NormElectronEM();


	// Setup plasma data
	PlasmaData pdata(argc,argv,mynorms);

	pdata.ny = 2;
	pdata.nz = 2;

	// This is a 1D, Electrostatic Problem
	pdata.ndimensions = 1;
	pdata.nVelocity = 3;
	pdata.iEM = 1;

//	pdata.Lx = 1.0;

	pdata.setup();

	//nptcls_node = ndevices*pdata.nptcls_gpu + nprocs * pdata.nptcls_cpu;

	//MPI_Allreduce(&nptcls_node,&nptcls_total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);



//	Coupled1D_ES* LOsolver = new Coupled1D_ES;
//	Coupled1D_ES_JFNK* LOsolver = new Coupled1D_ES_JFNK;
//	AmpereSolver* LOsolver = new AmpereSolver;
	HoLoInterface* LOsolver = new HoLoInterface(&pdata);
	Weibel_Initializer* initializer = new Weibel_Initializer(&pdata);



	ImplicitPIC* simulation = new ImplicitPIC(&pdata,LOsolver,initializer);

	simulation->simulate();







	MPI_Finalize();
}
