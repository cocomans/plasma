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
#include "../Field_Solvers/ConstantSolver.h"
#include "../Field_Solvers/HoLoInterface2D/HoLoInterface2D.h"
#include "../Problem_Initializers/EWave2D_Initializer.h"
#include "../Normalizations/NormIslandEM.h"




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



	Normalization* mynorms = new NormIslandEM();

	// Setup plasma data
	PlasmaData pdata(argc,argv,mynorms);
	pdata.nVelocity = 3;
	pdata.ndimensions = 2;
	pdata.iEM = 1;

	pdata.nz = 2;


	pdata.setup();



	HoLoInterface2D* LOsolver = new HoLoInterface2D(&pdata);
	EWave2D_Initializer* initializer = new EWave2D_Initializer(&pdata);


	ImplicitPIC* simulation = new ImplicitPIC(&pdata,LOsolver,initializer);

	simulation->simulate();







	MPI_Finalize();
}
