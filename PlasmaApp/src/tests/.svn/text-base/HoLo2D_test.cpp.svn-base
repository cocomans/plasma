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

	double alpha = 0.01;

	//nptcls_node = ndevices*pdata.nptcls_gpu + nprocs * pdata.nptcls_cpu;

	//MPI_Allreduce(&nptcls_node,&nptcls_total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);



//	Coupled1D_ES* LOsolver = new Coupled1D_ES;
//	Coupled1D_ES_JFNK* LOsolver = new Coupled1D_ES_JFNK;
//	AmpereSolver* LOsolver = new AmpereSolver;
	HoLoInterface2D* LOsolver = new HoLoInterface2D(&pdata);
	Weibel_Initializer* initializer = new Weibel_Initializer(&pdata);
	OutputCtrl* output_ctrl = pdata.output_ctrl;

	Teuchos::RCP<EpVecWrapper> X;
	Teuchos::RCP<EpVecWrapper> Xold;

	X = LOsolver -> X;
	Xold = LOsolver -> Xold;

	double kx = 2.0*pi_const/pdata.Lx;
	double Te = 1.0;
	ImplicitPIC* sim = new ImplicitPIC(&pdata,LOsolver,initializer);

	using namespace EM2D3V;
	for(int k=0;k<pdata.ny;k++)
	for(int j=0;j<pdata.ny;j++)
		for(int i=0;i<pdata.nx;i++)
		{
			double xc = (i+0.5)*pdata.dxdi + pdata.xmin;
			double xf = (i)*pdata.dxdi + pdata.xmin;

			double yc = (j+0.5)*pdata.dydi + pdata.ymin;
			double yf = j*pdata.dydi + pdata.ymin;

			sim->fields_next->getPhi(i,j,k) = -alpha*cos(kx*xc)/(kx*kx);
			sim->fields_next->getE(i,j,k,0) = alpha*sin(kx*xf)/(kx);
			sim->fields_next->getE(i,j,k,1) = 0;
			sim->fields_next->getE(i,j,k,2) = 0;



			sim->moments->moments_next->get_val(i,j,k,0,HOMoments_charge) = alpha*cos(kx*xc) + 1.0;
			sim->moments->moments_next->get_val(i,j,k,1,HOMoments_charge) = 1.0;

			X->Get(px_e,i,j) = 0.0;
			X->Get(py_e,i,j) = 0.0;
			X->Get(pz_e,i,j) = 0.0;

			X->Get(px_i,i,j) = 0.0;
			X->Get(py_i,i,j) = 0.0;
			X->Get(pz_i,i,j) = 0.0;

			sim->moments->moments_next->get_val(i,j,k,0,HOMoments_S2xx) = 2.0*Te/(3.0*pdata.mspecies[0]);
			sim->moments->moments_next->get_val(i,j,k,0,HOMoments_S2yy) = 2.0*Te/(3.0*pdata.mspecies[0]);
			sim->moments->moments_next->get_val(i,j,k,0,HOMoments_S2zz) = 2.0*Te/(3.0*pdata.mspecies[0]);

			sim->moments->moments_next->get_val(i,j,k,1,HOMoments_S2xx) = 2.0*Te/(3.0*pdata.mspecies[1]);
			sim->moments->moments_next->get_val(i,j,k,1,HOMoments_S2yy) = 2.0*Te/(3.0*pdata.mspecies[1]);
			sim->moments->moments_next->get_val(i,j,k,1,HOMoments_S2zz) = 2.0*Te/(3.0*pdata.mspecies[1]);




			sim->fields_next->getChi(i,j,k) = 0.0;
			sim->fields_next->getA(i,j,k,0) = 0.0;
			sim->fields_next->getA(i,j,k,1) = 0.0;
			sim->fields_next->getA(i,j,k,2) = 0.0;
			sim->fields_next->getB(i,j,k,0) = 0.0;
			sim->fields_next->getB(i,j,k,1) = 0.0;
			sim->fields_next->getB(i,j,k,2) = 0.0;









		}


	for(int i=0;i<pdata.nsteps;i++)
	{

	}







	MPI_Finalize();
}
