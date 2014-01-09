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
#include "../RunData.h"
//#include "../Field_Solvers/Coupled1D_ES.h"
//#include "../Field_Solvers/Coupled1D_ES_JFNK.h"
#include "../Field_Solvers/AmpereSolver.h"
#include "../Problem_Initializers/TwoStream_Initializer.h"
#include "../Field_Solvers/HoLoInterface/HoLoInterface.h"
#include "../Field_Solvers/HoLoInterface2D/HoLoInterface2D.h"
#include "../Util/webui/sse-server.h"
#include "../Util/webui/sseDataStream.h"
#include "../OutputCtrl.h"



int main(int argc,char** argv)
{
	int rc,myid,num_nodes;


	system("numactl --show");


	int provided;
	rc = MPI_Init(&argc,&argv);

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
//	iControlServer = 1;
	PlasmaData pdata(argc,argv);

	sseServer* ControlServer= pdata.output_ctrl->server;
	WebControlSystem* control= pdata.output_ctrl->control;

//	// Setup the input parser
//	InputParser* parser = new InputParser(argc,argv);
//	sseDataStream* stream = new sseDataStream(&pdata,4,0);
//	int exit_control = 0;
//
//
//	char* ROOT;
//	int PORTi;
//	parser->GetParam("SERVER_PORT_NUMBER","-p",PORTi,32000);
//
//	std::stringstream tmp;
//	tmp << PORTi;
//	char PORTc[6];
//
//	sprintf(PORTc,"%s",tmp.str().c_str());
//	WebControlSystem* control = new WebControlSystem(parser);
//	ControlServer = new sseServer(ROOT,PORTc,stream,control,&exit_control);


	while(!(control->CheckiRun()))
	{
		sleep(1);
	}

	MPI_Finalize();
}
