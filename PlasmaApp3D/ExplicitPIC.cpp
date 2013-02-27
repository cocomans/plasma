#include "ExplicitPIC.h"
#include "HOMoments.h"
#include "PlasmaData.h"
#include "ParticleList.h"
#include "ParticleListCPU.h"
#include "ParticleListGPU.cuh"
#include "FieldDataCPU.h"
#include "FieldDataGPU.cuh"
#include "ProblemInitializer.h"
#include "FieldSolver.h"
#include "ParallelInfo.h"


ExplicitPIC::ExplicitPIC(PlasmaData* pdata_in, FieldSolver* LOsolver_in, ProblemInitializer* initializer_in,
						ParallelInfo** pinfo,int ndevices_in)
{
	nthreads = pinfo[0]->nthreads;
	myid_mpi = pinfo[0]->myid_mpi;
	n_nodes = pinfo[0]->n_nodes;
	ndevices = ndevices_in;

	pll_info = pinfo;
	pdata = pdata_in;

	LOsolver = LOsolver_in;
	initializer = initializer_in;


	// Allocate space for object pointers
	fields = (FieldData**)malloc((nthreads+1)*sizeof(FieldData*));
	particles = (ParticleList**)malloc(nthreads*sizeof(ParticleList*));
	moments = (HOMoments**)malloc((nthreads+1)*sizeof(HOMoments*));

	fields[0] = new FieldDataCPU;
	fields[0] -> allocate(pdata);

	moments[0] = new HOMoments(pdata);
}

void ExplicitPIC::simulate()
{

	// Set the number of OMP threads
	//omp_set_num_threads(nthreads);



	int tid = 0;
	int nptcls;
	int numprocs = omp_get_num_procs();

	//printf("nthreads = %i with tid = %i\n",nthreads,tid);

	pll_info[tid+1] = new ParallelInfo;
	*(pll_info[tid+1]) = *(pll_info[0]);
	ParallelInfo* myinfo = pll_info[tid+1];

	myinfo->tid = tid;

	if(tid < ndevices)
		myinfo->device_type = 1;
	else
		myinfo->device_type = 0;


	// Allocate Particle Lists, HOMoments, and Fields
	if(myinfo->device_type == 0)
	{
		// Decvice is a CPU
		particles[tid] = new ParticleListCPU();
		fields[tid+1] = new FieldDataCPU();

		moments[tid+1] = (HOMoments*)malloc(numprocs*sizeof(HOMoments));
		for(int i=0;i<numprocs;i++)
		{
			moments[tid+1][i] = *(new HOMoments(pdata));
		}

		nptcls = pdata->nptcls_cpu;

	}
	else if(myinfo->device_type == 1)
	{

#ifndef NO_CUDA
		// Device is a GPU
		particles[tid] = new ParticleListGPU();
		fields[tid+1] = new FieldDataGPU();
		moments[tid+1] = new HOMoments(pdata);

		nptcls = pdata->nptcls_gpu;
#endif
	}

	printf("Allocating Particles\n");
	particles[tid] -> allocate(pdata,nptcls);
	printf("Allocating Fields\n");
	fields[tid+1] -> allocate(pdata);




	// Initialize fields and particles
	printf("Initializing Particles\n");
	initializer -> initialize_particles(particles[tid],moments[tid+1],myinfo);

	// Reduce Moments
	printf("Reducing Moments\n");
	moments[tid+1] -> mpi_reduce(moments,myinfo);

	// Initialize Fields
	printf("Initializing Fields\n");
	initializer -> initialize_fields(fields,moments,myinfo);

	moments[0] -> init_plot();
	fields[0] -> init_plot();

	if(myinfo->myid_mpi == 0){
	moments[0] -> plot(pdata->nz/2,0,0,HOMoments_currentx);
	}
	getchar();

	pdata->time_done = 0;

	// Main time loop
	for(int istep=0;istep<pdata->nsteps;istep++)
	{


		// Move the particles
		//printf("Pushing Particles\n");
		particles[tid] -> push(pdata,fields[tid+1],moments[tid+1]);

		// Reduce Moment quantities
		moments[tid+1] -> mpi_reduce(moments,myinfo);
		if(myinfo->myid_mpi == 0){
			moments[0] -> reset_plot();
			moments[0] -> plot(pdata->nz/2,0,0,HOMoments_currentx);
		}

		// Field Solve
		if(myinfo->tid == 0)
		{
			if(myinfo->myid_mpi == 0)
			{

				LOsolver->solve(pdata,fields[0],moments[0]);
			}
		}


		// Broadcast Fields
		fields[tid+1] -> broadcast(fields,myinfo);

		if(myinfo->myid_mpi == 0){
			fields[0] -> reset_plot();
			fields[0] -> plot(pdata,pdata->nz/2,0,0,0);
		}

		// Plot Fields and Moments
		initializer->check_step(particles[tid],moments,moments_old,fields,myinfo);

		pdata->time_done += pdata->dt;
	}

	getchar();

	moments[0]->close_plot();
	fields[0] -> close_plot();







}
