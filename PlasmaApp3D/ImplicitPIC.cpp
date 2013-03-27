/*-------------------------------------------------------------------------*/
/**
	@file		ImplicitPIC.cpp
	@author	J. Payne
	@date		12/21/2012


*/
/*--------------------------------------------------------------------------*/
#include "ImplicitPIC.h"
#include "HOMoments.h"
#include "PlasmaData.h"
#include "ParticleList.h"
#include "ParticleListCPU.h"
#include "ParticleListCPUSorted.h"
#include "ParticleListGPU.cuh"
#include "ParticleListGPUSimple.cuh"
#include "ParticleListMIC.h"
#include "FieldDataCPU.h"
#include "FieldDataCPU2D.h"
#include "FieldDataGPU.cuh"
#include "ProblemInitializer.h"
#include "FieldSolver.h"
#include "ParallelInfo.h"
#include <gnuplot_i.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "Util/mkpath.h"
#include "Util/OutputPage.h"
#include "ShapeFunctions.h"




ImplicitPIC::ImplicitPIC(PlasmaData* pdata_in, FieldSolver* LOsolver_in, ProblemInitializer* initializer_in,
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
	fields_half = (FieldData**)malloc((nthreads+1)*sizeof(FieldData*));
	particles_old = (ParticleList**)malloc(nthreads*sizeof(ParticleList*));
	particles_next = (ParticleList**)malloc(nthreads*sizeof(ParticleList*));
	moments = (HOMoments**)malloc((nthreads+1)*sizeof(HOMoments*));
	if(pdata->ndimensions == 2)
	{
		fields_half[0] = new FieldDataCPU2D();
		fields_next = new FieldDataCPU2D();
		fields_old = new FieldDataCPU2D();
	}

	else
	{
		fields_half[0] = new FieldDataCPU();
		fields_next = new FieldDataCPU();
		fields_old = new FieldDataCPU();
	}

	fields_half[0] -> allocate(pdata);
	fields_next -> allocate(pdata);
	fields_old -> allocate(pdata);

	moments[0] = new HOMoments(pdata);
	moments_old = new HOMoments(pdata);

}

/*-------------------------------------------------------------------------*/
/**
	@brief Perform an Implicit PIC simulation

	In this method ParticleLists, FieldData, and HOMoments are chosen based on
	device type, allocated, initialized, and advanced in a time step loop.
	Currently this method uses a Crank-Nicholson discritization and a Piccard
	iteration scheme to converege the non-linear interactino between the HO and
	LO systems for every time step.

	Once all time steps have been completed results are saved and timing
	data is both printed and saved.
*/
/*--------------------------------------------------------------------------*/
void ImplicitPIC::simulate()
{

	// Set the number of OMP threads




	int tid = 0;
	int nptcls;
	int numprocs;

	CPUSpeed myspeed;

	pll_info[tid+1] = new ParallelInfo;
	numprocs = pdata->num_cores;
	pll_info[0]->nthreads = numprocs;
	*(pll_info[tid+1]) = *(pll_info[0]);
	ParallelInfo* myinfo = pll_info[tid+1];
	omp_set_num_threads(numprocs);


	setvbuf( stdout, NULL, _IONBF, 0 );
	if(myinfo->myid_mpi == 0) printf("nthreads = %i with tid = %i\n",numprocs,tid);



	myinfo->tid = tid;

	myinfo->device_type = pdata->device_type;


	// Allocate Particle Lists, HOMoments, and Fields
	if(myinfo->device_type == 0)
	{
		// Decvice is a CPU

		switch(pdata->iParticleListCPU)
		{
		case 0:
			particles_old[tid] = new ParticleListMIC();
			particles_next[tid] = new ParticleListMIC();
			break;
		case 1:
			particles_old[tid] = new ParticleListCPUSorted();
			particles_next[tid] = new ParticleListCPUSorted();
			break;
		default:
			particles_old[tid] = new ParticleListCPU();
			particles_next[tid] = new ParticleListCPU();
			break;
		}

		if(pdata->ndimensions == 2)
		{
			fields_half[tid+1] = new FieldDataCPU2D();
		}

		else
		{
			fields_half[tid+1] = new FieldDataCPU();
		}


		moments[tid+1] = (HOMoments*)malloc(numprocs*sizeof(HOMoments));
		for(int i=0;i<numprocs;i++)
		{
			moments[tid+1][i] = *(new HOMoments(pdata));
		}

	}
	else if(myinfo->device_type == 1)
	{

#ifndef NO_CUDA
		myinfo->nthreads = 1;
		// Device is a GPU
		particles_old[tid] = new ParticleListGPU();
		particles_next[tid] = new ParticleListGPU();
		fields_half[tid+1] = new FieldDataCPU2D();
		moments[tid+1] = new HOMoments(pdata);
#endif

	}

	if(myinfo->myid_mpi == 0) printf("Allocating Particles\n");

	particles_old[tid] -> allocate(pdata,pdata->my_nptcls);
	particles_next[tid] -> allocate(pdata,pdata->my_nptcls);

	// since we don't care about resetting these values from picard to picard
	particles_next[tid] -> num_subcycles = particles_old[tid] -> num_subcycles;
	particles_next[tid] -> num_piccard = particles_old[tid] -> num_piccard;
	particles_next[tid] -> num_piccard2 = particles_old[tid] -> num_piccard2;
	particles_next[tid] -> piccard_timer = particles_old[tid] -> piccard_timer;
	particles_next[tid] -> accel_timer = particles_old[tid] -> accel_timer;
	particles_next[tid] -> tally_timer = particles_old[tid] -> tally_timer;
	particles_next[tid] -> crossing_timer = particles_old[tid] -> crossing_timer;
	particles_next[tid] -> dtau_est_timer = particles_old[tid] -> dtau_est_timer;
	particles_next[tid] -> tally_timer2 = particles_old[tid] -> tally_timer2;
	particles_next[tid] -> load_store_timer = particles_old[tid] -> load_store_timer;
	particles_next[tid] -> push_timer = particles_old[tid] -> push_timer;

//	if(myinfo->myid_mpi == 0)
//	{
		if(pdata->plot_flag)particles_old[tid] -> plot = gnuplot_init();

		particles_next[tid] -> plot = particles_old[tid] -> plot;
//	}
	if(myinfo->myid_mpi == 0) printf("Allocating Fields\n");

	fields_half[tid+1] -> allocate(pdata);

	MPI_Barrier(MPI_COMM_WORLD);


	// Initialize fields and particles
	if(myinfo->myid_mpi == 0) printf("Initializing Particles\n");
	MPI_Barrier(MPI_COMM_WORLD);
	initializer -> initialize_particles(particles_old[tid],moments[tid+1],myinfo);

	MPI_Barrier(MPI_COMM_WORLD);
	// Reduce Moments
	if(myinfo->myid_mpi == 0) printf("Reducing Moments\n");
	MPI_Barrier(MPI_COMM_WORLD);
	moments[tid+1] -> mpi_reduce(moments,myinfo);

	if(myinfo->myid_mpi == 0) printf("Applying Weights\n");
	if(tid == 0)
	{
		if((myinfo->myid_mpi == 0)||(pdata->lo_all))
		{
			moments[0]->apply_weights();
		}
	}

	// Initialize Fields
	if(myinfo->myid_mpi == 0) printf("Initializing Fields\n");
	initializer -> initialize_fields(fields_half,moments,myinfo);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myinfo->myid_mpi == 0) printf("Broadcasting Fields\n");
	fields_half[tid+1] -> broadcast(fields_half,myinfo);

	MPI_Barrier(MPI_COMM_WORLD);
	HOMoments charge_plot;


	if(tid == 0)
	{
		if((myinfo->myid_mpi == 0)||(pdata->lo_all))
		{
			if(myinfo->myid_mpi == 0) printf("Copying Fields\n");
			fields_old -> copy_from(fields_half[0]);

			moments_old -> copy_from(moments[0]);
			if(myinfo->myid_mpi == 0) printf("Init LO\n");
			LOsolver->init(pdata,fields_old,fields_half[0],fields_next,moments_old,moments[0]);

			fields_next -> copy_from(fields_old);

			LOsolver->solve(pdata,fields_next,moments[0]);
		}
	}


	// Plotting Stuff
	if((myinfo->myid_mpi == 0)&&(pdata->plot_flag)){
		printf("Plotting Stuff\n");
		moments[tid+1] -> init_plot();
		fields_half[0] -> init_plot();
		charge_plot = *(moments[0]);
		moments[0] -> init_plot();
		fields_old -> init_plot();
		fields_next -> init_plot();
		charge_plot.init_plot();
		moments[0] -> plot(pdata->nz/2,0,0,HOMoments_currentxyz);
		charge_plot.plot(pdata->nz/2,0,0,HOMoments_charge);
		fields_old -> plot(pdata,pdata->nz/2,0,1,3);

		//getchar();
	}





	double push_time = 0;
	double pushcomm_time = 0;

	double nsteps_temp = 0;
	double nsteps_world = 0;
	int npiccard_total = 0;
	double npiccard_total2 = 0;
	nsubsteps_total = 0;
	nsteps_node = 0;
//	getchar();

	MPI_Barrier(MPI_COMM_WORLD);




	// Main time loop
	if(myinfo->myid_mpi == 0) printf("Entering Primary Time Loop\n");
	total_timer.start();
	for(int istep=0;istep<pdata->nsteps;istep++)
	{
		LOSolve_timer.start();
		// Guess the value of the fields for the next step
		// but only the root node needs to do this
		if(tid == 0)
		{
			if((myinfo->myid_mpi == 0)||(pdata->lo_all))
			{
				fields_next -> copy_from(fields_old);
				moments_old -> copy_from(moments[0]);
				LOsolver -> update_solution();




			}
		}
		LOSolve_timer.stop();

		// Iterate between fields and particles
		float residual = 2.0*1.0e-5;
		int picard_iter = 0;


//		particles_next[tid] -> copy_stats_from(particles_old[tid]);

		while(picard_iter < 40)
		{
			picard_iter++;
			step_timer.start();

//			printf("This is the %i'th outer Picard iteration\n",picard_iter);









			// Root node calculates the half value of the fields
			if(tid == 0)
			{
				if((myinfo->myid_mpi == 0)||(pdata->lo_all))
				{
					LOSolve_timer.start();
					Average_Fields(fields_old,fields_next,(fields_half[0]));
					LOSolve_timer.stop();
				//	fields_half[0] -> reset_plot();
				//	fields_half[0] ->  plot(pdata,pdata->nz/2,0,0,0);
				}
				else
				{
					// Set the values of the new particles to the old particles
					particles_next[tid] -> copy_from(particles_old[tid]);
				}
			}

			HOSolve_timer.start();

			MPI_Barrier(MPI_COMM_WORLD);
			// Broadcast the half value of the fields
//			printf("BroadCasting Fields\n");
			if(!pdata->lo_all)
			{
				Comm_timer.start();
				fields_half[0] -> broadcast(fields_half,myinfo);
				Comm_timer.stop();
			}


//			printf("Copying Particles\n");
			if((myinfo->myid_mpi == 0)||(pdata->lo_all))
				particles_next[tid] -> copy_from(particles_old[tid]);
//			MPI_Barrier(MPI_COMM_WORLD);
			//MPI_Barrier(MPI_COMM_WORLD);


			//MPI_Barrier(MPI_COMM_WORLD);
			push_timer.start();
			//MPI_Barrier(MPI_COMM_WORLD);
			// Move the particles
			if((myinfo->myid_mpi == 0))printf("Pushing Particles\n");
			nsteps_node += particles_next[tid] -> push(pdata,fields_half[0],moments[tid+1]);


			//MPI_Barrier(MPI_COMM_WORLD);
			push_timer.stop();
			//MPI_Barrier(MPI_COMM_WORLD);

		//	printf("Node %i clock speed: %e \%\n",pdata->mynode,push_timer.getFrequency());

			// Reduce Moment quantities
			MPI_Barrier(MPI_COMM_WORLD);
			if((myinfo->myid_mpi == 0))printf("Reducing Moments\n");
			Comm_timer.start();
			moments[tid+1] -> mpi_reduce(moments,myinfo);
			Comm_timer.stop();
			//MPI_Barrier(MPI_COMM_WORLD);
			HOSolve_timer.stop();

			if((myinfo->myid_mpi == 0))printf("Applying Weights\n");
			if(tid == 0)
				{
					if((myinfo->myid_mpi == 0)||(pdata->lo_all))
					{
						moments[0]->apply_weights();
					}
				}

			if(myinfo->myid_mpi == 0){
			//moments[0] -> reset_plot();
			//moments[0] -> plot(pdata->nz/2,0,HOMoments_currentx);
			}

			// Root Node calculates residual
			if((myinfo->tid == 0))
			{
				if((myinfo->myid_mpi == 0)||(pdata->lo_all))
				{

					residual = LOsolver->calc_residual(pdata,fields_next,fields_old,moments[0]);

				}
			}

			// Broadcast residual to all other nodes
			if(!pdata->lo_all)
				MPI_Bcast(&residual,1,MPI_FLOAT,0,MPI_COMM_WORLD);

			if((myinfo->myid_mpi == 0))printf("residual = %e\n",residual);




			// Field Solve
			LOSolve_timer.start();
			if((myinfo->myid_mpi == 0)) printf("Solving LO System\n");
			if(myinfo->tid == 0)
			{
				if((myinfo->myid_mpi == 0)||(pdata->lo_all))
				{
					fields_next -> copy_from(fields_old);

					LOsolver->solve(pdata,fields_next,moments[0]);

					if(myinfo->myid_mpi == 0){
					//fields_next -> reset_plot();
				//	fields_next -> plot(pdata,pdata->nz/2,0,0,0);
					}
				}
			}

			LOSolve_timer.stop();





			step_timer.stop();

			// If residual is below tolerance, exit the loop
			if( (residual <= 1.0e-9)||(picard_iter >= 40))
				break;



		} /* while(residual > pdata->epsilon_r) */

		npiccard_total += picard_iter;
		npiccard_total2 += (picard_iter)*(picard_iter);
		npiccard_outer = npiccard_total;

		if(myinfo->myid_mpi == 0){

			printf(" %i",istep);
			//fflush(stdout);



			if(pdata->plot_flag)
			{
				fields_next -> reset_plot();
				fields_next -> plot(pdata,0,0,1,3);
				moments[0] -> reset_plot();
				moments[0] -> plot(pdata->nz/2,0,0,HOMoments_currentxyz);

				charge_plot.reset_plot();
				charge_plot.plot(pdata->nz/2,0,0,HOMoments_charge);
			}


			//printf("Doing %e / %e particle steps per second for %e steps\n",1.0e3*nsteps_total/HOSolve_timer.get_cummulative(),
				//	1.0e3*nsteps_total/push_timer.get_cummulative(),nsteps_total);


		}

		// Swap the new particles and fields with the old
//		FieldData* fields_temp = fields_old;
//		fields_old = fields_next;
//		fields_next = fields_temp;

		ParticleList** particles_temp = particles_old;
		particles_old = particles_next;
		particles_next = particles_temp;
		if((myinfo->myid_mpi == 0)||(pdata->lo_all))
		{

			Average_Fields(fields_old,fields_next,(fields_half[0]));
			fields_old -> copy_from(fields_next);
		//	fields_half[0] -> reset_plot();
		//	fields_half[0] ->  plot(pdata,pdata->nz/2,0,0,0);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// Plot Fields and Moments
		initializer->check_step(particles_old[tid],moments,moments_old,&fields_next,myinfo);
	} /* for(int istep=0;istep<pdata->nsteps;istep++) */

	nsteps_temp = nsteps_node;
	npiccard_outer = npiccard_total;
	npiccard_outer2 = npiccard_total2;
	pdata->npiccard_outer = npiccard_total;



	MPI_Allreduce(&nsteps_temp,&nsteps_world,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	nsubsteps_total = nsteps_world;

	MPI_Barrier(MPI_COMM_WORLD);
	total_timer.stop();


	initializer->finish(particles_old[tid],moments,moments_old,&fields_next,myinfo);

//	getchar();

		save_timing();

//		getchar();
	//moments[0]->close_plot();
	//fields_next -> close_plot();








}

double2 calc_perf(CPUTimer* timer,double nsubsteps_total,int num_nodes)
{
	double2 result; // result, std

	double temp = 1.0e6*timer->get_cummulative()/num_nodes;
	double temp2 = temp*temp*num_nodes;

	double t1,t2;

	MPI_Allreduce(&temp,&t1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	MPI_Allreduce(&temp2,&t2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	result.x = t1/nsubsteps_total;
	result.y = sqrt(fabs(t1*t1 - t2))/nsubsteps_total;

	return result;

}

double2 calc_particle_perf(PlasmaData* pdata,ParticleList* plist,int itime,int idevice)
{
	double2 result; // result, std

	int num_devices;
	if(idevice == 1)
		num_devices = pdata->ngpus;
	else
		num_devices = pdata->num_nodes - pdata->ngpus;

	double temp;
	if(plist->device_type == idevice)
		temp = plist->get_cummulative_time(itime)/num_devices;
	else
		temp = 0;

	double temp2 = temp*temp;

	double t1,t2;

	MPI_Allreduce(&temp,&t1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	MPI_Allreduce(&temp2,&t2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	result.x = t1;
	result.y = sqrt(fabs(t1*t1 - t2));

	return result;

}

std::string ImplicitPIC::calc_timing(double perf_weight,const char* legend)
{
	double2 HOSolve_perf,push_perf,communication_perf,LOSolve_perf,total_perf,step_perf;

	double temp;

	double steps = 1.0e6/(perf_weight*pdata->num_nodes);

	HOSolve_perf = calc_perf(&HOSolve_timer,steps,pdata->num_nodes);
	push_perf = calc_perf(&push_timer,steps,pdata->num_nodes);
	communication_perf = calc_perf(&Comm_timer,steps,pdata->num_nodes);
	step_perf = calc_perf(&step_timer,steps,pdata->num_nodes);

	//communication_perf = HOSolve_perf - push_perf;

	if(pdata->lo_all)
		LOSolve_perf = calc_perf(&LOSolve_timer,steps,pdata->num_nodes);
	else{
		LOSolve_perf.x = perf_weight*LOSolve_timer.get_cummulative()*pdata->num_nodes;
		LOSolve_perf.y = 0;
	}
	total_perf.x = perf_weight*total_timer.get_cummulative()*pdata->num_nodes;
	total_perf.y = 0;

	char tempc[128];
	std::string result("");

	sprintf(tempc,"Performance in %s\n",legend);
	result += tempc;
	sprintf(tempc,"HOSolve step took: %e +- %e\n",HOSolve_perf.x,HOSolve_perf.y);
	result += tempc;
	sprintf(tempc,"push step took: %e +- %e\n",push_perf.x,push_perf.y);
	result += tempc;
	sprintf(tempc,"Communication took: %e +- %e\n",communication_perf,communication_perf.y);
	result += tempc;
	sprintf(tempc,"LOSolve step took: %e +- %e\n",LOSolve_perf.x,LOSolve_perf.y);
	result += tempc;
	sprintf(tempc,"Total step took: %e +- %e\n",total_perf.x,total_perf.y);
	result += tempc;


	return result;
}


void ImplicitPIC::save_timing()
{

	particles_old[0] -> subcycle_stats(pdata);
	particles_old[0] -> piccard_stats(pdata);

	double2 HOSolve_perf,push_perf,communication_perf,LOSolve_perf,total_perf,step_perf;

	double perf_weight = 1.0e6/(nsubsteps_total*pdata->num_nodes);

	double temp;

//	temp = HOSolve_timer.get_cummulative();
//
//	MPI_Allreduce(&temp,&HOSolve_perf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//
//	HOSolve_perf = perf_weight*HOSolve_perf;
//
//	temp = push_timer.get_cummulative();
//
//	MPI_Allreduce(&temp,&push_perf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//
//	push_perf = perf_weight*push_perf;
//
//	temp = Comm_timer.get_cummulative();
//
//	MPI_Allreduce(&temp,&communication_perf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//
//	communication_perf = perf_weight*communication_perf;
//
//	temp = step_timer.get_cummulative();
//
//	MPI_Allreduce(&temp,&step_perf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//
//	step_perf = perf_weight*step_perf;

	HOSolve_perf = calc_perf(&HOSolve_timer,nsubsteps_total,pdata->num_nodes);
	push_perf = calc_perf(&push_timer,nsubsteps_total,pdata->num_nodes);
	communication_perf = calc_perf(&Comm_timer,nsubsteps_total,pdata->num_nodes);
	step_perf = calc_perf(&step_timer,nsubsteps_total,pdata->num_nodes);



	//communication_perf = HOSolve_perf - push_perf;

	if(pdata->lo_all)
		LOSolve_perf = calc_perf(&LOSolve_timer,nsubsteps_total,pdata->num_nodes);
	else{
		LOSolve_perf.x = perf_weight*LOSolve_timer.get_cummulative()*pdata->num_nodes;
		LOSolve_perf.y = 0;
	}
	total_perf.x = perf_weight*total_timer.get_cummulative()*pdata->num_nodes;
	total_perf.y = 0;


	if(pdata->mynode == 0){
	printf("Sim Params:\n");
	printf("Num Cores: %i\n",pdata->num_cores);
	printf("Num Nodes: %i\n",pdata->num_nodes);
	printf("CPU Vec Length: %i\n",pdata->cpu_vec_length);
	printf("Num Ptcls per Node: %i\n",pdata->nptcls);
	printf("Nx : %i\n",pdata->nx);
	printf("Total Subcycle Steps: %e\n", (double)nsubsteps_total);
	printf("Total number of outer piccard: %i\n",pdata->npiccard_outer);

	printf("Performance in ns/particle-substep\n");
	printf("HOSolve step took: %f +- %f\n",HOSolve_perf.x,HOSolve_perf.y);
	printf("push step took: %f +- %f\n",push_perf.x,push_perf.y);
	printf("Communication took: %f +- %f\n",communication_perf,communication_perf.y);
	printf("LOSolve step took: %f +- %f\n",LOSolve_perf.x,LOSolve_perf.y);
	printf("Total step took: %f +- %f\n",total_perf.x,total_perf.y);

	}

	std::string piccard_perf = calc_timing(1.0e6/(npiccard_outer*pdata->num_nodes),"ns/outer piccard");

	if(pdata->mynode == 0)
		printf("%s\n",piccard_perf.c_str());

	MPI_Barrier(MPI_COMM_WORLD);
	printf("Node %i took %e ms to push %e steps\n",pdata->mynode,push_timer.get_cummulative(),(double)nsteps_node);


	int num_gpus;
	int igpu = 0;

	if(pdata->device_type == 1)
		igpu = 1;

	MPI_Reduce(&igpu,&num_gpus,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

	double np = npiccard_outer/((double)pdata->nsteps);
	double np2 = npiccard_outer2/((double)pdata->nsteps);
	double npiccard_std = sqrt(fabs(np*np - np2));
//	MPI_Barrier(MPI_COMM_WORLD);
//	printf("Node %i took %e ms to HO\n",pdata->mynode,HOSolve_timer.get_cummulative());
//




//	char filename[60];
//	FILE* fp;
//
//	int nparams = 19;
//	int ntimes;
//	char* run_params[nparams];
//	char* time_names[ntimes];
//
//
//
//	run_params[0] = "dt";
//	run_params[1] = "nsteps";
//	run_params[2] = "nx";
//	run_params[3] = "ny";
//	run_params[4] = "nz";
//	run_params[5] = "n_electrons";
//	run_params[6] = "n_ions";
//	run_params[7] = "nptcls_gpu";
//	run_params[8] = "nptcls_cpu";
//	run_params[9] = "ndimensions";
//
//	run_params[10] = "Lx";
//	run_params[11] = "Ly";
//	run_params[12] = "Lz";
//	run_params[13] = "Problem_Type";
//	run_params[14] = "mass_ratio";
//	run_params[15] = "cpu_vec_length";
//	run_params[16] = "num_cores";
//	run_params[17] = "num_nodes";
//
//	run_params[18] = "Num_SubCycle_Steps";
//
//
//
//	time_names[0] = "HOSolve_time";
//	time_names[1] = "Push_time";
//	time_names[2] = "Communication_time";
//	time_names[3] = "LOSolve_time";
//	time_names[4] = "Step_time";
//	time_names[5] = "Total_time";

	int nparticle_times = 10;
	double2 particle_times[nparticle_times];
	for(int i=0;i<nparticle_times;i++)
	{
		particle_times[i] = calc_particle_perf(pdata,particles_old[0],i,0);
	}

	if(pdata->mynode == 0)
	{
		char filename[60];
		FILE* fp;

		int nparams = 13;
		int ntimes = 16;
		char* run_params[nparams];
		char* time_names[ntimes];

		double2 run_times[ntimes];
		double2 adjusted_times[ntimes];



		run_params[0] = "dt";
		run_params[1] = "nsteps";
		run_params[2] = "nsubsteps_total";
		run_params[3] = "nptcls_total";
		run_params[4] = "ncells";
		run_params[5] = "npiccard_outer";
		run_params[6] = "npiccard_std";
		run_params[7] = "vector_length";
		run_params[8] = "num_cores";

		run_params[9] = "num_nodes";
		run_params[10] = "num_gpus";
		run_params[11] = "lo-all";
		run_params[12] = "OutputID";



		time_names[0] = "HOSolve_time";
		time_names[1] = "Push_time";
		time_names[2] = "Comm_time";
		time_names[3] = "LOSolve_time";
		time_names[4] = "Step_time";
		time_names[5] = "Total_time";
		time_names[6] = "Ppicard_time";
		time_names[7] = "Accel_time";
		time_names[8] = "Tally_time";
		time_names[9] = "Crossing_time";
		time_names[10] = "Dtau_est_time";
		time_names[11] = "ChargeS2Tally_time";
		time_names[12] = "PLoadStore_time";
		time_names[13] = "PPiccardOther_time";
		time_names[14] = "Push2_time";
		time_names[15] = "PushOther_time";



		adjusted_times[0] = HOSolve_perf;
		adjusted_times[1] = push_perf;
		adjusted_times[2] = communication_perf;
		adjusted_times[3] = LOSolve_perf;
		adjusted_times[4] = step_perf;
		adjusted_times[5] = total_perf;





		for(int i=0;i<ntimes;i++)
		{
			run_times[i].x = adjusted_times[i].x / (perf_weight*pdata->num_nodes) ;
			run_times[i].y = adjusted_times[i].y / (perf_weight*pdata->num_nodes) ;

		}

		for(int i=0;i<nparticle_times;i++)
		{
			run_times[i+6] = particle_times[i];
			adjusted_times[i+6].x = run_times[i+6].x*1.0e6/nsubsteps_total;
			adjusted_times[i+6].y = run_times[i+6].y*1.0e6/nsubsteps_total;
		}

		// Check the output directory
		mkpath("./benchmarks",0777);

		// Setup the filename
		sprintf(filename,"./benchmarks/benchmark%i.dat",pdata->runid);

		// Check to see if the file exists
		fp = fopen(filename,"r+");

		// If the file doesn't exist yet, create it and write the top line
		if(fp == NULL)
		{

			fp = fopen(filename,"w");
			char header[nparams*16+35*(ntimes)];
			for(int i=0;i<nparams;i++)
			{
				fprintf(fp,"%s,",run_params[i]);
			}

			for(int i=0;i<ntimes;i++)
			{
				fprintf(fp,"%s(ms), %s(ms), %s(ns), %s(ns),",time_names[i],"std","adjusted","std");
			}

			fprintf(fp,"\n");
		}

		fclose(fp);

		fp = fopen(filename,"a");

		char lineout[nparams*16+35*(ntimes)];

		std::string nameout("");
		nameout += pdata->SimName;
		nameout += "/";
		nameout += pdata->output_id;

		fprintf(fp,"%f,%i,%e,"
				   "%i,%i,"
				   "%i,%e,"
				   "%i,%i,%i,"
				   "%i,%i,",
					pdata->dt,pdata->nsteps,(double)nsubsteps_total,
					pdata->nptcls_species[0]+pdata->nptcls_species[1],pdata->nx,
					npiccard_outer,npiccard_std,
					pdata->cpu_vec_length,pdata->num_cores,pdata->num_nodes,
					num_gpus,pdata->lo_all);

		fprintf(fp,"%s,",nameout.c_str());


		for(int i=0;i<ntimes;i++)
		{
			fprintf(fp,"%f,%f,%f,%f,",run_times[i].x,run_times[i].y,adjusted_times[i].x,adjusted_times[i].y);
		}

		fprintf(fp,"\n");


		fclose(fp);

		printf("\n");
		for(int i=0;i<ntimes;i++)
		{
			char temp[30];
			char temp2[20];

			sprintf(temp,"%s runtime: ",time_names[i]);
			sprintf(temp2,"%*f",(40-strlen(temp)),adjusted_times[i].x);
			printf("%s %s(ns)\n",temp,temp2);
		}

	}

	output();
}

OutputPage ImplicitPIC::get_subcycle_dists()
{
	OutputPage subcycle_dists;
	/* Subcycle output data */
	int ndists = pdata->nspecies * pdata->ndevices;
	double dists[ndists][NSUBCYCLE_BINS+2];
	double* dist_temp;
	double null_dist[NSUBCYCLE_BINS+2];
	int nnodes[ndists];

	double4 my_sub_stats = particles_old[0]->subcycle_stats(pdata);
	double4 substats[pdata->num_nodes];
	double tempstats[4] = {my_sub_stats.x,my_sub_stats.y,my_sub_stats.z,my_sub_stats.w};
	double substats_t[4*pdata->num_nodes];
	printf("substats: %f, %f, %f, %f\n",tempstats[0],tempstats[1],tempstats[2],tempstats[3]);


	MPI_Barrier(MPI_COMM_WORLD);

	// Gather all of the subcycle statistics
	MPI_Gather(tempstats,4,MPI_DOUBLE,substats_t,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(pdata->mynode == 0)
	for(int i=0;i<pdata->num_nodes;i++)
	{
		substats[i] = make_double4(substats_t[4*i],substats_t[4*i+1],substats_t[4*i+2],substats_t[4*i+3]);
		printf("substats[%i]: %f, %f, %f, %f\n",i,substats[i].x,substats[i].y,substats[i].z,substats[i].w);

	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0;i<NSUBCYCLE_BINS+2;i++)
		null_dist[i] = 0;

	// get the distributions
	dist_temp = subcycle_dist(pdata,particles_old[0]);


	// get reduce everything
	printf("Reducing Subcycle Distributions\n");
	for(int i = 0;i<pdata->nspecies;i++)
		for(int j=0;j<pdata->ndevices;j++)
		{
			int node_count = 0;
			double* dist_reduce = null_dist;
			if(pdata->my_species == i && pdata->device_type == j)
			{
				dist_reduce = dist_temp;
				node_count = 1;
			}

			MPI_Reduce(dist_reduce,dists[j + pdata->ndevices*i],NSUBCYCLE_BINS+2,
					MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(&node_count,nnodes+j + pdata->ndevices*i,1,
					MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		}

	if(pdata->mynode == 0)
	{
		printf("Normalizing Subcycle Distributions\n");
		// Normalize everything
		for(int i = 0;i<ndists;i++)
			for(int k=0;k<NSUBCYCLE_BINS+2;k++)
				dists[i][k] /= (double)nnodes[i];

		for(int i=0;i<pdata->nspecies;i++)
		{
			double mins = dists[pdata->ndevices*i][NSUBCYCLE_BINS];
			double maxs = dists[pdata->ndevices*i][NSUBCYCLE_BINS+1];
			for(int j=0;j<pdata->ndevices;j++)
			{
				mins = fmin(mins,dists[j+pdata->ndevices*i][NSUBCYCLE_BINS]);
				maxs = fmax(maxs,dists[j+pdata->ndevices*i][NSUBCYCLE_BINS+1]);
			}

			for(int j=0;j<pdata->ndevices;j++)
			{
				dists[j+pdata->ndevices*i][NSUBCYCLE_BINS] = mins;
				dists[j+pdata->ndevices*i][NSUBCYCLE_BINS+1] = maxs;
			}
		}

		if(pdata->plot_flag){
		// Plot distributions to gnuplot.
		gnuplot_ctrl* plots[ndists];
		for(int i=0;i<ndists;i++)
		{
			float tempys[NSUBCYCLE_BINS];
			float tempxs[NSUBCYCLE_BINS];
			double mins = dists[i][NSUBCYCLE_BINS];
			double maxs = dists[i][NSUBCYCLE_BINS+1];
			double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);

			if(!isnan(dists[i][0])){
			for(int k=0;k<NSUBCYCLE_BINS;k++)
			{
				tempys[k] = dists[i][k];
				tempxs[k] = k*dsdi+mins;
			}

			plots[i] = gnuplot_init();

			gnuplot_plot_xy(plots[i],tempxs,tempys,NSUBCYCLE_BINS);
			}


		}
		}

		char temp[128];

		// Do the Node stats
		printf("Printing Node Subcycle Stats\n");
		subcycle_dists.nextline() = "#Node, Average Subs, std, min, max";
		for(int i=0;i<pdata->num_nodes;i++)
		{
			sprintf(temp,"#%i, %f, %f, %f, %f",
					i,substats[i].x,substats[i].y,substats[i].z,substats[i].w);
			subcycle_dists.nextline() = temp;
		}

		// Do the distributions
		printf("Printing Subcycle Distributions\n");
		for(int i=0;i<pdata->nspecies;i++)
		{
			std::string legend("#NSubcycles");
			for(int j=0;j<pdata->ndevices;j++)
			{
				sprintf(temp,",device:%i",j);
				legend += temp;
			}
			subcycle_dists.nextline() = legend;
			sprintf(temp,"#ispecies: %i",i);
			subcycle_dists.nextline() = temp;

			// Do the arrays
			double mins = dists[pdata->ndevices*i][NSUBCYCLE_BINS];
			double maxs = dists[pdata->ndevices*i][NSUBCYCLE_BINS+1];
			double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);


			for(int k=0;k<NSUBCYCLE_BINS;k++)
			{
				std::string &line = subcycle_dists.nextline();
				sprintf(temp,"%f",k*dsdi + mins);
				line = temp;

				for(int j=0;j<pdata->ndevices;j++)
				{
					double count = dists[j + pdata->ndevices*i][k];

					sprintf(temp,", %e",count);
		//			printf("line[%i] = %e\n",k,count);
					line += temp;
				}



			}

		}
	}

	return subcycle_dists;
}

OutputPage ImplicitPIC::get_piccard_dists()
{
	OutputPage piccard_dists;
	/* Subcycle output data */
	int ndists = pdata->nspecies * pdata->ndevices;
	double dists[ndists][NSUBCYCLE_BINS+2];
	double* dist_temp;
	double null_dist[NSUBCYCLE_BINS+2];
	int nnodes[ndists];

	double4 my_sub_stats = particles_old[0]->piccard_stats(pdata);
	double4 substats[pdata->num_nodes];
	double tempstats[4] = {my_sub_stats.x,my_sub_stats.y,my_sub_stats.z,my_sub_stats.w};
	double substats_t[4*pdata->num_nodes];
	printf("piccardstats: %f, %f, %f, %f\n",tempstats[0],tempstats[1],tempstats[2],tempstats[3]);


	MPI_Barrier(MPI_COMM_WORLD);

	// Gather all of the subcycle statistics
	MPI_Gather(tempstats,4,MPI_DOUBLE,substats_t,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(pdata->mynode == 0)
	for(int i=0;i<pdata->num_nodes;i++)
	{
		substats[i] = make_double4(substats_t[4*i],substats_t[4*i+1],substats_t[4*i+2],substats_t[4*i+3]);
		printf("substats[%i]: %f, %f, %f, %f\n",i,substats[i].x,substats[i].y,substats[i].z,substats[i].w);

	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0;i<NSUBCYCLE_BINS+2;i++)
		null_dist[i] = 0;

	// get the distributions
	dist_temp = piccard_dist(pdata,particles_old[0]);


	// get reduce everything
	if(pdata->mynode == 0)printf("Reducing Piccard Distributions\n");
	for(int i = 0;i<pdata->nspecies;i++)
		for(int j=0;j<pdata->ndevices;j++)
		{
			int node_count = 0;
			double* dist_reduce = null_dist;
			if(pdata->my_species == i && pdata->device_type == j)
			{
				dist_reduce = dist_temp;
				node_count = 1;
			}

			MPI_Reduce(dist_reduce,dists[j + pdata->ndevices*i],NSUBCYCLE_BINS+2,
					MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(&node_count,nnodes+j + pdata->ndevices*i,1,
					MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		}

	if(pdata->mynode == 0)
	{
		printf("Normalizing Piccard Distributions\n");
		// Normalize everything
		for(int i = 0;i<ndists;i++)
			for(int k=0;k<NSUBCYCLE_BINS+2;k++)
				dists[i][k] /= (double)nnodes[i];

		for(int i=0;i<pdata->nspecies;i++)
		{
			double mins = dists[pdata->ndevices*i][NSUBCYCLE_BINS];
			double maxs = dists[pdata->ndevices*i][NSUBCYCLE_BINS+1];
			for(int j=0;j<pdata->ndevices;j++)
			{
				mins = fmin(mins,dists[j+pdata->ndevices*i][NSUBCYCLE_BINS]);
				maxs = fmax(maxs,dists[j+pdata->ndevices*i][NSUBCYCLE_BINS+1]);
			}

			for(int j=0;j<pdata->ndevices;j++)
			{
				dists[j+pdata->ndevices*i][NSUBCYCLE_BINS] = mins;
				dists[j+pdata->ndevices*i][NSUBCYCLE_BINS+1] = maxs;
			}
		}

		// Plot distributions to gnuplot.
		if(pdata->plot_flag){
		gnuplot_ctrl* plots[ndists];
		for(int i=0;i<ndists;i++)
		{
			float tempys[NSUBCYCLE_BINS];
			float tempxs[NSUBCYCLE_BINS];
			double mins = dists[i][NSUBCYCLE_BINS];
			double maxs = dists[i][NSUBCYCLE_BINS+1];
			double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);
			if(!isnan(dists[i][0])){
			for(int k=0;k<NSUBCYCLE_BINS;k++)
			{
				tempys[k] = dists[i][k];
				tempxs[k] = k*dsdi+mins;
			}

			plots[i] = gnuplot_init();

			gnuplot_plot_xy(plots[i],tempxs,tempys,NSUBCYCLE_BINS);

			}

		}
		}

		char temp[128];

		// Do the Node stats
		printf("Printing Node Piccard Stats\n");
		piccard_dists.nextline() = "#Node, Average Iters, std, min, max";
		for(int i=0;i<pdata->num_nodes;i++)
		{
			sprintf(temp,"#%i, %f, %f, %f, %f",
					i,substats[i].x,substats[i].y,substats[i].z,substats[i].w);
			piccard_dists.nextline() = temp;
		}

		// Do the distributions
		printf("Printing Subcycle Distributions\n");
		for(int i=0;i<pdata->nspecies;i++)
		{
			std::string legend("#NPiccard");
			for(int j=0;j<pdata->ndevices;j++)
			{
				sprintf(temp,",device:%i",j);
				legend += temp;
			}
			piccard_dists.nextline() = legend;
			sprintf(temp,"#ispecies: %i",i);
			piccard_dists.nextline() = temp;

			// Do the arrays
			double mins = dists[pdata->ndevices*i][NSUBCYCLE_BINS];
			double maxs = dists[pdata->ndevices*i][NSUBCYCLE_BINS+1];
			double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);


			for(int k=0;k<NSUBCYCLE_BINS;k++)
			{
				std::string &line = piccard_dists.nextline();
				sprintf(temp,"%f",k*dsdi + mins);
				line = temp;
				for(int j=0;j<pdata->ndevices;j++)
				{
					double count = dists[j + pdata->ndevices*i][k];

					sprintf(temp,", %e",count);
		//			printf("line[%i] = %e\n",k,count);
					line += temp;
				}



			}

		}
	}

	return piccard_dists;
}




void ImplicitPIC::output()
{
	OutputPage subcycle_page;
	OutputPage piccard_page;

	subcycle_page = get_subcycle_dists();
	piccard_page = get_piccard_dists();


	// Get subcycle and piccard distributions
	if(pdata->mynode == 0)
	{
		// First setup the output path name
		printf("Setting Output Path Name\n");
		char output_path[128];
		sprintf(output_path,"./output/%s/%s",pdata->SimName,pdata->output_id);

		// Make the output path if it doesn't exist
		printf("Checking Output Path\n");
		mkpath(output_path,0777);

		char filename[128];
		sprintf(filename,"%s/SubcycleDist.dat",output_path);

		printf("Writing Subcycle Distributions\n");
		subcycle_page.writepage(filename);
		sprintf(filename,"%s/PiccardDist.dat",output_path);
		printf("Writing Subcycle Distributions\n");
		piccard_page.writepage(filename);

	}


}




double* subcycle_dist(PlasmaData* pdata,ParticleList* plist)
{
	// allocate distribution array
	double* dist = (double*)malloc((NSUBCYCLE_BINS+2)*sizeof(double));
	int* num_subcycles_temp = (int*)malloc(plist->nptcls*sizeof(int));

	if(plist->device_type == 1){
#ifndef NO_CUDA
	CUDA_SAFE_CALL(cudaMemcpy(num_subcycles_temp,plist->num_subcycles,plist->nptcls*sizeof(int),
			cudaMemcpyDeviceToHost));
#endif
	}else
		memcpy(num_subcycles_temp,plist->num_subcycles,plist->nptcls*sizeof(int));

	double4 stats = plist->subcycle_stats(pdata);
	double mins_l = stats.z; // local min
	double maxs_l = stats.w; // local max
	double limits[2];
	double mins,maxs;
	// Get the global min and max for the species / device combo
	for(int i = 0;i<pdata->nspecies;i++)
	{
			double limits_t[2] = {500,0};
			if(plist->ispecies == i)
			{
				limits_t[0] = fmax(mins_l,0);
				limits_t[1] = maxs_l;
			}

			MPI_Allreduce(limits_t,limits,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
			MPI_Allreduce(limits_t+1,limits+1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

			if(plist->ispecies == i)
			{
				mins = limits[0];
				maxs = limits[1];
			}
	}



	dist[NSUBCYCLE_BINS] = mins; // Minimum number of subcycles
	dist[NSUBCYCLE_BINS+1] = maxs; // Maximum number of subcycles

	for(int i=0;i<NSUBCYCLE_BINS;i++)
		dist[i] = 0;

	double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);
	double dids = 1.0/dsdi;

	double scale = 1.0/((double)plist->nptcls);
	double scale2 = 1.0/pdata->npiccard_outer;

	for(int i=0;i<plist->nptcls;i++)
	{
		double delta = ((num_subcycles_temp[i]*scale2-mins)*dids);
		int isub = floor(delta);
		delta = delta - isub;
		for(int j=0;j<2;j++)
		{
			int i_out;
			double xp;
			double temp;

			i_out = isub + j;
			xp = j - delta;

			dist[std::min(std::max(i_out,0),NSUBCYCLE_BINS-1)] += S1_shape(xp)*scale;

		}

	}

	free(num_subcycles_temp);

	return dist;

}



double* piccard_dist(PlasmaData* pdata,ParticleList* plist)
{
	// allocate distribution array
	double* dist = (double*)malloc((NSUBCYCLE_BINS+2)*sizeof(double));
	int* num_subcycles_temp = (int*)malloc(plist->nptcls*sizeof(int));
	double* num_piccard_temp = (double*)malloc(plist->nptcls*sizeof(double));

	if(plist->device_type == 1){
#ifndef NO_CUDA
	CUDA_SAFE_CALL(cudaMemcpy(num_subcycles_temp,plist->num_subcycles,plist->nptcls*sizeof(int),
			cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(num_piccard_temp,plist->num_piccard,plist->nptcls*sizeof(double),
			cudaMemcpyDeviceToHost));
#endif

	}else{
		memcpy(num_subcycles_temp,plist->num_subcycles,plist->nptcls*sizeof(int));
		memcpy(num_piccard_temp,plist->num_piccard,plist->nptcls*sizeof(double));
	}

	double4 stats = plist->piccard_stats(pdata);
	double mins_l = stats.z; // local min
	double maxs_l = stats.w; // local max
	double limits[2];
	double mins,maxs;
	// Get the global min and max for the species / device combo
	for(int i = 0;i<pdata->nspecies;i++)
	{
			double limits_t[2] = {500,0};
			if(plist->ispecies == i)
			{
				limits_t[0] = fmax(mins_l,0);
				limits_t[1] = maxs_l;
			}

			MPI_Allreduce(limits_t,limits,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
			MPI_Allreduce(limits_t+1,limits+1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

			if(plist->ispecies == i)
			{
				mins = limits[0];
				maxs = limits[1];
			}
	}



	dist[NSUBCYCLE_BINS] = mins; // Minimum number of subcycles
	dist[NSUBCYCLE_BINS+1] = maxs; // Maximum number of subcycles

	for(int i=0;i<NSUBCYCLE_BINS;i++)
		dist[i] = 0;

	double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);
	double dids = 1.0/dsdi;

	double scale = 1.0/((double)plist->nptcls);
	double scale2 = 1.0/pdata->npiccard_outer;

	for(int i=0;i<plist->nptcls;i++)
	{
		double delta = ((num_piccard_temp[i]/num_subcycles_temp[i]-mins)*dids);
		int isub = floor(delta);
		delta = delta - isub;
		for(int j=0;j<2;j++)
		{
			int i_out;
			double xp;
			double temp;

			i_out = isub + j;
			xp = j - delta;

			dist[std::min(std::max(i_out,0),NSUBCYCLE_BINS-1)] += S1_shape(xp)*scale;

		}

	}

	free(num_subcycles_temp);
	free(num_piccard_temp);

	return dist;

}




































