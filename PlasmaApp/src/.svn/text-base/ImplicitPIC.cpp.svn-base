/*-------------------------------------------------------------------------*/
/**
	@file		ImplicitPIC.cpp
	@author	J. Payne
	@date		12/21/2012


*/
/*--------------------------------------------------------------------------*/
#include "ImplicitPIC.h"
#include "NodeHOMoments.h"
#include "NodeFieldData.h"
#include "NodeParticleList.h"
#include "PlasmaData.h"
#include "HOMomentsCPU.h"
#include "FieldDataCPU.h"
#include "ProblemInitializer.h"
#include "FieldSolver.h"
#include "ParallelInfo.h"
#include <gnuplot_i.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "Util/mkpath.h"
#include "Util/OutputPage.h"
#include "RunData.h"
#include "OutputCtrl.h"


double& ImplicitPIC::nsubsteps_total(){return pdata->rstatus->nsubsteps_total;}
/// Total number of subcycle steps for individual node
double& ImplicitPIC::nsteps_node(){return pdata->rstatus->nsteps_node;}
/// Total number of outer piccard iterations
int& ImplicitPIC::npiccard_outer(){return pdata->rstatus->npiccard_outer;}
/// Squared number of piccard iterations
double& ImplicitPIC::npiccard_outer2(){return pdata->rstatus->npiccard_outer2;}

int ImplicitPIC::myid_mpi(){return pdata->node_info->rank_g;}

ImplicitPIC::ImplicitPIC(PlasmaData* pdata_in, FieldSolver* LOsolver_in,
		ProblemInitializer* initializer_in)
{



	pdata = pdata_in;

	LOsolver = LOsolver_in;
	initializer = initializer_in;



	rstatus = pdata->rstatus;

	if(myid_mpi() == 0)
		debug_log.print("Setting up output control\n");
	output_ctrl = pdata->output_ctrl;




	fields_old = new FieldDataCPU();
	fields_next = new FieldDataCPU();
	fields_half = new NodeFieldData();

	particles_old = new NodeParticleList();
	particles_next = new NodeParticleList();

	moments = new NodeHOMoments();

	if(myid_mpi() == 0)
		debug_log.print("Allocating Particles\n");
	particles_old -> allocate(pdata);
	particles_next -> allocate(pdata,particles_old);

	if(myid_mpi() == 0)
		debug_log.print("Allocating Fields\n");
	fields_old -> allocate(pdata);
	fields_next -> allocate(pdata);
	fields_half -> allocate(pdata);

	if(myid_mpi() == 0)
		debug_log.print("Allocating Moments\n");
	moments -> allocate(pdata);





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
	// Initialize fields and particles
	if(myid_mpi() == 0) debug_log.print("Initializing Particles\n");
	initializer -> initialize_particles(particles_old,moments);
	particles_next -> copy_from(particles_old);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid_mpi() == 0)debug_log.print("Finished initializing particles\n");


	MPI_Barrier(MPI_COMM_WORLD);
	// Reduce Moments
	if(myid_mpi() == 0) debug_log.print("Reducing Moments\n");

	moments -> reduce();

	if(myid_mpi() == 0) debug_log.print("Applying Weights\n");

	moments -> apply_weights();

	moments->UpdateMoments();



	// Initialize Fields
	if(myid_mpi() == 0) printf("Initializing Fields\n");

	initializer -> initialize_fields(fields_old,moments);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid_mpi() == 0) printf("Broadcasting Fields\n");

	MPI_Barrier(MPI_COMM_WORLD);


	fields_next -> copy_from(fields_old);


	if(myid_mpi() == 0) debug_log.print("Copying Fields\n");

	if(myid_mpi() == 0) debug_log.print("Init LO\n");
	LOsolver->init(pdata,fields_old,fields_half,fields_next,moments);


	fields_next -> copy_from(fields_old);

	fields_half -> average(fields_old,fields_next);


	if((myid_mpi() == 0)||(pdata->rdata->lo_all))
	{
		LOsolver->calc_residual(pdata,fields_next,fields_old,moments);
	}

	// Plotting Stuff
	debug_log.print("Plotting Stuff\n");
	output_ctrl->InitialOutput(fields_half,fields_old,fields_next,moments,particles_old);

	debug_log.print("Finished Plotting Stuff\n");


	printf("Finished Coprocessing Initialization\n");




	double push_time = 0;
	double pushcomm_time = 0;

	double nsteps_temp = 0;
	double nsteps_world = 0;
	int npiccard_total = 0;
	double npiccard_total2 = 0;
	nsubsteps_total() = 0;
	nsteps_node() = 0;
//	getchar();

	MPI_Barrier(MPI_COMM_WORLD);




	// Main time loop
	if(myid_mpi() == 0) debug_log.print("Entering Primary Time Loop\n");
	output_ctrl->total_timer->start();
	for(int istep=0;istep<pdata->nsteps;istep++)
	{
		pdata->rstatus->istep = istep+1;
		output_ctrl->LOSolve_timer->start();
		// Guess the value of the fields for the next step
		moments->UpdateMoments();
		fields_old -> copy_from(fields_next);

//		if((myid_mpi() == 0)||(pdata->rdata->lo_all))
//		{
//			if((myid_mpi() == 0))debug_log.print("Field Solve\n");
//			LOsolver -> InitStepSolve(pdata,fields_next,moments);
//		}


		output_ctrl->LOSolve_timer->stop();

		// Iterate between fields and particles
		float residual = 2.0*1.0e-5;
		int picard_iter = 0;

		pdata->rstatus->nsubsteps_step  = 0;
		if(myid_mpi() == 0) debug_log.print("Entering Outer Picard Loop\n");
		while(picard_iter < pdata->nOuterPicard_max)
		{
			Outer_log.print("Timestep %i, Outer Picard Iteration %i\n",istep,picard_iter);
			picard_iter++;
			pdata->rstatus->ipiccard_outer = picard_iter;
			output_ctrl->step_timer->start();


			// Root node calculates the half value of the fields
			if((myid_mpi() == 0)||(pdata->rdata->lo_all))
			{

				if((myid_mpi() == 0))debug_log.print("Average Fields\n");
				output_ctrl->LOSolve_timer->start();
				fields_half -> average(fields_old,fields_next);
				output_ctrl->LOSolve_timer->stop();

			}
			else
			{
				// Set the values of the new particles to the old particles
				particles_next -> copy_from(particles_old);
			}


			output_ctrl -> HOSolve_timer ->start();


			// Broadcast the half value of the fields
//			printf("BroadCasting Fields\n");
			output_ctrl -> Comm_timer->start();
			fields_half -> broadcast();
			output_ctrl -> Comm_timer->stop();


//			printf("Copying Particles\n");
			if((myid_mpi() == 0)||(pdata->rdata->lo_all))
				particles_next -> copy_from(particles_old);




			output_ctrl->push_timer->start();
			// Move the particles
			moments->set_vals(0);
			if((myid_mpi() == 0))debug_log.print("Pushing Particles\n");
			pdata->rstatus->nsubsteps_op = particles_next -> push(fields_half,moments);

			nsteps_node() += pdata->rstatus->nsubsteps_op;
			pdata->rstatus->nsubsteps_step += pdata->rstatus->nsubsteps_op;
			output_ctrl->push_timer->stop();

			// Reduce Moment quantities
			output_ctrl->HOSolve_timer->stop();

			if((myid_mpi() == 0))debug_log.print("Reducing Moments\n");
			output_ctrl->Comm_timer->start();
			moments -> reduce();
			output_ctrl->Comm_timer->stop();
			//MPI_Barrier(MPI_COMM_WORLD);
			output_ctrl->LOSolve_timer->start();
			if((myid_mpi() == 0))debug_log.print("Applying Weights\n");
			if((myid_mpi() == 0)||(pdata->rdata->lo_all))
			{
				moments->apply_weights();
			}


			// Root Node calculates residual
			if((myid_mpi() == 0)||(pdata->rdata->lo_all))
			{
				residual = LOsolver->calc_residual(pdata,fields_next,fields_old,moments);
			}


			if((myid_mpi() == 0)||(pdata->rdata->lo_all))
			{
				if((myid_mpi() == 0))debug_log.print("Field Solve\n");
				LOsolver -> solve(pdata,fields_next,moments);
			}
			output_ctrl->LOSolve_timer->stop();

			// Broadcast residual to all other nodes
			if(!pdata->rdata->lo_all)
				MPI_Bcast(&residual,1,MPI_FLOAT,0,MPI_COMM_WORLD);

			if((myid_mpi() == 0))debug_log.print("residual = %e\n",residual);





			output_ctrl->step_timer->stop();

			// If residual is below tolerance, exit the loop
			if( (residual <= pdata->tol_outer_picard)||(picard_iter >= pdata->nOuterPicard_max))
				break;



		} /* while(residual > pdata->epsilon_r) */

//		if((myid_mpi() == 0)||(pdata->rdata->lo_all))
//		{
//			output_ctrl->LOSolve_timer->start();
////			if((myid_mpi() == 0))debug_log.print("Field Copy\n");
////			fields_next -> copy_from(fields_old);
////			if((myid_mpi() == 0))debug_log.print("Field Solve\n");
////			LOsolver -> solve(pdata,fields_next,moments);
////			if((myid_mpi() == 0))debug_log.print("Average Fields\n");
//			fields_half -> average(fields_old,fields_next);
//			output_ctrl->LOSolve_timer->stop();
//
//		}
//		fields_half -> broadcast();

		npiccard_total += picard_iter;
		npiccard_total2 += (picard_iter)*(picard_iter);
		npiccard_outer() = npiccard_total;
		pdata->rstatus->npiccard_outer_last = picard_iter;
		pdata->rstatus->npiccard_outer_avg = npiccard_outer()/(istep+1.0);

		if(myid_mpi() == 0)
			printf(" %i",istep);

		// Plot Fields and Moments
		if((myid_mpi() == 0))debug_log.print("Step Outputs\n");

		output_ctrl->StepOutput(moments,particles_next,fields_old,fields_next);
		initializer->check_step(particles_next,moments,fields_old,fields_next);
		// Coprocessing Fields and Moments


		// Swap the new particles with the old
//		NodeParticleList* particles_temp = particles_old;
//		particles_old = particles_next;
//		particles_next = particles_temp;

		particles_old->copy_from(particles_next);


		npiccard_outer() = npiccard_total;
		npiccard_outer2() = npiccard_total2;

	} /* for(int istep=0;istep<pdata->nsteps;istep++) */





	nsteps_temp = nsteps_node();
	MPI_Allreduce(&nsteps_temp,&nsteps_world,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	nsubsteps_total() = nsteps_world;

	MPI_Barrier(MPI_COMM_WORLD);
	output_ctrl->total_timer->stop();


	initializer->finish(particles_old,moments,fields_half);

	output_ctrl->FinalOutput(moments,particles_old,fields_old,fields_next);



//	getchar();


//		getchar();
	//moments[0]->close_plot();
	//fields_next -> close_plot();

	debug_log.print("Total Subcycle Steps: %e\n", (double)nsubsteps_total());
	debug_log.print("Total Outer Picard Steps: %e\n", (double)npiccard_outer());

	output_ctrl->exit_control = 1;
	sleep(5);





}























