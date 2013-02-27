/*-------------------------------------------------------------------------*/
/**
	@file		ImplicitPIC.h
	@author	J. Payne
	@date		1/08/2012
	@brief	Declares the ImplicitPIC Class, a class that controls the execution
	flow of the entire simulation.

*/
/*--------------------------------------------------------------------------*/
#ifndef IMPLICIT_PIC_H
#define IMPLICIT_PIC_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include "CPUTimer.h"

class PlasmaData;
class ProblemInitializer;
class ParticleList;
class FieldData;
class ParallelInfo;
class FieldSolver;
class HOMoments;
class OutputPage;

double* subcycle_dist(PlasmaData* pdata,ParticleList* plist);
double* piccard_dist(PlasmaData* pdata,ParticleList* plist);

/*-------------------------------------------------------------------------*/
/**
	@class ImplicitPIC ImplicitPIC.h
	@author	J. Payne
	@date		1/08/2012
	@brief	Class to control the execution flow of the implicit algorithm.

	This class calls the initialization of particles and fields, contains the primary
	time step loop, the HO-LO piccard iteration, and final cleanup and output
	routines.

	A flow diagram of this method is shown below.
	\image html sim_overview.png

*/
/*--------------------------------------------------------------------------*/
class ImplicitPIC
{
public:

	/*-------------------------------------------------------------------------*/
	/**
		@brief ImplicitPIC Constructor
		@param[in] pdata_in Simulation information and parameters
		@param[in] LOsolver_in LO solver to be used.
		@param[in] initializer_in ProblemInitializer to use in simulation.
		@param[in] pinfo Parallel configuration information
		@param[in] ndevices_in number of compute devices (deprecated)

		@todo Move particle list and field data allocation and initialization
		into constructor. Get rid of deprecated OpenMP device handeling.

		@todo add support for multiple species per mpi task.

	*/
	/*--------------------------------------------------------------------------*/
	ImplicitPIC(PlasmaData* pdata_in, FieldSolver* LOsolver_in, ProblemInitializer* initializer_in,
			ParallelInfo** pinfo,int ndevices_in);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Empty ImplicitPIC constructor
	*/
	/*--------------------------------------------------------------------------*/
	ImplicitPIC(void){};

	/*-------------------------------------------------------------------------*/
	/**
		\deprecated initialization takes place in constructor

	*/
	/*--------------------------------------------------------------------------*/
	void init(void);

	/*-------------------------------------------------------------------------*/
	/**
		\deprecated replaced by simulate
	*/
	/*--------------------------------------------------------------------------*/
	void run(void);

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
	void simulate(void);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Print out and save timing data and particle statistics

		Gather timing and particle statistics from all nodes and calculate
		global stats.
	*/
	/*--------------------------------------------------------------------------*/
	void save_timing();
	/*-------------------------------------------------------------------------*/
	/**
		@brief Write stats, fields, and moments to disk

	*/
	/*--------------------------------------------------------------------------*/
	void output();
	/*-------------------------------------------------------------------------*/
	/**
		@brief Compute particle subcycle statistics

		@result Page of strings to be written to disk
	*/
	/*--------------------------------------------------------------------------*/
	OutputPage get_subcycle_dists();
	/*-------------------------------------------------------------------------*/
	/**
		@brief Compute particle piccard iteration statistics

		@result Page of strings to be written to disk
	*/
	/*--------------------------------------------------------------------------*/
	OutputPage get_piccard_dists();

	std::string calc_timing(double perf_weight,const char* legend);


	/// Field Data for time \f$\mathcal{F}^{t}\f$
	FieldData* fields_old;
	/// Field Data for time \f$\mathcal{F}^{t+1}\f$
	FieldData* fields_next;
	/// Field Data for time \f$\mathcal{F}^{t+1/2}\f$
	FieldData** fields_half;
	/// Particle Data for time \f$\mathcal{P}^{t}\f$
	ParticleList** particles_old;
	/// Particle Data for time \f$\mathcal{P}^{t+1}\f$
	ParticleList** particles_next;
	/// Particle Data for time \f$\mathcal{M}^{t+1}\f$
	HOMoments**	moments;
	/// Particle Data for time \f$\mathcal{M}^{t}\f$
	HOMoments*	moments_old;

	/// Simulation Information
	PlasmaData* pdata;
	/// Parallel Configuration Information
	ParallelInfo** pll_info;

	/// LO system and field solver
	FieldSolver* LOsolver;
	/// Problem initializer used to setup initial particle distribution and fields
	ProblemInitializer* initializer;

	/// Timer for particle push
	CPUTimer push_timer;
	/// Timer for entire HO system
	CPUTimer HOSolve_timer;
	/// Timer for MPI Communication
	CPUTimer Comm_timer;
	/// Timer for total run time
	CPUTimer total_timer;
	/// Timer for LO system solve
	CPUTimer LOSolve_timer;
	/// Timer for single time step
	CPUTimer step_timer;
	/// Total number of subcycle steps for all nodes
	long long int nsubsteps_total;
	/// Total number of subcycle steps for individual node
	long long int nsteps_node;
	/// Total number of outer piccard iterations
	int npiccard_outer;
	/// Squared number of piccard iterations
	double npiccard_outer2;

	/// MPI rank
	int myid_mpi;
	/// number of compute devices
	int ndevices;
	/// number of OpenMP threads to use
	int nthreads;
	/// number of MPI tasks
	int n_nodes;
};











#endif /* IMPLICIT_PIC_H */
