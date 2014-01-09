/*-------------------------------------------------------------------------*/
/**
	@file		ImplicitPIC.h
	@author	J. Payne
	@date		04/22/2013
	@brief	Declares the ImplicitPIC Class, a class that controls the execution
	flow of the entire simulation.

*/
/*--------------------------------------------------------------------------*/
#ifndef IMPLICIT_PIC_H
#define IMPLICIT_PIC_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <gnuplot_i.h>
#include <mpi.h>
#include "CPUTimer.h"

class PlasmaData;
class ParallelInfo;
class RunStatus;
class NodeHOMoments;
class NodeParticleList;
class NodeFieldData;
class FieldDataCPU;
class OutputCtrl;

class ProblemInitializer;
class FieldSolver;
class OutputPage;



/*-------------------------------------------------------------------------*/
/**
	@class ImplicitPIC ImplicitPIC.h
	@author	J. Payne
	@date		04/22/2013
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
	ImplicitPIC(PlasmaData* pdata_in, FieldSolver* LOsolver_in,
			ProblemInitializer* initializer_in);

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
	FieldDataCPU* fields_old;
	/// Field Data for time \f$\mathcal{F}^{t+1}\f$
	FieldDataCPU* fields_next;
	/// Field Data for time \f$\mathcal{F}^{t+1/2}\f$
	NodeFieldData* fields_half;

	/// Particle Data for time \f$\mathcal{P}^{t}\f$
	NodeParticleList* particles_old;
	/// Particle Data for time \f$\mathcal{P}^{t+1}\f$
	NodeParticleList* particles_next;

	/// HO Moments for times \f$\mathcal{M}^{t}\f$ and \f$\mathcal{M}^{t+1}\f$
	NodeHOMoments* moments;



	/// LO system and field solver
	FieldSolver* LOsolver;
	/// Problem initializer used to setup initial particle distribution and fields
	ProblemInitializer* initializer;

	/// Simulation Information
	PlasmaData* pdata;
	RunStatus* rstatus;
	OutputCtrl* output_ctrl;

	/// Total number of subcycle steps for all nodes
	double& nsubsteps_total();
	/// Total number of subcycle steps for individual node
	double& nsteps_node();
	/// Total number of outer piccard iterations
	int& npiccard_outer();
	/// Squared number of piccard iterations
	double& npiccard_outer2();

	/// MPI rank
	int myid_mpi();
};











#endif /* IMPLICIT_PIC_H */
