/*-------------------------------------------------------------------------*/
/**
	@file	RunData.h
	@author	J. Payne
	@date		04/22/2013
	@brief	Declares the RunData class. Also defines simulation precision,
	several macro calls, constants, and utility functions.


*/
/*--------------------------------------------------------------------------*/
#ifndef Run_DATA_H
#define Run_DATA_H


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include <gnuplot_i.h>
#include <mpi.h>

class RunData
{
public:

	char output_id[16];
	char* SimName; // Name of simulation e.g. IonAcoustic, TwoStream, etc
	char* outpath;
	int runid;
	bool plot_flag;
	bool lo_all; // do the lo solve on all nodes, no field broadcast.


};

class RunStatus
{
public:

	RunStatus()
	{
		npiccard_outer = 0;
		npiccard_outer2 = 0;
		nsubsteps_total = 0;
		nsteps_node = 0;
		istep = 0;
		nsubsteps_total = 0;
	}

	int istep;
	int ipiccard_outer;
	int nsubsteps_op;
	int nsubsteps_step;

	double npiccard_outer_avg;
	int npiccard_outer_last;

	int npiccard_outer;
	double npiccard_outer2;

	double nsubsteps_total;
	double nsteps_node;


};



#endif /* Run_DATA_H */
