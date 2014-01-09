/*-------------------------------------------------------------------------*/
/**
	@file		FieldSolver.h
	@author	J. Payne
	@date		1/04/2012
	@brief	Declares the FieldSolver class, a pure virtual class that takes
	in \f$\rho^{HO,t+1}\f$, \f$\vec{j}^{HO,t+1/2}\f$, and \f$\bar{S}^{HO,t+1}\f$,
	and outputs the updated field values. The children of this class can carry out
	this operation in different ways resulting in constant solvers, simple ampere
	solvers, Coupled HO-LO solvers, EM-solvers, etc...

*/
/*--------------------------------------------------------------------------*/
#ifndef FIELD_SOLVER_H
#define FIELD_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include <gnuplot_i.h>
#include <mpi.h>
#include <omp.h>
#include "PlasmaData.h"
#include "Util/webui/LogStream.h"


class FieldDataCPU;
class NodeFieldData;
class NodeHOMoments;
class ParallelInfo;

/*-------------------------------------------------------------------------*/
/**
	@class FieldSolver FieldSolver.h
	@author	J. Payne
	@date		12/20/2012
	@brief The FieldSolver class is a pure virtual class that takes
	in \f$\rho^{HO,t+1}\f$, \f$\vec{j}^{HO,t+1/2}\f$, and \f$\bar{S}^{HO,t+1}\f$,
	and outputs the updated field values. The children of this class can carry out
	this operation in different ways resulting in constant solvers, simple ampere
	solvers, Coupled HO-LO solvers, EM-solvers, etc...

*/
/*--------------------------------------------------------------------------*/
class FieldSolver
{
public:

	virtual ~FieldSolver(){};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Initializes the field solver
		@param[in] pdata Plasma Data in.
		@param[in] fields_old pointer to field values at current time
		@oaram[in[ fields_half pointer to field values at half time
		@param[in] fields_next pointer to field values at next time
		@param[in] moments_old pointer to current HO Moment values
		@param[in] moments_next pointer to next time HO Moment values.

		For some instances of this class this function will do nothing, for
		others it will initialize the LO moment system.

	*/
	/*--------------------------------------------------------------------------*/
	virtual void init(PlasmaData* pdata,
			FieldDataCPU* fields_old,
			NodeFieldData* fields_half,
			FieldDataCPU* fields_next,
			NodeHOMoments* moments){}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Solve the LO system and return the field values.
		@param[in] pdata Plasma Data in
		@oaram[in,out] fields results of the LO solver are placed here.
		@param[in] moments pointer to HO moments for next time step.


	*/
	/*--------------------------------------------------------------------------*/
	virtual void solve(PlasmaData* pdata,
			FieldDataCPU* fields, //output
			NodeHOMoments* moments){};

	virtual void InitStepSolve(PlasmaData* pdata,
			FieldDataCPU* fields, //output
			NodeHOMoments* moments){};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return the residual of the LO system
		@param[in] pdata Plasma Data in
		@oaram[in] fields_next field values for next time step.
		@param[in] fields_old field values for the current time step.
		@param[in] moments pointer to HO moments for next time step.

		@result residual of the lo system.


	*/
	/*--------------------------------------------------------------------------*/
	virtual realkind calc_residual(PlasmaData* pdata,
			FieldDataCPU* fields_next,
			FieldDataCPU* fields_old,
			NodeHOMoments* moments){return 0.0f;};


	/*-------------------------------------------------------------------------*/
	/**
		@brief Replace the old Moment Solution with the current one,
		do any cleanup or swapping that takes place between time steps.


	*/
	/*--------------------------------------------------------------------------*/
	virtual void update_solution(void){};
};


#endif /* FIELD_SOLVER_H */
