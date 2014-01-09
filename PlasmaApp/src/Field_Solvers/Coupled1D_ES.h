#ifndef COUPLED1D_ES_H
#define COUPLED1D_ES_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "../cuda_defines.h"
#include <gnuplot_i.h>
#include <mpi.h>

#include "../FieldSolver.h"
#include "../PlasmaData.h"
#include "ConsistencyTermES.h"
#include "LOMomentsES.h"
#include "LOSolverSI.h"
#include "../FieldDataCPU.h"
#include "../HOMomentsCPU.h"
#include "NLResidual.h"



class Coupled1D_ES : public FieldSolver
{
public:

	//==================================================================
	//	Destructor
	//==================================================================

	~Coupled1D_ES();
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
	void init(PlasmaData* pdata,
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
	void solve(PlasmaData* pdata,
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
	realkind calc_residual(PlasmaData* pdata,
			FieldDataCPU* fields_next,
			FieldDataCPU* fields_old,
			NodeHOMoments* moments){return 0.0f;};


	/*-------------------------------------------------------------------------*/
	/**
		@brief Replace the old Moment Solution with the current one,
		do any cleanup or swapping that takes place between time steps.


	*/
	/*--------------------------------------------------------------------------*/
	void update_solution(void){};

	//==================================================================
	//	Member Variables
	//==================================================================
	PlasmaData* 	pData;
	FieldDataCPU* 		curInternalField;
	FieldDataCPU* 		fields_old;
	FieldDataCPU* 		fields_half;
	FieldDataCPU* 		fields_next;
	HOMomentsCPU* 		moments_old;
	HOMomentsCPU* 		moments_next;
	NLResidual* 	residual;
	LOSolverSI* 	solver;



};






#endif /* COUPLED1D_ES_H */
