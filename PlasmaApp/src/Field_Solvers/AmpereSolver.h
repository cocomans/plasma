#ifndef AMPERE_SOLVER_H
#define AMPERE_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "../cuda_defines.h"
#include <gnuplot_i.h>
#include <mpi.h>
#include "../FieldSolver.h"




class AmpereSolver : public FieldSolver
{
public:

	~AmpereSolver();

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
			NodeHOMoments* moments);

	void solve(PlasmaData* pdata,FieldDataCPU* fields,NodeHOMoments* moments);

	realkind calc_residual(PlasmaData* pdata,
			FieldDataCPU* fields_next,
			FieldDataCPU* fields_old,
			NodeHOMoments* moments);

	FieldDataCPU* fields_old;
	FieldDataCPU* fields_next;

};






#endif /* AMPERE_SOLVER_H */
