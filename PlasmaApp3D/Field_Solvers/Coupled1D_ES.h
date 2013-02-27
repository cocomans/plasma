#ifndef COUPLED1D_ES_H
#define COUPLED1D_ES_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>

#include "../FieldSolver.h"
#include "../PlasmaData.h"
#include "ConsistencyTermES.h"
#include "LOMomentsES.h"
#include "LOSolverSI.h"
#include "../FieldData.h"
#include "../HOMoments.h"
#include "NLResidual.h"



class Coupled1D_ES : public FieldSolver
{
public:

	//==================================================================
	//	Destructor
	//==================================================================

	~Coupled1D_ES();
	//==================================================================
	//	Member Function: init
	//==================================================================

	void init(PlasmaData* pdata,
			FieldData* fields_old,
			FieldData* fields_half,
			FieldData* fields_next,
			HOMoments* moments_old,
			HOMoments* moments_next);
	//==================================================================
	//	Member Function: solve
	//==================================================================
	void solve(	PlasmaData* pdata,
				FieldData* fields_next,
				HOMoments* moments_next);
	//==================================================================
	//	Member Function: calc_residual
	//==================================================================
	realkind calc_residual(	PlasmaData* pData,
							FieldData* fields_next,
							FieldData* fields_old,
							HOMoments* curHOMoments);

	//==================================================================
	//	Member Variables
	//==================================================================
	PlasmaData* 	pData;
	FieldData* 		curInternalField;
	FieldData* 		fields_old;
	FieldData* 		fields_half;
	FieldData* 		fields_next;
	HOMoments* 		moments_old;
	HOMoments* 		moments_next;
	NLResidual* 	residual;
	LOSolverSI* 	solver;



};






#endif /* COUPLED1D_ES_H */
