//========================================================================
//
//	AUTHOR: WILLIAM T. TAITANO
//	This is the header file for the coupled1D_ES solver using JFNK ( I may want to make a virtual function to switch
//	between SI and JFNK at some point...
//
//========================================================================


#ifndef COUPLED1D_ES_JFNK_H_
#define COUPLED1D_ES_JFNK_H_

//====================================
//	General header-files
//====================================
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
//#include <cuda.h>
//#include <cuda_runtime.h>
#include <gnuplot_i.h>
//#include <mpi.h>
//====================================
//	User-defined header-files
//====================================
#include "../FieldSolver.h"
#include "../PlasmaData.h"
#include "ConsistencyTermES.h"
#include "LOMomentsES.h"
#include "LOSolverJFNK.h"
#include "../FieldData.h"
#include "../HOMoments.h"
#include "NLResidual.h"

//====================================
//	Coupled1D_ES_JFNK inherits from FieldSolver
//====================================
class Coupled1D_ES_JFNK: public FieldSolver
{
public:
	//========================================================================
	//
	//	Destructor
	//
	//========================================================================
	~Coupled1D_ES_JFNK();
	//========================================================================
	//
	//	Member Function: init
	//
	//========================================================================

	void 	init(	PlasmaData* pData,
					FieldData* fields_old,
					FieldData* fields_half,
					FieldData* fields_next,
					HOMoments* moments_old,
					HOMoments* moments_next);
	//========================================================================
	//
	//	Member Function (virtual): solve
	//
	//========================================================================
	void	solve(	PlasmaData* pdata,
						FieldData* fields_next,
						HOMoments* moments_next);

	//==================================================================
	//
	//	Member Function: calc_residual
	//
	//==================================================================
	realkind calc_residual(	PlasmaData* pData,
							FieldData* fields_next,
							FieldData* fields_old,
							HOMoments* curHOMoments);
	//========================================================================
	//
	//	Member Variables
	//
	//========================================================================
	PlasmaData* 			pData;
	FieldData* 				curInternalField;
	FieldData* 				fields_old;
	FieldData* 				fields_half;
	FieldData* 				fields_next;
	HOMoments* 			moments_old;
	HOMoments* 			moments_next;
	NLResidual* 			residual;
	LOSolverJFNK* 		solver;
};

#endif /* COUPLED1D_ES_JFNK_H_ */
