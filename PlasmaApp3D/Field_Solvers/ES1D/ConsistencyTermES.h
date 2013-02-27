//=================================================================================
//
//	AUTHOR: WILLIAM T. TAITANO
//	This is a class for the consistency term
//
//=================================================================================
#ifndef ConsistencyTermES_H
#define ConsistencyTermES_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include "../PlasmaData.h"
#include "../FieldData.h"
#include "../HOMoments.h"

class ParallelInfo;

enum consistencyTerm
{
	consistencyterm_continuity = 0,
	consistencyterm_momentum = 1
};

class ConsistencyTermES {
public:
	//================================================
	//	Constructor
	//================================================
	ConsistencyTermES(PlasmaData* pDataIn);
	//================================================
	//	Destructor
	//================================================
	~ConsistencyTermES();
	//================================================
	//	Member Function: get_val
	//================================================
	realkind& get_val(const int ix, const int iy, const int iz,
				const int ispecies,enum consistencyTerm consistencyterm);
	//================================================
	//	Member Function: consistencyCalcContinuity
	//================================================
	void consistencyCalcContinuity(	int species, HOMoments* curHOMoments,
									HOMoments* oldHOMoments);
	//================================================
	//	Member Function: consistencyCalcMomentum
	//================================================
	void consistencyCalcMomentum(	int species, HOMoments* curHOMoments, HOMoments* oldHOMoments,
								 	FieldData* curFieldData, FieldData* oldFieldData);
	//================================================
	//	Member Variable
	//================================================
	PlasmaData* pData;
	realkind* 	gammaN;
	realkind* 	gammaNU;
};
#endif
