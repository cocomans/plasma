/*=========================================================================
                                                                                
Copyright (c) 2011, Los Alamos National Security, LLC

All rights reserved.

Copyright 2011. Los Alamos National Security, LLC. 
This software was produced under U.S. Government contract DE-AC52-06NA25396 
for Los Alamos National Laboratory (LANL), which is operated by 
Los Alamos National Security, LLC for the U.S. Department of Energy. 
The U.S. Government has rights to use, reproduce, and distribute this software. 
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,
EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
If software is modified to produce derivative works, such modified software 
should be clearly marked, so as not to confuse it with the version available 
from LANL.
 
Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
-   Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer. 
-   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution. 
-   Neither the name of Los Alamos National Security, LLC, Los Alamos National
    Laboratory, LANL, the U.S. Government, nor the names of its contributors
    may be used to endorse or promote products derived from this software 
    without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
                                                                                
=========================================================================*/

#ifndef NLResidual_h
#define NLResidual_h

#include "../PlasmaData.h"
#include "../HOMoments.h"
#include "../FieldData.h"
#include "ConsistencyTermES.h"
#include "LOMomentsES.h"

using namespace std;

class NLResidual {
public:
	//=================================================================
	//	Constructor
	//=================================================================
	NLResidual(PlasmaData* pDataIn);
	//=================================================================
	//	Destructor
	//=================================================================
	~NLResidual();
	//=================================================================
	// Calculate all residuals and maximum residuals
	//=================================================================
	double calculateNLResidual(
			HOMoments* curHOMoments,
			HOMoments* oldHOMoments,
			LOMomentsES* curLOMomentsES,
			LOMomentsES* oldLOMomentsES,
			FieldData* curFieldData,
			FieldData* oldFieldData,
			ConsistencyTermES* consistencyTermES
        	);
	//=================================================================
	// Non-linear residual function for Ampere Law equation
	//=================================================================
	void nlres_E(
			HOMoments* curHOMoments,
			HOMoments* oldHOMoments,
			LOMomentsES* curLOMomentsES,
			LOMomentsES* oldLOMomentsES,
			FieldData* curFieldData,
			FieldData* oldFieldData);
	//=================================================================
	// Non-linear residual function for continuity equation
	//=================================================================
	void nlres_r(
			int sIndx,
			HOMoments* curHOMoments,
			HOMoments* oldHOMoments,
			LOMomentsES* curLOMomentsES,
			LOMomentsES* oldLOMomentsES,
			FieldData* curFieldData,
			FieldData* oldFieldData,
			ConsistencyTermES* consistencyTermES);
	//=================================================================
	// Non-linear residual function for momentum equation
	//=================================================================
	void nlres_ru(
			int sIndx,
			HOMoments* curHOMoments,
			HOMoments* oldHOMoments,
			LOMomentsES* curLOMomentsES,
			LOMomentsES* oldLOMomentsES,
			FieldData* curFieldData,
			FieldData* oldFieldData,
			ConsistencyTermES* consistencyTermES);
	//=================================================================
	//	Calculate the maximum L2 norm for convergence purpose
	//=================================================================
	double maxL2(
			double* F,
			int size);
	//=================================================================
	//	Calculate the absolute value of maximum value
	//=================================================================
	double absoluteMax(
			double* F,
			int size);
	//=================================================================
	//	Print the residual L2
	//=================================================================
	void printResidual();

	//=================================================================
	//	Member Functions: Get functions
	//=================================================================
	double* 	getF_E()              	{ return this->F_E; }
	double** 	getF_r()             	{ return this->F_r; }
	double** 	getF_ru()            	{ return this->F_ru; }
	double* 	getMax_r()           	{ return this->maxResidual_r; }
	double* 	getMax_ru()           	{ return this->maxResidual_ru; }
	double 		getMax_E()             	{ return this->maxResidual_E; }
	double 		getTotMaxResidual()    	{ return this->totMaxResidual; }
	//=================================================================
	//	Member variables
	//=================================================================
	PlasmaData* pData;

	double** F_r;         // Non-linear residual for continuity equation
	double** F_ru;        // Non-linear residual for moment equation
	double* F_E;          // Non-linear residual for Ampere Law equation

	double* maxResidual_r;
	double* maxResidual_ru;
	double maxResidual_E;
	double totMaxResidual;
};

#endif
