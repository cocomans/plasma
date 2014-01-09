/*-------------------------------------------------------------------------*/
/**
	@file		HOMomentsCPU.h
	@author	J. Payne
	@date	04/22/2013
	@brief	Declares the HOMomentsCPU Class, a class that stores the HO moments
	tallied from the HO system and used in the LO system.

*/
/*--------------------------------------------------------------------------*/
#ifndef HOMoments_MIC_H
#define HOMoments_MIC_H

#include "HOMomentsCPU.h"


// Basically this lets us use a CPU particle list
// in the cases where we want to build without MIC support

#ifndef NO_MIC
#include "HOMomentsMIC1.h"
typedef HOMomentsMIC1 HOMomentsMIC;
#else
typedef HOMomentsCPU HOMomentsMIC;
#endif


#endif /* HOMoments_MIC_H */
