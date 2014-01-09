/*-------------------------------------------------------------------------*/
/**
	@file		HOMomentsCPU.h
	@author	J. Payne
	@date	04/22/2013
	@brief	Declares the HOMomentsCPU Class, a class that stores the HO moments
	tallied from the HO system and used in the LO system.

*/
/*--------------------------------------------------------------------------*/
#ifndef HOMoments_GPU_H
#define HOMoments_GPU_H



// Basically this lets us use a CPU particle list
// in the cases where we want to build without MIC support

#ifndef NO_CUDA
#include "HOMomentsGPUs.h"
typedef HOMomentsGPUs HOMomentsGPU;
#else
#include "HOMomentsCPU.h"
typedef HOMomentsCPU HOMomentsGPU;
#endif




#endif /* HOMoments_GPU_H */
