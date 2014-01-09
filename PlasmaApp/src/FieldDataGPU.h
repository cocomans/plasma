/*-------------------------------------------------------------------------*/
/**
	@file		FieldDataGPU.h
	@author	J. Payne
	@date		1/04/2012
	@brief	Declares the FieldDataGPU, a typedef based on another FieldData class
	depending on whether CUDA support is enabled or not

*/
/*--------------------------------------------------------------------------*/
#ifndef FIELD_DATA_GPU_H
#define FIELD_DATA_GPU_H

#include "FieldDataCPU.h"

// Basically this lets us use a CPU particle list
// in the cases where we want to build without MIC support

#ifndef NO_CUDA
#include "FieldDataGPUSimple.cuh"
typedef FieldDataGPUSimple FieldDataGPU;
#else
typedef FieldDataCPU FieldDataGPU;
#endif



#endif /* FIELD_DATA_GPU_H */
