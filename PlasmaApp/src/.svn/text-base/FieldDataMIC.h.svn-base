/*-------------------------------------------------------------------------*/
/**
	@file		FieldDataMIC.h
	@author	J. Payne
	@date		04/23/2013
	@brief	Declares the FieldDataMIC, a typedef based on another FieldData class
	depending on whether CUDA support is enabled or not

*/
/*--------------------------------------------------------------------------*/
#ifndef FIELD_DATA_MIC_H
#define FIELD_DATA_MIC_H

#include "FieldDataCPU.h"

// Basically this lets us use a CPU particle list
// in the cases where we want to build without MIC support

#ifndef NO_MIC
#include "FieldDataMICSimple.h"
typedef FieldDataMICSimple FieldDataMIC;
#else
typedef FieldDataCPU FieldDataMIC;
#endif



#endif /* FIELD_DATA_MIC_H */
