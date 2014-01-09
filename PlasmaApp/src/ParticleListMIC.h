/*-------------------------------------------------------------------------*/
/**
	@file		ParticleListMIC.h
	@author	J. Payne
	@date		1/09/2012
	@brief	Declares the ParticleListMIC type, which is a MIC based
	Particle list if compiled with MIC support, or a CPU particle list
	without MIC support


*/
/*--------------------------------------------------------------------------*/
#ifndef PARTICLE_LIST_MIC_H
#define PARTICLE_LIST_MIC_H

#include "ParticleListCPU.h"


// Basically this lets us use a CPU particle list
// in the cases where we want to build without MIC support

#ifndef NO_MIC
#include "ParticleListMICSimple.h"
typedef ParticleListMICSimple ParticleListMIC;
#else
typedef ParticleListCPU ParticleListMIC;
#endif

#endif /* PARTICLE_LIST_MIC_H */


