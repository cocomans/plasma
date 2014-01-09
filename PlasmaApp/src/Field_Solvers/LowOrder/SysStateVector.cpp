/*
 * SysStateVector.cpp
 *
 *  Created on: Dec 19, 2013
 *      Author: payne
 */
#include "SysStateVector.h"
#include "StateVectorT.h"
#include "Mesh.h"
namespace MultiScale
{

SysStateVector::SysStateVector(SimParams* _params)
{
	moments = HoLo::StateVectorGen(_params);
	particles = new HoLo::HighOrder::NodeParticleList();
	stress = new HoLo::StressTensor();

	// Allocate them
	moments->Allocate(_params->mesh);
	particles->Allocate(_params);
	stress->Allocate(_params);

}

// This is so that we can share timers
SysStateVector::SysStateVector(SimParams* _params,SysStateVector* _in)
{
	moments = HoLo::StateVectorGen(_params);
	particles = new HoLo::HighOrder::NodeParticleList();
	stress = new HoLo::StressTensor();

	// Allocate them
	moments->Allocate(_params->mesh,_in->moments);
	particles->Allocate(_params,_in->particles);
	stress->Allocate(_params,_in->stress);

}

} /* namespace MultiScale */
