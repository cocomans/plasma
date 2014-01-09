/*-------------------------------------------------------------------------*/
/**
	@file		HOMoments.h
	@author	J. Payne
	@date		1/04/2012
	@brief	Declares the HOMoments Class, a class that stores the HO moments
	tallied from the HO system and used in the LO system.

*/
/*--------------------------------------------------------------------------*/
#ifndef HOMoments_H
#define HOMoments_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include <gnuplot_i.h>
#include <mpi.h>

/*-------------------------------------------------------------------------*/
/**
	@enum HOMoments_moment
	Enum for choosing a moment to return

	\var HOMoments_moment::HOMoments_charge
	Return the value for charge

	\var HOMoments_moment::HOMoments_currentx
	Return the x component of the current

	\var HOMoments_moment::HOMoments_currenty
	Return the y component of the current

	\var HOMoments_moment::HOMoments_currentz
	Return the z component of the current

	\var HOMoments_moment::HOMoments_S2xx
	Return the xx component of the stress tensor

	\var HOMoments_moment::HOMoments_S2xy
	Return the xy component of the stress tensor

	\var HOMoments_moment::HOMoments_S2xz
	Return the xz component of the stress tensor

	\var HOMoments_moment::HOMoments_S2yy
	Return the yy component of the stress tensor

	\var HOMoments_moment::HOMoments_S2yz
	Return the yz component of the stress tensor

	\var HOMoments_moment::HOMoments_S2zz
	Return the zz component of the stress tensor

*/
/*--------------------------------------------------------------------------*/
enum HOMoments_moment
{
	HOMoments_charge = 0,
	HOMoments_currentx = 1,
	HOMoments_currenty = 2,
	HOMoments_currentz = 3,
	HOMoments_S2xx = 4,
	HOMoments_S2xy = 5,
	HOMoments_S2xz = 6,
	HOMoments_S2yy = 7,
	HOMoments_S2yz = 8,
	HOMoments_S2zz = 9,
	HOMoments_currentxyz = 10
};




#endif /* HOMoments_H */
