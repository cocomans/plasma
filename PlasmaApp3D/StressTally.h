#ifndef Sress_Tally_H
#define Sress_Tally_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "PlasmaData.h"

#ifdef GPU_CODE
#define FUNCTION_TYPE __attribute__((device))
#else
#define FUNCTION_TYPE __attribute__((host,device))
#endif

class StressTally
{
public:
	FUNCTION_TYPE
	StressTally(realkind* stress_in,
			int3 dims_in,
			realkind spacingx,realkind spacingy,realkind spacingz,
			int ndimensions_in);
	FUNCTION_TYPE
	StressTally();

	/*-------------------------------------------------------------------------*/
	/**
		@brief 1D1V tally function

		Tallies the x-component of the current for only the 1D1V case.
		@param[in] px Particle x position, logical space [0,1]
		@param[in] vx Particle x velocity, physical units.
		@param[in] ix_in Particle cell x-index
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$

	*/
	/*--------------------------------------------------------------------------*/
	FUNCTION_TYPE
	void tally1d1v(const realkind px,
			 const realkind vx,
			 const int ix_in,
			 const realkind scale);
	FUNCTION_TYPE
	void tally(const realkind px, const realkind py, const realkind pz,
			 const realkind vx, const realkind vy, const realkind vz,
			 const int ix, const int iy, const int iz,
			 const realkind scale);


	realkind* stress;

	int nx,ny,nz;
	int nptcls;
	int ndimensions;
	realkind dx,dy,dz;


};









#endif /* Charge_Tally_H */
