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
			int ndimensions_in,int nVel_in) :
			S2xx(stress_in),nx(dims_in.x),ny(dims_in.y),nz(dims_in.z),
			dx(spacingx),dy(spacingy),dz(spacingz),
			ndimensions(ndimensions_in),nVel(nVel_in)
			{}

	FUNCTION_TYPE
	StressTally(realkind* _S2xx,
							 realkind* _S2xy,
							 realkind* _S2xz,
							 realkind* _S2yy,
							 realkind* _S2yz,
							 realkind* _S2zz,
							 int _nx,int _ny,int _nz,
							 int _ndimensions,int _nvel) :
							 S2xx(_S2xx),S2xy(_S2xy),S2xz(_S2xz),
							 S2yy(_S2yy),S2yz(_S2yz),S2zz(_S2zz),
							 nx(_nx),ny(_ny),nz(_nz),ndimensions(_ndimensions),nVel(_nvel)
	{

	}

	FUNCTION_TYPE
	StressTally(realkind* _S2xx,
				 int _nx,int _ny,int _nz) :
				 S2xx(_S2xx),
				 nx(_nx),ny(_ny),nz(_nz),ndimensions(1)
	{

	}
	FUNCTION_TYPE
	StressTally():nx(0),ny(0),nz(0),ndimensions(1){}

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
	void tally1d3v(const realkind px,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,
			 const realkind scale);
	FUNCTION_TYPE
	void tally2d3v(const realkind px,const realkind py,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,const int iy_in,
			 const realkind scale);

	FUNCTION_TYPE
	void tally(const realkind px, const realkind py, const realkind pz,
			 const realkind vx, const realkind vy, const realkind vz,
			 const int ix, const int iy, const int iz,
			 const realkind scale);


	/// xx component of Stress \f$S_{xx}^{t+1}\f$ (2nd Moment)
	realkind* S2xx;
	/// xy component of Stress \f$S_{xy}^{t+1}\f$ (2nd Moment)
	realkind* S2xy;
	/// xz component of Stress \f$S_{xz}^{t+1}\f$ (2nd Moment)
	realkind* S2xz;
	/// yy component of Stress \f$S_{yy}^{t+1}\f$ (2nd Moment)
	realkind* S2yy;
	/// yz component of Stress \f$S_{yz}^{t+1}\f$ (2nd Moment)
	realkind* S2yz;
	/// zz component of Stress \f$S_{zz}^{t+1}\f$ (2nd Moment)
	realkind* S2zz;

	int ix0,iy0,iz0;
	const int nx,ny,nz;
	const int ndimensions;
	int nVel;
	realkind dx,dy,dz;


};









#endif /* Charge_Tally_H */
