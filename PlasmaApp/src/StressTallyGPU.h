#ifndef Sress_Tally_GPU_H
#define Sress_Tally_GPU_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include "PlasmaData.h"



class StressTallyGPU
{
public:
	__device__
	StressTallyGPU(float* stress_in,
			int3 dims_in,
			realkind spacingx,realkind spacingy,realkind spacingz,
			int ndimensions_in,int nVel_in) :
			S2xx(stress_in),S2xy(stress_in),S2xz(stress_in),
			 S2yy(stress_in),S2yz(stress_in),S2zz(stress_in),
			 nx(dims_in.x),ny(dims_in.y),nz(dims_in.z),
			dx(spacingx),dy(spacingy),dz(spacingz),
			ndimensions(ndimensions_in),nVel(nVel_in)
			{}

	__device__
	StressTallyGPU(float* _S2xx,
			float* _S2xy,
			float* _S2xz,
			float* _S2yy,
			float* _S2yz,
			float* _S2zz,
							 int _nx,int _ny,int _nz,
							 int _ndimensions,int _nvel) :
							 S2xx(_S2xx),S2xy(_S2xy),S2xz(_S2xz),
							 S2yy(_S2yy),S2yz(_S2yz),S2zz(_S2zz),
							 nx(_nx),ny(_ny),nz(_nz),ndimensions(_ndimensions),nVel(_nvel)
	{

	}

	__device__
	StressTallyGPU(float* _S2xx,
				 int _nx,int _ny,int _nz) :
				 S2xx(_S2xx),S2xy(_S2xx),S2xz(_S2xx),
				 S2yy(_S2xx),S2yz(_S2xx),S2zz(_S2xx),
				 nx(_nx),ny(_ny),nz(_nz),ndimensions(1)
	{

	}

	__device__
	StressTallyGPU():S2xx(NULL),S2xy(NULL),S2xz(NULL),
	 S2yy(NULL),S2yz(NULL),S2zz(NULL),nx(0),ny(0),nz(0),ndimensions(1){}

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
	__device__
	void tally1d1v(const realkind px,
			 const realkind vx,
			 const int ix_in,
			 const realkind scale);

	__device__
	void tally1d3v(const realkind px,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,
			 const realkind scale);
	__device__
	void tally2d3v(const realkind px,const realkind py,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,const int iy_in,
			 const realkind scale);

	__device__
	void tally(const realkind px, const realkind py, const realkind pz,
			 const realkind vx, const realkind vy, const realkind vz,
			 const int ix, const int iy, const int iz,
			 const realkind scale);


	/// xx component of Stress \f$S_{xx}^{t+1}\f$ (2nd Moment)
	float* const S2xx;
	/// xy component of Stress \f$S_{xy}^{t+1}\f$ (2nd Moment)
	float* const S2xy;
	/// xz component of Stress \f$S_{xz}^{t+1}\f$ (2nd Moment)
	float* const S2xz;
	/// yy component of Stress \f$S_{yy}^{t+1}\f$ (2nd Moment)
	float* const S2yy;
	/// yz component of Stress \f$S_{yz}^{t+1}\f$ (2nd Moment)
	float* const S2yz;
	/// zz component of Stress \f$S_{zz}^{t+1}\f$ (2nd Moment)
	float* const S2zz;

	int ix0,iy0,iz0;
	const int nx,ny,nz;
	const int ndimensions;
	int nVel;
	realkind dx,dy,dz;


};









#endif /* Sress_Tally_GPU_H */
