/*-------------------------------------------------------------------------*/
/**
	@file		CurrentTally.h
	@author	J. Payne
	@date		12/21/2012
	@brief	Class to control the tallying of the 1st high order moment, current

	This class is populated with a local current tally, This can encompass the
	entire domain or a subset of the domain. This class should be used for
	both GPU and cpu code.


	@note Payne, 12/21/2012 This class should probably be updated to use a PlasmaData
	based constructor.

*/
/*--------------------------------------------------------------------------*/
#ifndef Current_Tally_H
#define Current_Tally_H

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

/*-------------------------------------------------------------------------*/
/**
	@class CurrentTally CurrentTally.h
	@author	J. Payne
	@date		12/21/2012
	@brief	Class to control the tallying of the 1st high order moment, current

	\paragraph Detailed information:
	This class is populated with a local current tally, This can encompass the
	entire domain or a subset of the domain. This class should be used for
	both GPU and cpu code.

	Essentially this class handles the following equations:
	\f{eqnarray}{
	\overline{nu}^{t+1/2}_{i+1/2,j,k} = \frac{\omega}{\Delta x_{i+1/2}\Delta y_{j} \Delta z_{k}}\sum_p\frac{1}{\Delta t}\sum_{\nu}\Delta \tau^{\nu}v^{\nu+1/2}_{x,p}{\cal S}_{nu}\left(\vec{r}_{i+1/2,j,k} - \vec{r}^{\nu+1/2}_p \right)
	\label{eqn:accum_jx}
	\f}

	\f{eqnarray}{
	\overline{nv}^{t+1/2}_{i,j+1/2,k} = \frac{\omega}{\Delta x_{i}\Delta y_{j+1/2} \Delta z_{k}}\sum_p\frac{1}{\Delta t}\sum_{\nu}\Delta \tau^{\nu}v^{\nu+1/2}_{y,p}{\cal S}_{nv}\left(\vec{r}_{i,j+1/2,k} - \vec{r}^{\nu+1/2}_p \right)
	\label{eqn:accum_jy}
	\f}

	\f{eqnarray}{
	\overline{nw}^{t+1/2}_{i,j,k+1/2} = \frac{\omega}{\Delta x_{i}\Delta y_{j} \Delta z_{k+1/2}}\sum_p\frac{1}{\Delta t}\sum_{\nu}\Delta \tau^{\nu}v^{\nu+1/2}_{z,p}{\cal S}_{nw}\left(\vec{r}_{i,j,k+1/2} - \vec{r}^{\nu+1/2}_p \right)
	\label{eqn:accum_jz}
	\f}
	@note Payne, 12/21/2012 This class should probably be updated to use a PlasmaData
	based constructor.

*/
/*--------------------------------------------------------------------------*/
class CurrentTally
{
public:

	/*-------------------------------------------------------------------------*/
	/**
		@brief Constructs a CurrentTally object
		@param[in] currentx_in,currenty_in,currentz_in Pointers
		 to memory space where x, y, and z current tallies will be stored
		@param[in] dims_in size of the CurrentTally domain: nx, ny, nz
		@param[in] spacingx, Physical x space per cell: dxdi
		@param[in] spacingy, Physical y space per cell: dydi
		@param[in] spacingz, Physical z spcae per cell: dzdi
		@param[in] ndimensions_in Number of spatial dimensions to be used.

	*/
	/*--------------------------------------------------------------------------*/
	FUNCTION_TYPE
	CurrentTally(realkind* currentx_in,
			realkind* currenty_in,
			realkind* currentz_in,
			int3 dims_in,
			realkind spacingx,realkind spacingy,realkind spacingz,
			int ndimensions_in);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Empty Constructor
	*/
	/*--------------------------------------------------------------------------*/
	FUNCTION_TYPE
	CurrentTally();

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
	/*-------------------------------------------------------------------------*/
	/**
		@brief 1D2V tally function

		Tallies the x and y components of the current for only the 1D2V case.
		@param[in] px Particle x positions, logical space [0,1]
		@param[in] vx,vy Particle x,y velocities, physical units.
		@param[in] ix_in Particle cell x-index
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$

	*/
	/*--------------------------------------------------------------------------*/
	FUNCTION_TYPE
	void tally1d2v(const realkind px,
			 const realkind vx,const realkind vy,
			 const int ix_in,
			 const realkind scale);
	/*-------------------------------------------------------------------------*/
	/**
		@brief 1D3V tally function

		Tallies the x, y, and z components of the current for only the 1D3V case.
		@param[in] px Particle x position, logical space [0,1]
		@param[in] vx,vy,vz Particle x,y,z velocities, physical units.
		@param[in] ix_in Particle cell x-index
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$

	*/
	/*--------------------------------------------------------------------------*/
	FUNCTION_TYPE
	void tally1d3v(const realkind px,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,
			 const realkind scale);
	/*-------------------------------------------------------------------------*/
	/**
		@brief 2D2V tally function

		Tallies the x and y components of the current for only the 2D2V case.
		@param[in] px,py Particle x,y positions, logical space [0,1]
		@param[in] vx,vy Particle x,y velocities, physical units.
		@param[in] ix_in,iy_in Particle cell x,y indicies
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$

	*/
	/*--------------------------------------------------------------------------*/
	FUNCTION_TYPE
	void tally2d2v(const realkind px,const realkind py,
			 const realkind vx,const realkind vy,
			 const int ix_in,const int iy_in,
			 const realkind scale);
	/*-------------------------------------------------------------------------*/
	/**
		@brief 2D3V tally function

		Tallies the x, y, and z components of the current for only the 2D3V case.
		@param[in] px,py Particle x,y positions, logical space [0,1]
		@param[in] vx,vy,vz Particle x,y,z velocities, physical units.
		@param[in] ix_in,iy_in Particle cell x,y indicies
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$

	*/
	/*--------------------------------------------------------------------------*/
	FUNCTION_TYPE
	void tally2d3v(const realkind px,const realkind py,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,const int iy_in,
			 const realkind scale);
	/*-------------------------------------------------------------------------*/
	/**
		@brief 3D3V tally function

		Tallies the x, y, and z components of the current for only the 3D3V case.
		@param[in] px,py,pz Particle x,y,z positions, logical space [0,1]
		@param[in] vx,vy,vz Particle x,y,z velocities, physical units.
		@param[in] ix_in,iy_in,iz_in Particle cell x,y,z indicies
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$

	*/
	/*--------------------------------------------------------------------------*/
	FUNCTION_TYPE
	void tally3d3v(const realkind px,const realkind py,const realkind pz,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,const int iy_in,const int iz_in,
			 const realkind scale);

	/*-------------------------------------------------------------------------*/
	/**
		@deprecated Replaced by dimensonality specific cases.
	*/
	/*--------------------------------------------------------------------------*/
	FUNCTION_TYPE
	void tally(const realkind px, const realkind py, const realkind pz,
			 const realkind vx, const realkind vy, const realkind vz,
			 const int ix, const int iy, const int iz,
			 const realkind scale);


	realkind* currentx;
	realkind* currenty;
	realkind* currentz;

	/*! \name Dimensions of the domain covered by this object
	*/
	//@{
	int nx,ny,nz;
	//@}

	/// number of spatial dimensions
	int ndimensions;
	/*! \name Cell spacings
	*/
	//@{
	realkind dx,dy,dz;
	//@}



};









#endif /* Current_Tally_H */
