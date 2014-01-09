/*-------------------------------------------------------------------------*/
/**
	@file		ChargeTallyCPU.h
	@author	J. Payne
	@date		04/29/2013
	@brief	Class to control the tallying of the 0th high order moment

	This class is populated with a local charge tally, This can encompass the
	entire domain or a subset of the domain. This class should be used for
	both GPU and cpu code.

	@todo Implement spatial and velocity dimension specific versions of the tally
	method.

	@note Payne, 12/20/2012 This class should probably be updated to use a PlasmaData
	based constructor.

*/
/*--------------------------------------------------------------------------*/
#ifndef Charge_Tally_CPU_H
#define Charge_Tally_CPU_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "PlasmaData.h"

/*-------------------------------------------------------------------------*/
/**
	@class ChargeTally ChargeTally.h
	@author	J. Payne
	@date		12/20/2012
	@brief	Class to control the tallying of the 0th high order moment

	This class is populated with a local charge tally, This can encompass the
	entire domain or a subset of the domain. This class should be used for
	both GPU and CPU code.
	This class encompases the following equation:\f{eqnarray}{
	n^{t+1}_{i,j,k} 	= \frac{\omega}{\Delta x_i \Delta y_j \Delta z_k}\sum_{p}{\cal S}_n \left(\vec{r}_{i,j,k} - \vec{r}^{t+1}_{p} \right)
	\label{eqn:accum_rho}
	\f}

	@todo Implement spatial and velocity dimension specific versions of the tally
	method.

	@note Payne, 12/20/2012 This class should probably be updated to use a PlasmaData
	based constructor.

*/
/*--------------------------------------------------------------------------*/
class ChargeTallyCPU
{
public:

	/*-------------------------------------------------------------------------*/
	/**
		@brief Constructs a ChargeTally object
		@param[in] charge_in Pointer to memory space where tallies will be stored
		@param[in] dims_in size of the ChargeTally domain: nx, ny, nz
		@param[in] spacingx, Physical x space per cell: dxdi
		@param[in] spacingy, Physical y space per cell: dydi
		@param[in] spacingz, Physical z spcae per cell: dzdi
		@param[in] ndimensions_in Number of spatial dimensions to be used.

	*/
	/*--------------------------------------------------------------------------*/
	ChargeTallyCPU(realkind* charge_in,
			int3 dims_in,
			realkind spacingx,realkind spacingy,realkind spacingz,
			int ndimensions_in);


	/*-------------------------------------------------------------------------*/
	/**
		@brief Constructs an empty ChargeTally object
	*/
	/*--------------------------------------------------------------------------*/
	ChargeTallyCPU();


	/*-------------------------------------------------------------------------*/
	/**
		@brief 1D1V tally function

		Tallies the x-component of the current for only the 1D1V case.
		@param[in] px Particle x position, logical space [0,1]
		@param[in] ix_in Particle cell x-index
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$

	*/
	/*--------------------------------------------------------------------------*/
	void tally1d(const realkind px,
			 const int ix_in,
			 const realkind scale);

	/*-------------------------------------------------------------------------*/
	/**
		@brief 2D2V tally function

		Tallies the x and y components of the current for only the 2D2V case.
		@param[in] px,py Particle x,y positions, logical space [0,1]
		@param[in] ix_in,iy_in Particle cell x,y indicies
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$

	*/
	/*--------------------------------------------------------------------------*/
	void tally2d(const realkind px,const realkind py,
			 const int ix_in,const int iy_in,
			 const realkind scale);
	/*-------------------------------------------------------------------------*/
	/**
		@brief 3D3V tally function

		Tallies the x, y, and z components of the current for only the 3D3V case.
		@param[in] px,py,pz Particle x,y,z positions, logical space [0,1]
		@param[in] ix_in,iy_in,iz_in Particle cell x,y,z indicies
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$

	*/
	/*--------------------------------------------------------------------------*/
	void tally3d(const realkind px,const realkind py,const realkind pz,
			 const int ix_in,const int iy_in,const int iz_in,
			 const realkind scale);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Tally the 0th moment using the inputed particle properties
		@param[in] px Particle x position, logical space [0,1]
		@param[in] py Particle y position, logical space [0,1]
		@param[in] pz Particle z position, logical space [0,1]
		@param[in] ix Particle cell x-index
		@param[in] iy Particle cell y-index
		@param[in] iz Particle cell z-index
		@param[in] scale Particle \f$\Delta\tau^k / \Delta t\f$


	*/
	/*--------------------------------------------------------------------------*/
	void tally(const realkind px, const realkind py, const realkind pz,
			 const int ix, const int iy, const int iz,
			 const realkind scale);


	/// pointer to memory space in which tallies will be stored
	realkind* charge;

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









#endif /* Charge_Tally_CPU_H */
