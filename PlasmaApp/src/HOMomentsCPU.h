/*-------------------------------------------------------------------------*/
/**
	@file		HOMomentsCPU.h
	@author	J. Payne
	@date	04/22/2013
	@brief	Declares the HOMomentsCPU Class, a class that stores the HO moments
	tallied from the HO system and used in the LO system.

*/
/*--------------------------------------------------------------------------*/
#ifndef HOMoments_CPU_H
#define HOMoments_CPU_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include <gnuplot_i.h>
#include <mpi.h>
#include "PlasmaData.h"
#include "HOMoments.h"

/*-------------------------------------------------------------------------*/
/**
	@class HOMomentsCPU HOMomentsCPU.h
	@author	J. Payne
	@date		1/04/2012
	@brief	Class to store and transport the moments accumulated by the HO system.
	Moments are accumulated by the StressTally, ChargeTally, and CurrentTally
	objects, and stored by the HOMoment object. After a particle push the HOMoments
	class reduces the moments across all nodes.

	Currently moments are stored in a 4D array, 3 spatial dimensions and 1 species
	dimensions.

	\todo Add all components of the stress tensor.

*/
/*--------------------------------------------------------------------------*/
class HOMomentsCPU
{
public:

	/*-------------------------------------------------------------------------*/
	/**
		@brief Full Constructor. Constructs an HOMoments object and allocates
		arrays for the various moments on the specified device.
		@param[in] pdata_in pointer to simulation information
		@param[in] device_type_in device to allocate moment storage on. 0=CPU 1=GPU

	*/
	/*--------------------------------------------------------------------------*/
	HOMomentsCPU(PlasmaData* pdata_in);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Empty constructor

	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	HOMomentsCPU(){};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Empty Destructor.

	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	~HOMomentsCPU()
	{
	}

	__host__
	void Free()
	{
		free(all_data);
	}

	realkind* all_data;

	/// charge \f$\rho^{t+1}\f$ (0th Moment)
	realkind* charge;
	/// x component of current \f$j_x^{t+1/2}\f$ (1st Moment)
	realkind* currentx;
	/// y component of current \f$j_y^{t+1/2}\f$ (1st Moment)
	realkind* currenty;
	/// z component of current \f$j_z^{t+1/2}\f$ (1st Moment)
	realkind* currentz;
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

	/// Simulation information
	PlasmaData* pdata;

	/// type of device this moment is stored on. 0 for CPU 1 for GPU, 2 for MIC
	int device_type;

	int nx,ny,nz,nspecies;

	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a reference to an element of the specified moment
		@param[in] ix,iy,iz cell index of the desired element
		@param[in] ispecies Which species's moment to return.
		@param[in] moment Which moment to return

		@result Reference to one element of the specified moment array.

	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	realkind& get_val(const int ix, const int iy, const int iz,
			const int ispecies,enum HOMoments_moment moment);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a constant reference to an element of the specified moment
		@param[in] ix,iy,iz cell index of the desired element
		@param[in] ispecies Which species's moment to return.
		@param[in] moment Which moment to return

		@result Constant Reference to one element of the specified moment array.

	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	const realkind& get_val(const int ix, const int iy, const int iz,
			const int ispecies,enum HOMoments_moment moment)const;

	__host__ __device__
	enum HOMoments_moment get_enum(int i)
	{

		enum HOMoments_moment result;
		switch(i)
		{
		case 0:
			result = HOMoments_charge;
			break;
		case 1:
			result = HOMoments_currentx;
			break;
		case 2:
			result = HOMoments_currenty;
			break;
		case 3:
			result = HOMoments_currentz;
			break;
		case 4:
			result = HOMoments_S2xx;
			break;
		case 5:
			result = HOMoments_S2xy;
			break;
		case 6:
			result = HOMoments_S2xz;
			break;
		case 7:
			result = HOMoments_S2yy;
			break;
		case 8:
			result = HOMoments_S2yz;
			break;
		case 9:
			result = HOMoments_S2zz;
			break;
		default:
			break;
		}

		return result;
	}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Calculate and return the total particle kinetic energy from
		the stress tensor.

		@result Total system kinetic energy

		This function evaluates the following equation in order to return the
		total plasma kinetic energy.

		\f{eqnarray}{
		KE^{t+1} = \frac{1}{\Delta x\Delta y}\sum_{i,j,k,l}\left(S_{xx,i,j,k,l}^{HO,t+1}
					+S_{yy,i,j,k,l}^{HO,t+1}
					+S_{zz,i,j,k,l}^{HO,t+1}\right)
		\f}

	*/
	/*--------------------------------------------------------------------------*/
	__host__
	realkind evaluate_energy(void);

	__host__
	realkind evaluate_energy(int iS);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Copy all moments from the source to this

		@param[in] src Pointer to source moments, all values will be copied from
		 these moments
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void copy_from(HOMomentsCPU* src);

	__host__
	void copy_to(HOMomentsCPU* dst);

//	/*-------------------------------------------------------------------------*/
//	/**
//		@brief Copy all moments from the source to this
//
//		@param[in] src Pointer to source moments, all values will be copied from
//		 these moments
//	*/
//	/*--------------------------------------------------------------------------*/
//	__host__
//	void copy_from(HOMomentsGPU* src);
//
//	/*-------------------------------------------------------------------------*/
//	/**
//		@brief Copy all moments from the source to this
//
//		@param[in] src Pointer to source moments, all values will be copied from
//		 these moments
//	*/
//	/*--------------------------------------------------------------------------*/
//	__host__
//	void copy_from(HOMomentsMIC* src);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Plot a slice of a specific moment for a specific species

		@param[in] position Placement of slice plane on unused axis.
		@param[in] plane Slice plane to use, 0=xy, 1=xz, 2=yz
		@param[in] ispecies Species of the moments to plot.
		@param[in] moment Moment to plot.

	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void set_vals(realkind val_in);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Apply particle and volume weights to all elements of all moments

		This applies static weighting factors to convert moments from physical
		quantities to weighted densities. Also applies bilinear filtering to all
		moment quantities.

	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void apply_weights(void);



};









#endif /* HOMoments_CPU_H */
