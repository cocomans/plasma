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
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include "PlasmaData.h"

class ParallelInfo;

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

/*-------------------------------------------------------------------------*/
/**
	@class HOMoments HOMoments.h
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
class HOMoments
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
	HOMoments(PlasmaData* pdata_in,int device_type_in = 0);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Empty constructor

	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	HOMoments(){};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Empty Destructor.

	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	~HOMoments()
	{
	}

	__host__
	void Free()
	{
		free(charge);
		free(currentx);
		free(currenty);
		free(currentz);
		free(S2xx);
		free(S2xy);
		free(S2xz);
		free(S2xz);
		free(S2yz);
		free(S2zz);
	}

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

	/// GNUplot handle
	gnuplot_ctrl* plot_handle;

	/// type of device this moment is stored on. 0 for CPU 1 for GPU
	int device_type;

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

	/*-------------------------------------------------------------------------*/
	/**
		@brief Reduce all OpenMP copies of the moments to a single copy

		@param[in] tid OpenMP thread ID.

	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void reduce(int tid);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Reduce HO moments across all openMP threads and mpi nodes

		@param[in] all_moments Pointer to list containing all OpenMP copies of moments.
		@param[in] pll_info Parallel configuration information. \deprecated

	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void mpi_reduce(HOMoments** all_moments,ParallelInfo* pll_info);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Copy all moments from the source to this

		@param[in] src Pointer to source moments, all values will be copied from
		 these moments
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void copy_from(HOMoments* src);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Check charge conservation

		@param[in] moments_old Moments from time \f$t\f$, calling moment is at
		time \f$t+1\f$

		This method performs the following charge conservation check.

		\f{eqnarray}{
		CC^{t+1} = \frac{1}{n_x n_y n_z}\sum_{i,j,k,l}\frac{n^{t+1}_{i,j,k,l}-n^{t}_{i,j,k,l}}{\Delta t}
		+ \frac{\overline{nu}^{t+1/2}_{i+1/2,j,k,l}-\overline{nu}^{t+1/2}_{i-1/2,j,k,l}}{\Delta x}
		+ \frac{\overline{nv}^{t+1/2}_{i,j+1/2,k,l}-\overline{nv}^{t+1/2}_{i,j-1/2,k,l}}{\Delta y}
		+ \frac{\overline{nw}^{t+1/2}_{i,j,k+1/2,l}-\overline{nw}^{t+1/2}_{i,j,k-1/2,l}}{\Delta z}
		\f}
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	realkind check_charge(HOMoments* moments_old);
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
	void plot(int position, int plane, int ispecies,enum HOMoments_moment moment);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Simplified interface for other plot function. Plots only electron
		moments

		@param[in] position Placement of slice plane on unused axis.
		@param[in] plane Slice plane to use, 0=xy, 1=xz, 2=yz
		@param[in] moment Moment to plot.

	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void plot(int position, int plane = 0,enum HOMoments_moment moment = HOMoments_charge){plot(position,plane,0,moment);}
	/*-------------------------------------------------------------------------*/
	/**
		@brief Set all elements of all moments to the input value.

		@param[in] val_in value to set all array elements to.

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
	/*-------------------------------------------------------------------------*/
	/**
		@brief Reset moment plot.
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void reset_plot();
	/*-------------------------------------------------------------------------*/
	/**
		@brief Initialize moment plot.
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void init_plot();
	/*-------------------------------------------------------------------------*/
	/**
		@brief Close moment plot.
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void close_plot();



};









#endif /* HOMoments_H */
