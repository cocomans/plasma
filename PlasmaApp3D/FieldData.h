/*-------------------------------------------------------------------------*/
/**
	@file		FieldData.h
	@author	J. Payne
	@date		12/21/2012
	@brief	Declares the FieldData class, a pure virtual class used to encapsulate
	the device specific and optimization specific implementations for field data
	interpolation.


*/
/*--------------------------------------------------------------------------*/
#ifndef FIELD_DATA_H
#define FIELD_DATA_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include "PlasmaData.h"

class MomentSolution;
class ParallelInfo;


typedef struct {
	realkind vals[6]; // Ex, Ey, Ez, Bx, By, Bz
} FieldValues;

/*-------------------------------------------------------------------------*/
/**
	@enum FieldData_deriv
	Enum for choosing derivative of a field data interpolation.

	\var FieldData_deriv::FieldData_deriv_f
	Return the result of the function itself, no derivitive.

	\var FieldData_deriv::FieldData_deriv_dfdx
	Return the 1st derivitive in the x direction.

	\var FieldData_deriv::FieldData_deriv_dfdy
	Return the 1st derivitive in the y direction.

	\var FieldData_deriv::FieldData_deriv_dfdz
	Return the 1st derivitive in the z direction.
*/
/*--------------------------------------------------------------------------*/
enum FieldData_deriv
{
	FieldData_deriv_f = 0,
	FieldData_deriv_dfdx = 1,
	FieldData_deriv_dfdy = 2,
	FieldData_deriv_dfdz = 3
};

/*-------------------------------------------------------------------------*/
/**
	@class FieldData PlasmaData.h
	@author	J. Payne
	@date		12/20/2012
	@brief	Class to control interpolation of the field information
	to the High Order problem

	This class is a pure virtual class that is overloaded with device and
	optimization specific implementations. This parent class merely declares
	the interface for allocating, plotting, interpolating, and setting the
	field data.

	The members declared in this class may be overwritten by children classes
	in order to encompass things such as non-uniform mesh and other optimizations.

	Every child of this class should hold to the following discritization and
	interpolation rules.

	First the Electric and Magnetic fields are chosen to live on a Yee mesh \cite yee1966.
	\image html yeemesh.png

	** Insert Figures here **

	Second, in order for energy conservation the shape functions for the electric
	field interpolation must depend on the interpolation used for the current tally.

	\f{eqnarray}{
	E^{l}\left(\vec{r}_p\right) = \sum_{i,j,k}E^l_{i_l,j_l,k_l}{\cal S}_{E^l}\left(\vec{r}_{i_l,j_l,k_l} - \vec{r}_{p} \right)
	\f}

	With \f$l\f$ indicating the component of \f$E\f$, \f$i_l\f$, \f$j_l\f$, and \f$k_l\f$ are the
	cell index positions where \f$E^l\f$ is stored and are equivilent to \f$f_l = f+1/2\delta_{fl}\f$.
	Where \f${\cal S}_{E^l}\left(\vec{r}_{i_l,j_l,k_l} - \vec{r}_{p} \right)\f$ is
	\f{eqnarray}{
	{\cal S}_{E^l}\left(\vec{r}_{i_l,j_l,k_l} - \vec{r}_{p} \right) = \prod_{m}S^{1+\delta_{lm}}\left(r^m_{i,j,k} - r^m_p \right)
	\f}

	@todo Add more information on Yee mesh, figures, etc. Also add more information on
	shape functions.

*/
/*--------------------------------------------------------------------------*/
class FieldData
{
public:


	__host__ __device__
	virtual ~FieldData(){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Allocate space for the various fields
		@param[in] nx_in,ny_in,nz_in x, y and z field dimensions
		@param[in] nspecies_in number of particle species in simulation.

		@note Never call this alone, only call this from the other allocate
		function.

	*/
	/*--------------------------------------------------------------------------*/
	virtual void allocate(int nx_in, int ny_in, int nz_in, int nspecies_in){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Allocate space for the various fields using plasma data
		@param[in] pdata_in pointer to plasma data, contains all relevent info.

	*/
	/*--------------------------------------------------------------------------*/
	virtual void allocate(PlasmaData* pdata_in){};


//	virtual float3 get_accel(float3 x, float3 v, int3 icell);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Return cell spacing for a given x, y, and z position
		@param[in] x,y,z input position
		@result cell spacing for the given position.

	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual float3 get_divs(realkind x, realkind y, realkind z){};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return global cell spacing
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual float3 get_divs(void){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Return the interpolated electric field value for the input position
		x, y, z, icellx, icelly, icellz
		@param[in] x, y, z Particle position, logical space [0,1]
		@param[in] icellx, icelly, icellz  Particle cell index
		@param[in] icomponent Which component of the electric field to return
		@param[in] ideriv Which derivitive of the electric field to return

		@result Value of the requested component of the electric field at the input
		position.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual realkind intrpE(realkind x, realkind y, realkind z,
				int icellx, int icelly, int icellz,
				const int icomponent, const enum FieldData_deriv ideriv){return 0.0f;};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return the interpolated magnetic field value for the input position
		x, y, z, icellx, icelly, icellz
		@param[in] x, y, z Particle position, logical space [0,1]
		@param[in] icellx, icelly, icellz  Particle cell index
		@param[in] icomponent Which component of the magnetic field to return
		@param[in] ideriv Which derivitive of the magnetic field to return

		@result Value of the requested component of the magnetic field at the input
		position.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual realkind intrpB(realkind x, realkind y, realkind z,
				int icellx, int icelly, int icellz,
				const int icomponent, const enum FieldData_deriv ideriv){return 0.0f;};

	__host__ __device__
	virtual void intrpAccel(realkind x, realkind y, realkind z,
					realkind vx,realkind vy, realkind vz,
					int icellx, int icelly, int icellz,
/* Outputs */		realkind &accelx,realkind& accely, realkind& accelz){};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a reference to the electric field value stored at the input index
		@param[in] ix,iy,iz mesh index of value
		@param[in] icomponent Which component of the electric field to return

		@result Reference to the storage location of the requested component of the
		Electric field.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual realkind& getE(int ix, int iy, int iz,int icomponent){return *((realkind*)NULL);};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a reference to the magnetic field value stored at the input index
		@param[in] ix,iy,iz mesh index of value
		@param[in] icomponent Which component of the magnetic field to return

		@result Reference to the storage location of the requested component of the
		magnetic field.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual realkind& getB(int ix, int iy, int iz,int icomponent){return *((realkind*)NULL);};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Broadcast all values for all fields from the root node to all
		other MPI node.

		@param[in] all_fields contains pointers to all fields that need to be filled from
		broadcasted information.
		@param[in] myInfo parallel information \deprecated Parallel information contained
		in plasma data.

		Basically node 0 broadcasts the values contained in all_fields[0] to the
		all_fields[0] on every other node. The other nodes then populate any additional
		local copies of the field information.
	*/
	/*--------------------------------------------------------------------------*/
	virtual void broadcast(FieldData** all_fields,ParallelInfo* myInfo){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Copy all field information from the source field

		@param[in] src Source field, all values will be copied from this field
	*/
	/*--------------------------------------------------------------------------*/
	virtual void copy_from(FieldData* src){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Evaulate the total potential energy stored in the Electric and
		Magnetic fields

		@result total energy stored in all fields.
	*/
	/*--------------------------------------------------------------------------*/
	virtual realkind evaluate_energy(void){return 0.0f;};

	/*-------------------------------------------------------------------------*/
	/**
		@brief initialize field plots.
	*/
	/*--------------------------------------------------------------------------*/
	virtual void init_plot(){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief reset field plots
	*/
	/*--------------------------------------------------------------------------*/
	virtual void reset_plot(){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief close field plots
	*/
	/*--------------------------------------------------------------------------*/
	virtual void close_plot(){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Plot the stored fields in 2D.

		@param[in] pdata Inputed plasma data
		@param[in] position Position of slice plane
		@param[in] plane Slice plane, 0=xy, 1=xz, 2=yz
		@param[in] Fplot 0=Efield, 1=Bfield
		@param[in] Component to plot 0=x, 1=y, 2=z
	*/
	/*--------------------------------------------------------------------------*/
	virtual void plot(PlasmaData* pdata,int position, int plane,
			const int Fplot,const int icomponent){};

	/** plot handle for gnuplot plotting */
	gnuplot_ctrl* plot_handle;

	/** Face electric field x-component*/
	realkind* Ex;
	/** Face electric field y-component*/
	realkind* Ey;
	/** Face electric field z-component*/
	realkind* Ez;

	/** Edge magnetic field x-component*/
	realkind* Bx;	// Face magnetic field
	/** Edge magnetic field y-component*/
	realkind* By;	// Face magnetic field
	/** Edge magnetic field z-component*/
	realkind* Bz;	// Face magnetic field

	FieldValues* data;

	/** charge to mass ratio for each species */
	realkind* q2m; // (q_s/q_e)/(m_s/m_e)

	/** \name Relevent dimensions */
	//@{
	int nx,ny,nz,nspecies;
	//@}
	/** allocated size for each field */
	int alloc_size;
	/** number of spatial dimensions used in simulation */
	int ndimensions;
	/*! \name Cell spacings
	*/
	//@{
	realkind dx,dy,dz;
	//@}
	/** Simulation information */
	PlasmaData* pdata;

	int FieldType;


};



#endif /* FIELD_DATA_H */
