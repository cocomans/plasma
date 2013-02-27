/*-------------------------------------------------------------------------*/
/**
	@file		FieldDataGPU.cuh
	@author	J. Payne
	@date		12/21/2012
	@brief	Declares the FieldDataGPU class, a child of the FieldData class.
	FieldDataGPU handles the gpu specific storage and interpolation for the
	field data on the GPU. Further optimized versions of this class should be
	children of this class.

*/
/*--------------------------------------------------------------------------*/
#ifndef FIELD_DATA_GPU_H
#define FIELD_DATA_GPU_H

#include "FieldData.h"

class MomentSolution;



/*-------------------------------------------------------------------------*/
/**
	@class FieldDataGPU FieldDataGPU.cuh FieldData.h
	@author	J. Payne
	@date		12/20/2012
	@brief	Class to control interpolation of the field information
	to the High Order problem for the GPU version of the solver. This class handles
	the storage, maintainence, and interpolation of the field information for the High
	Order problem on the GPU.

	As a recap from FieldData.h here are the interpolation and storage rules for this
	class.

	First the Electric and Magnetic fields are chosen to live on a Yee mesh \cite yee1966.

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
class FieldDataGPU : public FieldData
{
public:

	__host__ __device__
	~FieldDataGPU()
	{

	}
	__host__
	void FieldFree(void)
	{
		free(Ex);
		free(Ey);
		free(Ez);

		free(Bx);
		free(By);
		free(Bz);

		free(q2m);
	}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Allocate space on the GPU for the various fields
		@param[in] nx_in,ny_in,nz_in x, y and z field dimensions
		@param[in] nspecies_in number of particle species in simulation.

		@note Never call this alone, only call this from the other allocate
		function.

	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void allocate(int nx_in, int ny_in, int nz_in, int nspecies_in);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Interface to allocate using plasma data
		@param[in] pdata_in pointer to plasma data, contains all relevent info.

	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void allocate(PlasmaData* pdata);
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
	realkind intrpE(realkind x, realkind y, realkind z,
				int icellx, int icelly, int icellz,
				const int icomponent, const enum FieldData_deriv ideriv);
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
	realkind intrpB(realkind x, realkind y, realkind z,
				int icellx, int icelly, int icellz,
				const int icomponent, const enum FieldData_deriv ideriv);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Broadcast all values for all fields from the root node to all
		other MPI node.

		@param[in] all_fields contains pointers to all fields that need to be filled from
		broadcasted information.
		@param[in] myInfo parallel information \deprecated Parallel information contained
		in plasma data.

		\note Not implemented for this child class.
	*/
	/*--------------------------------------------------------------------------*/
	void broadcast(FieldData** all_fields,ParallelInfo* myInfo);


	/* Single Precision Functions */
	template<int icomponent, enum FieldData_deriv ideriv> __host__ __device__
	realkind intrpET(realkind x, realkind y, realkind z,
			int icellx, int icelly, int icellz);

	template<int icomponent, enum FieldData_deriv ideriv> __host__ __device__
	realkind intrpBT(realkind x, realkind y, realkind z,
			int icellx, int icelly, int icellz);

//	float3 get_accel(float3 x, float3 v, int3 icell);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a reference to the electric field value stored at the input index
		@param[in] ix,iy,iz mesh index of value
		@param[in] icomponent Which component of the electric field to return

		@result Reference to the storage location of the requested component of the
		Electric field.
	*/
	/*--------------------------------------------------------------------------*/
	template<int icomponent> __host__ __device__
	realkind& getET(int ix, int iy, int iz);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a reference to the magnetic field value stored at the input index
		@param[in] ix,iy,iz mesh index of value
		@param[in] icomponent Which component of the magnetic field to return

		@result Reference to the storage location of the requested component of the
		magnetic field.
	*/
	/*--------------------------------------------------------------------------*/
	template<int icomponent> __host__ __device__
	realkind& getBT(int ix, int iy, int iz);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a reference to the electric field value stored at the input index
		hashed using a z-order curve.
		@param[in] ix,iy,iz mesh index of value
		@param[in] icomponent Which component of the electric field to return

		@result Reference to the storage location of the requested component of the
		electric field.
	*/
	/*--------------------------------------------------------------------------*/
	template<int icomponent> __host__ __device__
	realkind& getETz(int ix, int iy, int iz);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a reference to the magnetic field value stored at the input index
		hashed using a z-order curve.
		@param[in] ix,iy,iz mesh index of value
		@param[in] icomponent Which component of the magnetic field to return

		@result Reference to the storage location of the requested component of the
		magnetic field.
	*/
	/*--------------------------------------------------------------------------*/
	template<int icomponent> __host__ __device__
	realkind& getBTz(int ix, int iy, int iz);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Interface to the templated function getET() or getETz()
		@param[in] ix,iy,iz mesh index of value
		@param[in] icomponent Which component of the electric field to return

		@result Reference to the storage location of the requested component of the
		Electric field.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	realkind& getE(int ix, int iy, int iz,int icomponent);
	/*-------------------------------------------------------------------------*/
	/**
		@brief  Interface to the templated function getBT() or getBTz()
		@param[in] ix,iy,iz mesh index of value
		@param[in] icomponent Which component of the magnetic field to return

		@result Reference to the storage location of the requested component of the
		magnetic field.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	realkind& getB(int ix, int iy, int iz,int icomponent);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Copy all field information from the source field

		@param[in] src Source field, all values will be copied from this field
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void copy_from(FieldData* src);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Evaulate the total potential energy stored in the Electric and
		Magnetic fields

		@result total energy stored in all fields.
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	realkind evaluate_energy(void);
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return cell spacing for a given x, y, and z position
		@param[in] x,y,z input position
		@result cell spacing for the given position.

	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	float3 get_divs(realkind x, realkind y, realkind z)
	{
		return make_float3(dx,dy,dz);
	}
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return global cell spacing
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	float3 get_divs(void)
	{
		return make_float3(dx,dy,dz);
	}
	/*-------------------------------------------------------------------------*/
	/**
		@brief initialize field plots.
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void init_plot();
	/*-------------------------------------------------------------------------*/
	/**
		@brief reset field plots
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void reset_plot();
	/*-------------------------------------------------------------------------*/
	/**
		@brief close field plots
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	void close_plot();
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
	__host__
	void plot(PlasmaData* pdata,int position, int plane,
			const int Fplot,const int icomponent);


/*
	realkind* Ex;	// Face electric field
	realkind* Ey;	// Face electric field
	realkind* Ez;	// Face electric field

	realkind* Bx;	// Face magnetic field
	realkind* By;	// Face magnetic field
	realkind* Bz;	// Face magnetic field

	int nx,ny,nz;
	realkind dx,dy,dz;

*/
};





#endif /* FIELD_DATA_GPU_H */
