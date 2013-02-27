/*-------------------------------------------------------------------------*/
/**
	@file		ParticleList.h
	@author	J. Payne
	@date		1/09/2012
	@brief	Declares the pure virtual ParticleList class \f$\mathcal{P}^t\f$.


*/
/*--------------------------------------------------------------------------*/
#ifndef PARTICLE_LIST_H
#define PARTICLE_LIST_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
//#include <cuda.h>
//#include <cuda_runtime.h>
//#include <gnuplot_i.h>
#include "PlasmaData.h"



class FieldData;
class HOMoments;
class ParticleObj;
class ProblemInitializer;
class CPUTimer;

template<int N>
class ParticleObjN;


static int ParticleList_nfloats = 8;
static int ParticleList_nints = 3;

/*-------------------------------------------------------------------------*/
/**
	@class ParticleList ParticleList.h
	@author	J. Payne
	@date		1/09/2012
	@brief	Pure virtual class to store and manage particle data as well as provide an
	interface to the particle physics object.

	The ParticleList class is a pure virtual class that stores and manages all
	particle information. The children of this class will contain device and
	optimization specific storage and management solutions. Additionally this
	class serves as the interface to the particle physics object, which is
	universal to all devices.


*/
/*--------------------------------------------------------------------------*/
class ParticleList
{
public:
	/*-------------------------------------------------------------------------*/
	/**
		@brief Empty Constructor
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	ParticleList(){}
	/*-------------------------------------------------------------------------*/
	/**
		@brief Destructor
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual ~ParticleList(){};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Allocate a particle list using simulation information and
		a requested number of particles.
		@param[in] pdata Simulation information
		@param[in] nptcls_in requested number of particles for allocation
	*/
	/*--------------------------------------------------------------------------*/
	virtual void allocate(PlasmaData* pdata,int nptcls_in){};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Interface to the particle physics object
		@param[in] pdata Simulation information
		@param[in] fields Field Data
		@param[in,out] momments HO moments.

		@return total number of particle subcycles executed in push() call
	*/
	/*--------------------------------------------------------------------------*/
	virtual long long int push(PlasmaData* pdata, FieldData* fields, HOMoments* moments){return 0;};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Tally \f$n_{i,j,k}^{t+1}\f$
		@param[in] pdata Simulation information
		@param[in] momments HO moments.

	*/
	/*--------------------------------------------------------------------------*/
	virtual void accumulate_charge(PlasmaData* pdata, HOMoments* moments){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Initialize particle distribution using a Problem initializer and tally
		the initial HO moments.
		@param[in] initializer ProblemInitializer used to set the particle distribution
		@param[in,out] momments initial HO moments
	*/
	/*--------------------------------------------------------------------------*/
	virtual void init(ProblemInitializer* initializer, HOMoments* moments){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Sum up the kinetic energy of all particles in the list
		@param[in] pdata Simulation information in.
	*/
	/*--------------------------------------------------------------------------*/
	virtual realkind evaluate_energy(PlasmaData* pdata){return 0.0f;};

	// Diagnostic Methods
	/*-------------------------------------------------------------------------*/
	/**
		@brief Plot a sample of the particle distribution in x/vx space.
		@param[in] pdata Simulation information.

		\deprecated Use HOMoments evaluate energy.
	*/
	/*--------------------------------------------------------------------------*/
	virtual void plot_particles(PlasmaData* pdata){};

	virtual void copy_stats_from(const ParticleList* list_in){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Replace contents of this list with those in list_in
		@param[in] list_in Source particle list.
	*/
	/*--------------------------------------------------------------------------*/
	virtual void copy_from(const ParticleList* list_in){};

	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a reference to one float component of a particle
		@param[in] iptcl Particle ID
		@param[in] ifloat Float component to return.

		@result A reference to a float component of the particle iptcl.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual realkind& get_fvalue(int iptcl,int ifloat){return *(realkind*)NULL;};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a reference to one integer component of a particle
		@param[in] iptcl Particle ID
		@param[in] iint integer component to return.

		@result A reference to a integer component of the particle iptcl.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual int& get_ivalue(int iptcl,int iint){return *(int*)NULL;};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a constant reference to one float component of a particle
		@param[in] iptcl Particle ID
		@param[in] ifloat Float component to return.

		@result A constant reference to a float component of the particle iptcl.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual const realkind& get_fvalue(int iptcl,int ifloat)const{return *(realkind*)NULL;};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a constant reference to one integer component of a particle
		@param[in] iptcl Particle ID
		@param[in] iint integer component to return.

		@result A constant reference to a integer component of the particle iptcl.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	virtual const int& get_ivalue(int iptcl,int iint)const{return *(int*)NULL;};
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return the average number of subcycles per particle and the standard
		deviation. Also print out subcycle statistics.
		@param[in] pdata Simulation information.
	*/
	/*--------------------------------------------------------------------------*/
	virtual double4 subcycle_stats(PlasmaData* pdata){}
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return the average number of piccard iterations per particle and the standard
		deviation. Also print out piccard statistics.
		@param[in] pdata Simulation information.
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	virtual double4 piccard_stats(PlasmaData* pdata){};

	//template<int N>
	//ParticleList& operator=(const ParticleObjN<N>& particle){};

	// Methods to iterate over members
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a pointer to the requested float* member of the ParticleList
		@param[in] i member to return

		@result pointer to the requested float* member.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	realkind** get_float(const int i)
	const
	{

		realkind** result;
		switch(i)
		{
		case 0:
			result = (realkind**)&px;
			break;
		case 1:
			result = (realkind**)&py;
			break;
		case 2:
			result = (realkind**)&pz;
			break;
		case 3:
			result = (realkind**)&vx;
			break;
		case 4:
			result = (realkind**)&vy;
			break;
		case 5:
			result = (realkind**)&vz;
			break;
		case 6:
			result = (realkind**)&dt_finished;
			break;
		case 7:
			result = (realkind**)&num_piccard;
			break;
		case 8:
			result = (realkind**)&buffer;
			break;
		default:
			result = NULL;
			break;
		}

		return result;
	}
	/*-------------------------------------------------------------------------*/
	/**
		@brief Return a pointer to the requested integer* member of the ParticleList
		@param[in] i member to return

		@result pointer to the requested integer* member.
	*/
	/*--------------------------------------------------------------------------*/
	__host__ __device__
	int** get_int(const int i)
	const
	{
		int** result;

		switch(i)
		{
		case 0:
			result = (int**)&ix;
			break;
		case 1:
			result = (int**)&iy;
			break;
		case 2:
			result = (int**)&iz;
			break;
		case 3:
			result = (int**)&num_subcycles;
			break;
		default:
			result = NULL;
			break;
		}

		return result;
	}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Return the cummulative time for the requested timer
	*/
	/*--------------------------------------------------------------------------*/
	virtual double get_cummulative_time(int itimer){return 0.0;};

	/// x position within cell range [0:1]
	realkind* px;
	/// y position within cell range [0:1]
	realkind* py;
	/// z position within cell range [0:1]
	realkind* pz;
	/// x velocity
	realkind* vx;
	/// y velocity
	realkind* vy;
	/// z velocity
	realkind* vz;

	/// x cell index
	int* ix;
	/// y cell index
	int* iy;
	/// z cell index
	int* iz;

	/// completed portion of the time step
	realkind* dt_finished;
	/// cell cluster index (for sorting)
	int* cluster_id;
	/// buffer array used for sorting and various sums
	realkind* buffer;


	/// number of subcycles completed for each particle
	int* num_subcycles;

	/// number of particle-piccard iterations
	realkind* num_piccard;
	/// squared number of piccard iterations
	realkind* num_piccard2;

	/// species of particles stored in this list
	int ispecies;

	/// number of particles stored in this list
	int nptcls;
	/// number of particle slots allocated in this list
	int nptcls_allocated;

	/// Type of device that this list resides on
	int device_type; // 0 = cpu, 1 = gpu

	/// GNUPlot plotting handle
	gnuplot_ctrl* plot;

	CPUTimer* piccard_timer;
	CPUTimer* accel_timer;
	CPUTimer* tally_timer;
	CPUTimer* crossing_timer;
	CPUTimer* dtau_est_timer;

	CPUTimer* tally_timer2;
	CPUTimer* load_store_timer;

	CPUTimer* push_timer;

};









#endif /* PARTICLE_LIST_H */
