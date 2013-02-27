/*-------------------------------------------------------------------------*/
/**
	@file		ParticleListCPUSorted.h
	@author	J. Payne
	@date		2/25/2013
	@brief	Declares the ParticleListCPUSorted class \f$\mathcal{P}^t\f$, the CPU version
	child of the ParticleList class.


*/
/*--------------------------------------------------------------------------*/
#ifndef PARTICLE_LIST_CPU_SORTED_H
#define PARTICLE_LIST_CPU_SORTED_H

#include "ParticleListCPU.h"

class FieldData;
class HOMoments;
class ParticleObj;


//static int ParticleList_nfloats = 7;
//static int ParticleList_nints = 3;
/*-------------------------------------------------------------------------*/
/**
	@class ParticleListCPUSorted
	@author	J. Payne
	@date		2/25/2013
	@brief	Pure virtual class to store and manage particle data as well as provide an
	interface to the particle physics object.

	The ParticleList class is a pure virtual class that stores and manages all
	particle information. The children of this class will contain device and
	optimization specific storage and management solutions. Additionally this
	class serves as the interface to the particle physics object, which is
	universal to all devices.


*/
/*--------------------------------------------------------------------------*/
class ParticleListCPUSorted : public ParticleListCPU
{
public:

	__host__ __device__
	ParticleListCPUSorted();
	__host__ __device__
	~ParticleListCPUSorted();

	void allocate(PlasmaData* pdata,int nptcls_in);

	void init(ProblemInitializer* initializer, HOMoments* moments);

	realkind evaluate_energy(PlasmaData* pdata);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Interface to the particle physics object on the CPU
		@param[in] pdata Simulation information
		@param[in] fields Field Data
		@param[in,out] momments HO moments.

		@return total number of particle subcycles executed in push() call
	*/
	/*--------------------------------------------------------------------------*/
	long long int push(PlasmaData* pdata, FieldData* fields, HOMoments* moments);

	template<int VEC_LENGTH,const int nSpatial,const int nVel,const bool iEM>
	long long int pushT(PlasmaData* pdata, FieldData* fields, HOMoments* moments);

	template<int VEC_LENGTH>
	long long int push_interface(PlasmaData* pdata, FieldData* fields, HOMoments* moments);

	void sort(PlasmaData* pdata);

	void copy_stats_from(const ParticleList* list_in);

	void copy_from(const ParticleList* list_in);

	realkind& get_fvalue(int iptcl,int ifloat){return *((*get_float(ifloat)) + iptcl);};

	int& get_ivalue(int iptcl,int iint){return *((*get_int(iint)) + iptcl);};

	const realkind& get_fvalue(int iptcl,int ifloat)const{return *((*get_float(ifloat)) + iptcl);};

	const int& get_ivalue(int iptcl,int iint)const{return *((*get_int(iint)) + iptcl);};

	void CPUfree();

	double4 subcycle_stats(PlasmaData* pdata);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Return the average number of piccard iterations per particle and the standard
		deviation. Also print out piccard statistics.
		@param[in] pdata Simulation information.
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	double4 piccard_stats(PlasmaData* pdata);


	/*-------------------------------------------------------------------------*/
	/**
		@brief Return the cummulative time for the requested timer
	*/
	/*--------------------------------------------------------------------------*/
	double get_cummulative_time(int itimer);

	//void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	// Diagnostic Methods
	void plot_particles(PlasmaData* pdata);

	int num_cores;

	/// Particle index, used for sorting
	int* ptcl_index;





};





#endif /* PARTICLE_LIST_CPU_SORTED_H */
