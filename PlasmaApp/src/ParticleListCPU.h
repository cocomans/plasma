/*-------------------------------------------------------------------------*/
/**
	@file		ParticleListCPU.h
	@author	J. Payne
	@date		1/09/2012
	@brief	Declares the ParticleListCPU class \f$\mathcal{P}^t\f$, the CPU version
	child of the ParticleList class.


*/
/*--------------------------------------------------------------------------*/
#ifndef PARTICLE_LIST_CPU_H
#define PARTICLE_LIST_CPU_H
#include "PlasmaData.h"

class FieldDataCPU;
class HOMomentsCPU;
class ProblemInitializer;
class CPUTimer;

extern int ParticleListCPU_nfloats;
extern int ParticleListCPU_nints;

//static int ParticleList_nfloats = 7;
//static int ParticleList_nints = 3;
/*-------------------------------------------------------------------------*/
/**
	@class ParticleListCPU ParticleListCPU.h ParticleList.h PlasmaData.h
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
class ParticleListCPU
{
public:

	ParticleListCPU();

	~ParticleListCPU();

	void allocate(PlasmaData* pdata,int nptcls_in);

	void allocate(PlasmaData* pdata,ParticleListCPU* _old);

	void init(ProblemInitializer* initializer, HOMomentsCPU* moments, int offset);

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
	long long int push(PlasmaData* pdata, FieldDataCPU* fields, HOMomentsCPU* moments);

	template<int VEC_LENGTH,const int nSpatial,const int nVel,const bool iEM>
	long long int pushT(PlasmaData* pdata, FieldDataCPU* fields, HOMomentsCPU* moments);

	template<int VEC_LENGTH>
	long long int push_interface(PlasmaData* pdata, FieldDataCPU* fields, HOMomentsCPU* moments);

	void copy_stats_from(const ParticleListCPU* list_in);

	void copy_from(const ParticleListCPU* list_in);

	// Methods to iterate over members

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
	realkind* buffer64;
	int* buffer32;


	/// number of subcycles completed for each particle
	int* num_subcycles;

	/// number of particle-piccard iterations
	realkind* num_piccard;

	/// Original particle index
	int* iptcl_original;



	/// species of particles stored in this list
	int ispecies;

	/// number of particles stored in this list
	int nptcls;
	/// number of particle slots allocated in this list
	int nptcls_allocated;

	/// Type of device that this list resides on
	int device_type; // 0 = cpu, 1 = gpu


	CPUTimer* piccard_timer;
	CPUTimer* accel_timer;
	CPUTimer* tally_timer;
	CPUTimer* crossing_timer;
	CPUTimer* dtau_est_timer;

	CPUTimer* tally_timer2;
	CPUTimer* load_store_timer;

	CPUTimer* push_timer;

	CPUTimer* thread_timer;





};





#endif /* PARTICLE_LIST_CPU_H */
