#ifndef PARTICLE_LIST_GPU_SORTED_H
#define PARTICLE_LIST_GPU_SORTED_H


#include "FieldDataGPU.h"
#include "HOMomentsGPU.h"



class ParticleListCPU;
class ProblemInitializer;
class ClusterInfo;


//static int ParticleList_nfloats = 7;
//static int ParticleList_nints = 3;

class ParticleListGPUSorted
{
public:

	__host__ __device__
	ParticleListGPUSorted();

	__host__ __device__
	~ParticleListGPUSorted();

	void allocate(PlasmaData* pdata,int nptcls_in);

	void allocate(PlasmaData* pdata_in,ParticleListGPUSorted* _old);



	void init(ProblemInitializer* initializer, HOMomentsGPU* moments);

	realkind evaluate_energy(PlasmaData* pdata);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Tally \f$n_{i,j,k}^{t+1}\f$
		@param[in] pdata Simulation information
		@param[in] momments HO moments.

	*/
	/*--------------------------------------------------------------------------*/
	void accumulate_charge(PlasmaData* pdata, HOMomentsGPU* moments){};

	long long int push(PlasmaData* pdata, FieldDataGPU* fields, HOMomentsGPU* moments);

	long long int push_interface2(PlasmaData* pdata, FieldDataGPU* fields, HOMomentsGPU* moments);

	long long int push_interface(PlasmaData* pdata,
			FieldDataGPU* fields,
			HOMomentsGPU* moments);
	// Methods to iterate over members
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
		case 4:
			result = (int**)&nsubcycles_current;
			break;
		default:
			result = NULL;
			break;
		}

		return result;
	}


	void copy_from(const ParticleListGPUSorted* list_in);

	void copy_from(const ParticleListCPU* list_in);

	void copy_to(ParticleListCPU* list_in);


	__host__ __device__
	realkind& get_fvalue(int iptcl,int ifloat){return *((*get_float(ifloat)) + iptcl);};
	__host__ __device__
	int& get_ivalue(int iptcl,int iint){return *((*get_int(iint)) + iptcl);};
	__host__ __device__
	const realkind& get_fvalue(int iptcl,int ifloat)const{return *((*get_float(ifloat)) + iptcl);};
	__host__ __device__
	const int& get_ivalue(int iptcl,int iint)const{return *((*get_int(iint)) + iptcl);};

	void ReorderData(int* particleIDs,int nptcls_left_old);

	void SetupBlockingInfo(int nptcls_check);



	template<int nSpatial,int nVel,bool iEM>
	void SortedGPUPushH(PlasmaData* 			pdata,
						FieldDataGPU* 			fields,
						HOMomentsGPU* 				moments);

	void CPUfree();

	double get_cummulative_time(int itimer){return 0.0;};

	double4 subcycle_stats(PlasmaData* pdata);

	double4 piccard_stats(PlasmaData* pdata);

	//void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	// Diagnostic Methods
	void plot_particles(PlasmaData* pdata);

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




	/// species of particles stored in this list
	int ispecies;

	/// number of particles stored in this list
	int nptcls;
	/// number of particle slots allocated in this list
	int nptcls_allocated;

	/// Type of device that this list resides on
	int device_type; // 0 = cpu, 1 = gpu


	int* nsubcycles_thread;
	/// Number of subcycles completed for current time step
	int* nsubcycles_current;
	int* particleIDs_original; // originalID = particleIDs_original[currentID]
	//int* particleIDs_original_map; // currentID = particleIDs_original_map[originalID]
	int* particleIDs;

	/// Blocking information for thread blocks
	ClusterInfo* clusters;



	/// Array to store whether or not particle has finished, uses same space as buffer32
	int* ifinished;

	int blocksize;
	int gridsize;
	int nclusters;
	int ClusterSortDim;

	PlasmaData* pdata_d;






};

template<int nSpatial,int nVel> __global__
void write_cluster_ids(PlasmaData* pdata,
					ParticleListGPUSorted particles,
					int* ifinished,
					int nptcls_check);

template<int nSpatial,int nVel,bool iEM> __global__
void GPUBlockPush(PlasmaData* 				pdata,
					FieldDataGPU* 			fields,
					HOMomentsGPU* 				moments,
					ParticleListGPUSorted	particles,
					int* num_subcycles,
					int tally_size);


template<int nSpatial,int nVel,bool iEM>
void SimpleGPUPushH(PlasmaData* 			pdata,
					FieldDataGPU* 			fields,
					HOMomentsGPU* 				moments,
					ParticleListGPUSorted*	particles,
					int* num_subcycles);

void FieldCheckGPU(FieldDataGPU fields);

void MomentCheckGPU(HOMomentsGPU moments);

#endif /* PARTICLE_LIST_GPU_SORTED_H */
