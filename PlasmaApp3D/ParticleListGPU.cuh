#ifndef PARTICLE_LIST_GPU_H
#define PARTICLE_LIST_GPU_H

#include "ParticleList.h"


class FieldData;
class HOMoments;
class ParticleObj;
class FieldDataGPU;
class ClusterInfo;


//static int ParticleList_nfloats = 7;
//static int ParticleList_nints = 3;

class ParticleListGPU : public ParticleList
{
public:

	__host__ __device__
	ParticleListGPU();

	__host__ __device__
	~ParticleListGPU();

	void allocate(PlasmaData* pdata,int nptcls_in);

	void init(ProblemInitializer* initializer, HOMoments* moments);

	realkind evaluate_energy(PlasmaData* pdata);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Tally \f$n_{i,j,k}^{t+1}\f$
		@param[in] pdata Simulation information
		@param[in] momments HO moments.

	*/
	/*--------------------------------------------------------------------------*/
	void accumulate_charge(PlasmaData* pdata, HOMoments* moments){};

	long long int push(PlasmaData* pdata, FieldData* fields, HOMoments* moments);

	long long int push_interface2(PlasmaData* pdata, FieldDataGPU* fields, HOMoments* moments);

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


	void copy_from(const ParticleList* list_in);

	void copy_stats_from(const ParticleList* list_in){};

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

	long long int push_interface(PlasmaData* pdata,
			FieldDataGPU* fields,
			HOMoments* moments);

	template<int nSpatial,int nVel,bool iEM>
	void SortedGPUPushH(PlasmaData* 			pdata,
						FieldDataGPU* 			fields,
						HOMoments* 				moments);

	void CPUfree();

	double get_cummulative_time(int itimer){return 0.0;};

	double4 subcycle_stats(PlasmaData* pdata);

	double4 piccard_stats(PlasmaData* pdata);

	//void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	// Diagnostic Methods
	void plot_particles(PlasmaData* pdata);

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

	/// High order moments on the device (object in host memory)
	HOMoments* moments_d;
	/// Pointer to high order moments in constant memory (object in device memory)
	HOMoments* moments_d_cp;

	/// Field data on the GPU (object in host memory)
	FieldDataGPU* fields_d;
	/// Pointer to field data object in constant memory (object in device memory)
	FieldDataGPU* fields_d_cp;
	PlasmaData* pdata_d;

	/// This object in constant memory (object in device memory)
	ParticleListGPU* this_cp;
	int blocksize;
	int gridsize;
	int nclusters;





};

template<int nSpatial,int nVel> __global__
void write_cluster_ids(PlasmaData* pdata,
					ParticleListGPU particles,
					int* ifinished,
					int nptcls_check);

template<int nSpatial,int nVel,bool iEM> __global__
void GPUBlockPush(PlasmaData* 				pdata,
					FieldDataGPU* 			fields,
					HOMoments* 				moments,
					ParticleListGPU	particles,
					int* num_subcycles,
					int tally_size);





#endif /* PARTICLE_LIST_GPU_H */
