#ifndef PARTICLE_LIST_GPU_SIMPLE_H
#define PARTICLE_LIST_GPU_SIMPLE_H

#include "ParticleListGPU.cuh"


class FieldData;
class HOMoments;
class ParticleObj;


//static int ParticleList_nfloats = 7;
//static int ParticleList_nints = 3;

class ParticleListGPUSimple : public ParticleListGPU
{
public:
	__host__ __device__
	ParticleListGPUSimple(){};
	__host__ __device__
	~ParticleListGPUSimple(){};

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

	long long int push_interface3(PlasmaData* pdata, FieldDataGPU* fields, HOMoments* moments);

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

	void copy_from(const ParticleList* list_in);

	__host__ __device__
	realkind* get_float_ptr(const int i)
	const
	{

		realkind* result;
		switch(i)
		{
		case 0:
			result = (realkind*)px;
			break;
		case 1:
			result = (realkind*)py;
			break;
		case 2:
			result = (realkind*)pz;
			break;
		case 3:
			result = (realkind*)vx;
			break;
		case 4:
			result = (realkind*)vy;
			break;
		case 5:
			result = (realkind*)vz;
			break;
		case 6:
			result = (realkind*)dt_finished;
			break;
		case 7:
			result = (realkind*)buffer;
			break;
		default:
			result = NULL;
			break;
		}

		return result;
	}

	__host__ __device__
	int* get_int_ptr(const int i)
	const
	{
		int* result;

		switch(i)
		{
		case 0:
			result = (int*)ix;
			break;
		case 1:
			result = (int*)iy;
			break;
		case 2:
			result = (int*)iz;
			break;
		default:
			result = NULL;
			break;
		}

		return result;
	}

	__host__ __device__
	realkind& get_fvalue(int iptcl,int ifloat){return *((*get_float(ifloat)) + iptcl);};
	__host__ __device__
	int& get_ivalue(int iptcl,int iint){return *((*get_int(iint)) + iptcl);};
	__host__ __device__
	const realkind& get_fvalue(int iptcl,int ifloat)const{return *((*get_float(ifloat)) + iptcl);};
	__host__ __device__
	const int& get_ivalue(int iptcl,int iint)const{return *((*get_int(iint)) + iptcl);};

	void CPUfree();

	//void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	// Diagnostic Methods
	void plot_particles(PlasmaData* pdata);





};





#endif /* PARTICLE_LIST_GPU_SIMPLE_H */
