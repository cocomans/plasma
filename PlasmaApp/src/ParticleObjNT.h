#ifndef PARTICLE_OBJ_NT_H
#define PARTICLE_OBJ_NT_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "FieldData.h"
#include "typevecN.h"
#include "PlasmaData.h"
#ifndef GPU_CODE
#include "Util/CPUTimer.h"
#include <omp.h>
#endif






template<const int nSpatial,const int nVel>
class BlockLimiter;


template<const int N,const int nSpatial,const int nVel,const bool iEM,
		class ParticleListT,class FieldDataT,class CurrentTallyT,class PExitCheckT,
		class BoundaryCheckT>
class ParticleObjNT
{
public:

	__host__ __device__
	ParticleObjNT(int* pid_in){pid=pid_in;printf("Warning, no exit checker given\n");}

	__host__ __device__
	ParticleObjNT(int* pid_in,PExitCheckT* _exit_check):pid(pid_in),exit_checker(_exit_check){}

	__host__ __device__
	ParticleObjNT(){}

	__host__ __device__
	~ParticleObjNT(){}

	__device__
	ParticleObjNT& operator=(const ParticleListT& list_in)
	{

		for(int i=0;i<N;i++)
		{
			copy_in(list_in,i);
		}


		species = list_in.ispecies;


		return *this;
	}

	__device__
	void write_back_all(ParticleListT& list_in)
	{

		// write back the position values
		for(int i=0;i<nSpatial;i++)
		{
			for(int j=0;j<N;j++)
			list_in.get_fvalue(pid[j],i) = position(i)(j);

			for(int j=0;j<N;j++)
			list_in.get_ivalue(pid[j],i) = iposition(i)(j);
		}

		// write back the velocity values
		for(int i=0;i<nVel;i++)
		{
			for(int j=0;j<N;j++)
			list_in.get_fvalue(pid[j],i+3) = velocity(i)(j);
		}

		// write back dt_finished
		for(int j=0;j<N;j++)
		 list_in.get_fvalue(pid[j],6) = dt_finished(j);

		for(int j=0;j<N;j++)
		list_in.num_piccard[pid[j]] += npiccard(j);
//		list_in.num_piccard2[pid[iptcl]] += npiccard2(iptcl);
	}

	__device__
	void write_back(ParticleListT& list_in, const int iptcl)
	{

		// write back the position values
		for(int i=0;i<nSpatial;i++)
		{
			list_in.get_fvalue(pid[iptcl],i) = position(i)(iptcl);
			list_in.get_ivalue(pid[iptcl],i) = iposition(i)(iptcl);
		}

		// write back the velocity values
		for(int i=0;i<nVel;i++)
		{
			list_in.get_fvalue(pid[iptcl],i+3) = velocity(i)(iptcl);
		}

		// write back dt_finished
		 list_in.get_fvalue(pid[iptcl],6) = dt_finished(iptcl);

		list_in.num_piccard[pid[iptcl]] += npiccard(iptcl);
//		list_in.num_piccard2[pid[iptcl]] += npiccard2(iptcl);
//#ifndef GPU_CODE
//		printf("particle[%i] on cpu: %f, %i, %f\n",pid[iptcl],position(0)(iptcl),iposition(0)(iptcl),velocity(0)(iptcl));
//#else
//		printf("particle[%i] on gpu: %f, %i, %f\n",pid[iptcl],position(0)(iptcl),iposition(0)(iptcl),velocity(0)(iptcl));
//#endif
	}

	__device__
	void copy_in_all(const ParticleListT& list_in)
	{

		//printf("getting particle information for particle %i\n",pid[iptcl]);

		// copy in the position values
		for(int i=0;i<nSpatial;i++)
		{

			for(int j=0;j<N;j++)
			{
				position(i)(j)= list_in.get_fvalue(pid[j],i);
			}

			for(int j=0;j<N;j++)
			{
				iposition(i)(j)= list_in.get_ivalue(pid[j],i);
			}

//			position(i)(iptcl) = list_in.get_fvalue(pid[iptcl],i);
//			iposition(i)(iptcl) = list_in.get_ivalue(pid[iptcl],i);

		}

		// copy in the velocity values
		for(int i=0;i<nVel;i++)
		{
			for(int j=0;j<N;j++)
				velocity(i)(j) = list_in.get_fvalue(pid[j],i+3);
		}

		// copy in dt_finished
		for(int j=0;j<N;j++)
		dt_finished(j) = list_in.get_fvalue(pid[j],6);

		for(int j=0;j<N;j++)
		ifinished(j) = 0;
		for(int j=0;j<N;j++)
		npiccard(j) = 0;
//		npiccard2(iptcl) = 0;

	//	isubcycle(iptcl) = 0;

	}

	__device__
	void copy_in(const ParticleListT& list_in, const int iptcl)
	{

		//printf("getting particle information for particle %i\n",pid[iptcl]);
		int ptcl = pid[iptcl];
		// copy in the position values
		for(int i=0;i<nSpatial;i++)
		{

			realkind tempfloat = list_in.get_fvalue(ptcl,i);
			int tempint = list_in.get_ivalue(ptcl,i);
			position(i)(iptcl) = tempfloat;
			iposition(i)(iptcl) = tempint;

//			position(i)(iptcl) = list_in.get_fvalue(pid[iptcl],i);
//			iposition(i)(iptcl) = list_in.get_ivalue(pid[iptcl],i);

		}

		// copy in the velocity values
		for(int i=0;i<nVel;i++)
		{
			velocity(i)(iptcl) = list_in.get_fvalue(pid[iptcl],i+3);
		}

		// copy in dt_finished
		dt_finished(iptcl) = list_in.get_fvalue(pid[iptcl],6);

		ifinished(iptcl) = 0;

		npiccard(iptcl) = 0;

		species = list_in.ispecies;

//		printf("particle[%i] on cpu: %f, %i, %f\n",ptcl,position(0)(iptcl),iposition(0)(iptcl),velocity(0)(iptcl));

//		npiccard2(iptcl) = 0;

	//	isubcycle(iptcl) = 0;

	}

	__device__
	void copy_in_gpu(const ParticleListT& list_in, const int iptcl)
	{

		//printf("getting particle information for particle %i\n",pid[iptcl]);
		int ptcl = pid[iptcl];
		// copy in the position values
		for(int i=0;i<nSpatial;i++)
		{

//			realkind tempfloat = list_in.get_fvalue(ptcl,i);
//			int tempint = list_in.get_ivalue(ptcl,i);
//			position(i)(iptcl) = tempfloat;
//			iposition(i)(iptcl) = tempint;

			position(i)(iptcl) = list_in.get_fvalue(pid[iptcl],i);
			iposition(i)(iptcl) = list_in.get_ivalue(pid[iptcl],i);

		}

		// copy in the velocity values
		for(int i=0;i<nVel;i++)
		{
			velocity(i)(iptcl) = list_in.get_fvalue(pid[iptcl],i+3);
		}

		// copy in dt_finished
		dt_finished(iptcl) = list_in.get_fvalue(pid[iptcl],6);

		npiccard(iptcl) = 0;
		ifinished(iptcl) = 0;

		species = list_in.ispecies;
//		npiccard2(iptcl) = 0;

//		printf("particle[%i] on gpu: %f, %i, %f\n",ptcl,position(0)(iptcl),iposition(0)(iptcl),velocity(0)(iptcl));

	//	isubcycle(iptcl) = 0;

	}

	__device__ __attribute__((noinline))
	void push(PlasmaData* pdata,FieldDataT* fields,
			CurrentTallyT* currents,typevecN<int,N>& iter,const int nSubcycl_max);
	__device__
	typevecN<realkind,N> estimate_dtau(PlasmaData* pdata, FieldDataT* fields);

	__device__ __attribute__((noinline))
	void PicardIterate(PlasmaData* pdata, FieldDataT* fields, CurrentTallyT* currents,typevecN<realkind,N>& dtau0);

	__device__ __attribute__((noinline))
	void PicardIterate2(PlasmaData* pdata, FieldDataT* fields,
			CurrentTallyT* currents,typevecN<realkind,N>& dtau0);

	__device__
	typevecN<realkind,N> time_till_crossing(typevecN<realkind,N>& v,typevecN<realkind,N>& p,
					 typevecN<realkind,N>& dtau0,const realkind scale);
	__device__ __attribute__((noinline))
	typevecN<realkind,N> time_till_crossing2(
									typevecN<realkind,N>& vin_half,
									typevecN<realkind,N>& pin,
									typevecN<realkind,N>& dtau0,typevecN<realkind,N>& dtau_cur,
									const realkind scale);
	__device__
	void accumulate_current(PlasmaData* pdata, CurrentTallyT* currents,
			const typevecN<typevecN<realkind,N>,nSpatial> position_half,
			const typevecN<typevecN<realkind,N>,nVel> velocity_half,
			const typevecN<realkind,N>& dtau);

	__device__
	typevecN<typevecN<realkind,N>,3> calc_div_a(PlasmaData* pdata, FieldDataT* fields,
			const typevecN<typevecN<realkind,N>,3>& accel);
	__device__
	typevecN<typevecN<realkind,N>,nVel> get_accel(FieldDataT* fields,
								const typevecN<typevecN<realkind,N>,nSpatial>& x,
								const typevecN<typevecN<realkind,N>,nVel>& v);

	__device__
	typevecN<typevecN<realkind,N>,nVel> get_accel2(FieldDataT* fields,
											const typevecN<typevecN<realkind,N>,nSpatial>& x,
											typevecN<realkind,N> dtau);

	template<enum FieldData_deriv ideriv> __device__
	typevecN<typevecN<realkind,N>,nVel> get_B(FieldDataT* fields,
			const typevecN<typevecN<realkind,N>,nSpatial>& x,
			const typevecN<typevecN<int,N>,nSpatial>& ix);

	template<enum FieldData_deriv ideriv> __device__
	typevecN<typevecN<realkind,N>,nVel> get_E(FieldDataT* fields,
			const typevecN<typevecN<realkind,N>,nSpatial>& x,
			const typevecN<typevecN<int,N>,nSpatial>& ix);
	__device__
	bool check_exit(PlasmaData* pdata,const typevecN<int,N>& iter, const int nSubcycl_max);

	__device__
	void check_boundary();

	// Methods to iterate over members
	__device__
	typevecN<realkind,N>* get_float(const int& i)
	{
		typevecN<realkind,N>* result;

		switch(i)
		{
		case 0:
			result = &position(0);
			break;
		case 1:
			result = &position(1);
			break;
		case 2:
			result = &position(2);
			break;
		case 3:
			result = &velocity(0);
			break;
		case 4:
			result = &velocity(1);
			break;
		case 5:
			result = &velocity(2);
			break;
		case 6:
			result = &dt_finished;
			break;
		default:
			result = NULL;
			break;
		}

		return result;
	}

	__device__
	typevecN<int,N>* get_int(const int& i)
	{
		typevecN<int,N>* result;

		switch(i)
		{
		case 0:
			result = &iposition(0);
			break;
		case 1:
			result = &iposition(1);
			break;
		case 2:
			result = &iposition(2);
			break;
		default:
			result = NULL;
			break;
		}

		return result;
	}



	typevecN<typevecN<realkind,N>,nSpatial> position; // position within cell range [0:1]

	typevecN<typevecN<realkind,N>,nVel> velocity; // particle velocities


	typevecN<typevecN<int,N>,nSpatial> iposition; // position within cell range [0:1]

	typevecN<realkind,N> dt_finished; // completed portion of the time step

	typevecN<int,N> ifinished;

//	typevecN<realkind,N> isubcycle; // completed portion of the time step

	typevecN<int,N> npiccard; // number of piccard iterations
//	typevecN<int,N> npiccard2; // number of piccard iterations squared

	int* pid; // location of this particles data in the main list
	int species;

	/// Flexible exit check conditions
	PExitCheckT* exit_checker;

	PlasmaData* pdatal;

#ifndef GPU_CODE

	CPUTimer* piccard_timer;
	CPUTimer* accel_timer;
	CPUTimer* tally_timer;
	CPUTimer* crossing_timer;
	CPUTimer* dtau_est_timer;
#endif


};

template<const int ileft,const int nSpatial,const int nVel,const bool iEM,
		class ParticleListT,class FieldDataT,class CurrentTallyT,class PExitCheckT,
		class BoundaryCheckT> __attribute__((noinline)) __host__
void PushNT(PlasmaData* 		pdata,
				 FieldDataT* 		fields,
				 CurrentTallyT* 		current,
				 ParticleListT* 	plist,
				 int**				iter_array,
				 int**				iptcl,
				 int**				iptcl_new,
				 int&				nptcls_left,
				 int&				nptcl_done,
				 int&				nptcls_process,
				 long long int&				nSubSteps_done);

template<const int ileft,const int nSpatial,const int nVel,const bool iEM,
		class ParticleListT,class FieldDataT,class CurrentTallyT,class PExitCheckT,
		class BoundaryCheckT> __attribute__((noinline)) __host__
void shrink_pushT(PlasmaData* 		pdata,
				 FieldDataT* 		fields,
				 CurrentTallyT* 		current,
				 ParticleListT* 	plist,
				 int**				iter_array,
				 int**				iptcl,
				 int**				iptcl_new,
				 int&				nptcls_left,
				 int&				nptcl_done,
				 int&				nptcls_process,
				 long long int&				nSubSteps_done);









#endif /* PARTICLE_OBJ_H */
