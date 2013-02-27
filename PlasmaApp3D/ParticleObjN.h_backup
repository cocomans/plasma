#ifndef PARTICLE_OBJ_N_H
#define PARTICLE_OBJ_N_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "typevecN.h"
#include "FieldData.h"
#include "ParticleList.h"

class ParticleList;
class PlasmaData;
class HOMoments;
class CurrentTally;


template<const int N>
class ParticleObjN
{
public:

	__host__ __device__
	ParticleObjN(int* pid_in){pid=pid_in;}

	__host__ __device__
	ParticleObjN(){}

	__host__ __device__
	~ParticleObjN(){}

	__host__ __device__
	ParticleObjN& operator=(const ParticleList& list_in)
	{
		// copy realkind values
		for(int i=0;i<ParticleList_nfloats-1;i++)
		{
			for(int j=0;j<N;j++)
			{
				//(*get_float(i))(j) = *((list_in.get_float(i))[0]+pid[j]);

				(*get_float(i))(j) = list_in.get_fvalue(pid[j],i);
			}


		}

		// Copy int values
		for(int i=0;i<ParticleList_nints;i++)
		{
			for(int j=0;j<N;j++)
			{
				//(*get_int(i))(j) = *((list_in.get_int(i))[0]+pid[j]);

				(*get_int(i))(j) = list_in.get_ivalue(pid[j],i);
			}
		}

		for(int j=0;j<N;j++)
		{
			cluster_id(j) = list_in.cluster_id[pid[j]];
		}

		species = list_in.ispecies;


		return *this;
	}

	void write_back(ParticleList& list_in, const int iptcl)
	{
		// copy realkind values
		for(int i=0;i<ParticleList_nfloats-1;i++)
		{
			//printf("Old Value[%i] = %f\n",iptcl, *((list_in.get_float(i))[0]+pid[iptcl]));
			//printf("New Value[%i] =%f\n",iptcl,(*get_float(i))(iptcl));

			//*((list_in.get_float(i))[0]+pid[iptcl]) = (*get_float(i))(iptcl);

			list_in.get_fvalue(pid[iptcl],i) = (*get_float(i))(iptcl);


		}

		// Copy int values
		for(int i=0;i<ParticleList_nints;i++)
		{
			//*((list_in.get_int(i))[0]+pid[iptcl]) = (*get_int(i))(iptcl);

			list_in.get_ivalue(pid[iptcl],i) = (*get_int(i))(iptcl);
		}

		list_in.cluster_id[pid[iptcl]] = cluster_id(iptcl);


	}

	void copy_in(const ParticleList& list_in, const int iptcl)
	{
		// copy realkind values
		for(int i=0;i<ParticleList_nfloats-1;i++)
		{
			(*get_float(i))(iptcl) = list_in.get_fvalue(pid[iptcl],i);


			//(*get_float(i))(iptcl) = *((list_in.get_float(i))[0]+pid[iptcl]);
		}

		// Copy int values
		for(int i=0;i<ParticleList_nints;i++)
		{
			(*get_int(i))(iptcl) = list_in.get_ivalue(pid[iptcl],i);

			//(*get_int(i))(iptcl) = *((list_in.get_int(i))[0]+pid[iptcl]);
		}

		cluster_id(iptcl) = list_in.cluster_id[pid[iptcl]];

	}

	__attribute__((noinline))
	void push(PlasmaData* pdata,FieldData* fields,
			CurrentTally* currents,typevecN<int,N>& iter,const int nSubcycl_max);

	typevecN<realkind,N> estimate_dtau(PlasmaData* pdata, FieldData* fields);

	__attribute__((noinline))
	void PicardIterate(PlasmaData* pdata, FieldData* fields, CurrentTally*,typevecN<realkind,N>& dtau);

	typevecN<realkind,N> time_till_crossing(typevecN<realkind,N>& v,typevecN<realkind,N>& p,
					typevecN<realkind,N>& accel, typevecN<realkind,N>& dtau0,const realkind scale);

	void accumulate_current(PlasmaData* pdata, CurrentTally* currents,
			const typevecN<realkind,N> x_half,const typevecN<realkind,N> y_half,const typevecN<realkind,N> z_half,
			const typevecN<realkind,N> vx_half,const typevecN<realkind,N> vy_half,const typevecN<realkind,N> vz_half,
			const typevecN<realkind,N>& dtau);

	void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	typevecN<typevecN<realkind,N>,3> calc_div_a(PlasmaData* pdata, FieldData* fields,
			const typevecN<typevecN<realkind,N>,3>& accel);

	typevecN<typevecN<realkind,N>,3> get_accel(FieldData* fields,
								typevecN<realkind,N>& x,
								typevecN<realkind,N>& y,
								typevecN<realkind,N>& z,
								const typevecN<realkind,N>& vx,
								const typevecN<realkind,N>& vy,
								const typevecN<realkind,N>& vz);
	template<enum FieldData_deriv ideriv>
	typevecN<typevecN<realkind,N>,3> get_B(FieldData* fields,const typevecN<float3,N>& x);

	template<enum FieldData_deriv ideriv>
	typevecN<typevecN<realkind,N>,3> get_B(FieldData* fields,
			const typevecN<realkind,N>& x,
			const typevecN<realkind,N>& y,
			const typevecN<realkind,N>& z,
			const typevecN<int,N>& ix_out,
			const typevecN<int,N>& iy_out,
			const typevecN<int,N>& iz_out);

	template<enum FieldData_deriv ideriv>
	typevecN<typevecN<realkind,N>,3> get_E(FieldData* fields,const typevecN<float3,N>& x);

	template<enum FieldData_deriv ideriv>
	typevecN<typevecN<realkind,N>,3> get_E(FieldData* fields,
			const typevecN<realkind,N>& x,
			const typevecN<realkind,N>& y,
			const typevecN<realkind,N>& z,
			const typevecN<int,N>& ix_out,
			const typevecN<int,N>& iy_out,
			const typevecN<int,N>& iz_out);

	bool check_exit(PlasmaData* pdata,const typevecN<int,N>& iter, const int nSubcycl_max);

	void check_boundary();

	// Methods to iterate over members
	__device__ __host__
	typevecN<realkind,N>* get_float(const int& i)
	{
		typevecN<realkind,N>* result;

		switch(i)
		{
		case 0:
			result = &px;
			break;
		case 1:
			result = &py;
			break;
		case 2:
			result = &pz;
			break;
		case 3:
			result = &vx;
			break;
		case 4:
			result = &vy;
			break;
		case 5:
			result = &vz;
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

	__device__ __host__
	typevecN<int,N>* get_int(const int& i)
	{
		typevecN<int,N>* result;

		switch(i)
		{
		case 0:
			result = &ix;
			break;
		case 1:
			result = &iy;
			break;
		case 2:
			result = &iz;
			break;
		default:
			result = NULL;
			break;
		}

		return result;
	}



	typevecN<realkind,N> px; // x position within cell range [0:1]
	typevecN<realkind,N> py; // y position within cell range [0:1]
	typevecN<realkind,N> pz; // z position within cell range [0:1]

	typevecN<realkind,N> vx; // x velocity
	typevecN<realkind,N> vy; // y velocity
	typevecN<realkind,N> vz; // z velocity

	typevecN<int,N> ix; // x cell index
	typevecN<int,N> iy; // y cell index
	typevecN<int,N> iz; // z cell index

	typevecN<realkind,N> dt_finished; // completed portion of the time step

	typevecN<short int,N> cluster_id; // cell cluster index (for sorting)

	int* pid; // location of this particles data in the main list
	int species;


};

template<const int ileft> __attribute__((noinline))
void PushN(PlasmaData* 		pdata,
				 FieldData* 		fields,
				 CurrentTally* 		current,
				 ParticleList* 	plist,
				 int**				iter_array,
				 int**				iptcl,
				 int**				iptcl_new,
				 int&				nptcls_left,
				 int&				nptcl_done,
				 int&				nptcls_process,
				 long long int&				nSubSteps_done);

template<const int ileft> __attribute__((noinline))
void shrink_push(PlasmaData* 		pdata,
				 FieldData* 		fields,
				 CurrentTally* 		current,
				 ParticleList* 	plist,
				 int**				iter_array,
				 int**				iptcl,
				 int**				iptcl_new,
				 int&				nptcls_left,
				 int&				nptcl_done,
				 int&				nptcls_process,
				 long long int&				nSubSteps_done);


template<const int i>
class ploop {
public:
    ploop<i-1> x;
    ParticleObjN<i> particle;
};

template<>
class ploop<1> {
public:
	ParticleObjN<1> particle;
};








#endif /* PARTICLE_OBJ_H */
