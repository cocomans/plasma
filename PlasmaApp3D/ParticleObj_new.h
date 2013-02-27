#ifndef PARTICLE_OBJ_H
#define PARTICLE_OBJ_H

#include "ParticleList.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "typevecN.h"
#include "PlasmaData.h"

class FieldData;
class HOMoments;
class CurrentTally;





template<const int N>
class ParticleObj
{
public:

	__host__ __device__
	~ParticleObj();

	__host__ __device__
	ParticleObj(int* pid_in){pid=pid_in;}



	__host__ __device__
	ParticleObj<N>& operator=(const ParticleList& list_in)
	{
		// copy realkind values
		for(int i=0;i<ParticleList_nfloats-1;i++)
		{
			*get_float(i) = (list_in.get_float(i))[0]+*pid;

		}

		// Copy int values
		for(int i=0;i<ParticleList_nints;i++)
		{
			*get_int(i) = (list_in.get_int(i))[0]+*pid;
		}

		cluster_id = list_in.cluster_id+*pid;



		return *this;
	}

	void push2(void);

	void push(PlasmaData* pdata,FieldData* fields,
			CurrentTally* currents,const int nSubcycl_max)
	{

		dt_finished = 0.0;
		int iter = 0;


		// Begin Subcycle
		while(dt_finished < pdata->dt)
		{

			// Estimate dtau
			realkind dtau = estimate_dtau(pdata,fields);


			// Begin Picard - Crank-Nicholson
			ParticleObj<N> new_ptcl = PicardIterate(pdata,fields,dtau);


			// Check Cell crossing

			// Handle Cell crossing

			// Accumulate Current
			accumulate_current(pdata,currents);



			dt_finished += dtau;

			iter++;

			if(iter >= nSubcycl_max)
				break;
		}

	}

	typevecN<realkind,N> estimate_dtau(PlasmaData* pdata, FieldData* fields);

	ParticleObj<N> PicardIterate(PlasmaData* pdata, FieldData* fields, const typevecN<realkind,N> dtau);

	void accumulate_current(PlasmaData* pdata, CurrentTally* currents);

	void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	typevecN<float3,N> calc_div_a(PlasmaData* pdata, FieldData* fields,const typevecN<float3,N> accel);

	// Methods to iterate over members
	__device__ __host__
	typevecN<realkind,N>* get_float(const int i)
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
			result = NULL
			break;
		}

		return result;
	}

	__device__ __host__
	typevecN<int,N>* get_int(const int i)
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
			result = NULL
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



};




typedef ParticleObj<16> ParticleObjCPU2;




#endif /* PARTICLE_OBJ_H */
