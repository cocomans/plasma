#ifndef PARTICLE_OBJ_CPU_H
#define PARTICLE_OBJ_CPU_H

#include "ParticleList.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

class PlasmaData;
class FieldData;
class HOMoments;





class ParticleObjCPU : public ParticleObj
{
public:


	ParticleObjCPU(int* pid_in){pid=pid_in;}


	~ParticleObjCPU();


	ParticleObj& operator=(const ParticleList& list_in)
	{
		// copy float values
		for(int i=0;i<ParticleList_nfloats-1;i++)
		{
			*get_float(i) = (list_in.get_float(i))[0][*pid];
		}

		// Copy int values
		for(int i=0;i<ParticleList_nints;i++)
		{
			*get_int(i) = (list_in.get_int(i))[0][*pid];
		}

		cluster_id = list_in.cluster_id[*pid];



		return *this;
	}

	void push(PlasmaData* pdata,FieldData* fields,
			CurrentTally* currents,
			const int nSubcycle_max);

	float estimate_dtau(PlasmaData* pdata, FieldData* fields);

	ParticleObjCPU PicardIterate(PlasmaData* pdata, FieldData* fields, const float dtau);

	void accumulate_current(PlasmaData* pdata, CurrentTally* currents);

	void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	float3 calc_div_a(PlasmaData* pdata, FieldData* fields,const float3& accel);

};









#endif /* PARTICLE_OBJ_CPU_H */
