#ifndef PARTICLE_OBJ_N_H
#define PARTICLE_OBJ_N_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "PlasmaData.h"

class ParticleList;
class FieldData;
class HOMoments;
class CurrentTally;



class ParticleObj
{
public:


	ParticleObj(int* pid_in){pid=pid_in;}


	virtual ~ParticleObj();


	ParticleObj& operator=(const ParticleList& list_in)
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

	virtual void push(PlasmaData* pdata,FieldData* fields,
			CurrentTally* currents,const int nSubcycl_max);

	virtual realkind estimate_dtau(PlasmaData* pdata, FieldData* fields);

	virtual ParticleObj PicardIterate(PlasmaData* pdata, FieldData* fields, const realkind dtau);

	virtual void accumulate_current(PlasmaData* pdata, CurrentTally* currents);

	virtual void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	virtual float3 calc_div_a(PlasmaData* pdata, FieldData* fields,const float3& accel);

	// Methods to iterate over members
	realkind* get_float(const int& i)
	{
		realkind* result;

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

	int* get_int(const int& i)
	{
		int* result;

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



	realkind px; // x position within cell range [0:1]
	realkind py; // y position within cell range [0:1]
	realkind pz; // z position within cell range [0:1]

	realkind vx; // x velocity
	realkind vy; // y velocity
	realkind vz; // z velocity

	int ix; // x cell index
	int iy; // y cell index
	int iz; // z cell index

	realkind dt_finished; // completed portion of the time step

	short int cluster_id; // cell cluster index (for sorting)

	int* pid; // location of this particles data in the main list



};









#endif /* PARTICLE_OBJ_H */
