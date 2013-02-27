#ifndef PARTICLE_LIST_CPU_AOS_H
#define PARTICLE_LIST_CPU_AOS_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <string.h>
#include "ParticleListCPU.h"



//static int ParticleList_nfloats = 7;
//static int ParticleList_nints = 3;

struct CPUParticleData
{
	realkind px;
	realkind py;
	realkind pz;

	realkind vx;
	realkind vy;
	realkind vz;

	realkind dt_finished;

	int ix;
	int iy;
	int iz;
};

class ParticleListCPUAoS : public ParticleListCPU
{
public:

	~ParticleListCPUAoS()
	{
		nptcls_allocated = 0;
		nptcls = 0;
	}

	void allocate(PlasmaData* pdata,int nptcls_in);

	void copy_from(const ParticleList* list_in);

	void copy_from(const ParticleListCPUAoS* list_in);

	realkind& get_fvalue(int iptcl,int ifloat)
	{

		realkind* result = (realkind*)NULL;
		switch(ifloat)
		{
		case 0:
			return ((particles+iptcl)->px);
			break;
		case 1:
			return ((particles+iptcl)->py);
			break;
		case 2:
			return ((particles+iptcl)->pz);
			break;
		case 3:
			return ((particles+iptcl)->vx);
			break;
		case 4:
			return ((particles+iptcl)->vy);
			break;
		case 5:
			return ((particles+iptcl)->vz);
			break;
		case 6:
			return ((particles+iptcl)->dt_finished);
			break;
		default:
			result = NULL;
			break;
		}

		return *result;
	};

	int& get_ivalue(int iptcl,int iint)
	{

		int* result = (int*)NULL;
		switch(iint)
		{
		case 0:
			return ((particles+iptcl)->ix);
			break;
		case 1:
			return ((particles+iptcl)->iy);
			break;
		case 2:
			return ((particles+iptcl)->iz);
			break;
		default:
			result = NULL;
			break;
		}

		return *result;
	};

	const realkind& get_fvalue(int iptcl,int ifloat)
	const
	{

		realkind* result = (realkind*)NULL;
		switch(ifloat)
		{
		case 0:
			return ((particles+iptcl)->px);
			break;
		case 1:
			return ((particles+iptcl)->py);
			break;
		case 2:
			return ((particles+iptcl)->pz);
			break;
		case 3:
			return ((particles+iptcl)->vx);
			break;
		case 4:
			return ((particles+iptcl)->vy);
			break;
		case 5:
			return ((particles+iptcl)->vz);
			break;
		case 6:
			return ((particles+iptcl)->dt_finished);
			break;
		default:
			result = NULL;
			break;
		}

		return *result;
	};

	const int& get_ivalue(int iptcl,int iint)
	const
	{

		int* result = (int*)NULL;
		switch(iint)
		{
		case 0:
			return ((particles+iptcl)->ix);
			break;
		case 1:
			return ((particles+iptcl)->iy);
			break;
		case 2:
			return ((particles+iptcl)->iz);
			break;
		default:
			result = NULL;
			break;
		}

		return *result;
	};

	void CPUfree();

	//void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	// Diagnostic Methods
//	void plot_particles(PlasmaData* pdata);

	CPUParticleData* particles;





};





#endif /* PARTICLE_LIST_CPU_AOS_H */
