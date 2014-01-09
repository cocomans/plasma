/*-------------------------------------------------------------------------*/
/**
	@file	RunData.h
	@author	J. Payne
	@date		04/22/2013
	@brief	Declares the RunData class. Also defines simulation precision,
	several macro calls, constants, and utility functions.


*/
/*--------------------------------------------------------------------------*/

#ifndef NODEPARTICLELIST_H_
#define NODEPARTICLELIST_H_
#include "ParticleListCPU.h"
#include "ParticleListGPU.cuh"
#include "ParticleListMIC.h"

class ProblemInitializer;
class NodeFieldData;
class CPUTimer;


class NodeParticleList
{
public:

	void allocate(PlasmaData* _pdata);

	void allocate(PlasmaData* _pdata,NodeParticleList* _old);

	void init(ProblemInitializer* initializer,NodeHOMoments* moments);

	long long int push(NodeFieldData* fields,NodeHOMoments* moments);

	void pthreads_push(NodeFieldData* fields,NodeHOMoments* moments,int tid);


	void copy_from(NodeParticleList* src);

//	void subcycle_dist(double** xvals,double** dist);
//
//	void piccard_dist(double** xvals,double** dist);

	int getCPU_gpu(int i);

	int getCPU_mic(int i);

	PlasmaData* pdata;
	ParallelInfo** thread_info;

	long long int* nsubcycles_pushed;


	ParticleListCPU* cpu_particles;
	ParticleListGPU* gpu_particles;
	ParticleListMIC* mic_particles;

	CPUTimer* cpu_push_time;
	CPUTimer* gpu_push_time;


};


typedef struct{
NodeParticleList* particles;
NodeFieldData* fields;
NodeHOMoments* moments;
int tid;
} pthreads_push_data;

#endif /* NODEPARTICLELIST_H_ */
