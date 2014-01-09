/*
 * NodeParticleList.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: payne
 */

#include "omp.h"
#include "mpi.h"
#include "PlasmaData.h"
#include "NodeHOMoments.h"
#include "NodeParticleList.h"
#include "HOMomentsCPU.h"
#include "NodeFieldData.h"
#include "ParallelInfo.h"
#include "pthread.h"
#include "Util/CPUTimer.h"




void NodeParticleList::allocate(PlasmaData* _pdata)
{
	pdata = _pdata;
	thread_info = pdata->thread_info;

	int n_cpu_lists = pdata->nspecies+pdata->node_info->nGPU+pdata->node_info->nMIC;
	nsubcycles_pushed = (long long int*)malloc(n_cpu_lists*sizeof(long long int));
	omp_set_num_threads(pdata->node_info->nCPU);

	cpu_particles = (ParticleListCPU*)malloc(n_cpu_lists*sizeof(ParticleListCPU));

// This only sets up the HOMoments for the CPU, the other CPU HOMoments for
// the GPU and MIC must be set up when those device specific versions are set up
	printf("Allocating CPU particles\n");

	for(int i=0;i<n_cpu_lists;i++)
	{
		int idevice = thread_info[i] -> device_type;
		cpu_particles[i] = *(new ParticleListCPU());
		cpu_particles[i].allocate(pdata,pdata->nptcls_device[idevice]);
		cpu_particles[i].ispecies = thread_info[i] -> ispecies;
	}



	printf("Copying Particles to GPUn");

	if(pdata->node_info->nGPU > 0)
	{
		gpu_particles = (ParticleListGPU*)malloc(pdata->node_info->nGPU*sizeof(ParticleListGPU));

		for(int i=0;i<pdata->node_info->nGPU;i++)
		{
			CUDA_SAFE_CALL(cudaSetDevice(thread_info[pdata->node_info->nspecies+i]->gpu_info->igpu));
			gpu_particles[i] = *(new ParticleListGPU());
			gpu_particles[i].allocate(pdata,pdata->nptcls_device[1]);
		}
	}

	if(pdata->node_info->nMIC > 0)
	{
		mic_particles = (ParticleListMIC*)malloc(pdata->node_info->nMIC*sizeof(ParticleListMIC));

		for(int i=0;i<pdata->node_info->nMIC;i++)
		{
			mic_particles[i] = *(new ParticleListMIC());
			mic_particles[i].allocate(pdata,pdata->nptcls_device[2]);
		}
	}

	cpu_push_time = (CPUTimer*)malloc(pdata->nspecies*sizeof(CPUTimer));
	gpu_push_time = (CPUTimer*)malloc(pdata->node_info->nGPU*sizeof(CPUTimer));;

	for(int i=0;i<pdata->nspecies;i++)
		cpu_push_time[i] = *(new CPUTimer());

	for(int i=0;i<pdata->node_info->nGPU;i++)
		gpu_push_time[i] = *(new CPUTimer());


}

void NodeParticleList::allocate(PlasmaData* _pdata,NodeParticleList* _old)
{

	pdata = _pdata;
	thread_info = pdata->thread_info;
	*this = *_old;

	int n_cpu_lists = pdata->nspecies+pdata->node_info->nGPU+pdata->node_info->nMIC;

	nsubcycles_pushed = (long long int*)malloc(n_cpu_lists*sizeof(long long int));

	omp_set_num_threads(pdata->node_info->nCPU);

	cpu_particles = (ParticleListCPU*)malloc(n_cpu_lists*sizeof(ParticleListCPU));

// This only sets up the HOMoments for the CPU, the other CPU HOMoments for
// the GPU and MIC must be set up when those device specific versions are set up
	printf("Allocating CPU particles\n");

	for(int i=0;i<n_cpu_lists;i++)
	{
		cpu_particles[i] = *(new ParticleListCPU());
		cpu_particles[i].allocate(pdata,_old->cpu_particles+i);
		cpu_particles[i].ispecies = (_old->cpu_particles+i)->ispecies;
	}


	printf("Copying Particles to GPUn");

	if(pdata->node_info->nGPU > 0)
	{
		gpu_particles = (ParticleListGPU*)malloc(pdata->node_info->nGPU*sizeof(ParticleListGPU));

		for(int i=0;i<pdata->node_info->nGPU;i++)
		{
			CUDA_SAFE_CALL(cudaSetDevice(thread_info[pdata->node_info->nspecies+i]->gpu_info->igpu));
			gpu_particles[i] = *(new ParticleListGPU());
			gpu_particles[i].allocate(pdata,_old->gpu_particles+i);
		}
	}

	if(pdata->node_info->nMIC > 0)
	{
		mic_particles = (ParticleListMIC*)malloc(pdata->node_info->nMIC*sizeof(ParticleListMIC));

		for(int i=0;i<pdata->node_info->nMIC;i++)
		{
			mic_particles[i] = *(new ParticleListMIC());
			mic_particles[i].allocate(pdata,_old->mic_particles+i);
		}
	}

}

void NodeParticleList::init(ProblemInitializer* problem,NodeHOMoments* moments)
{
	int n_cpu_lists = pdata->nspecies+pdata->node_info->nGPU+pdata->node_info->nMIC;

	moments -> set_vals(0);

//#pragma omp for
	int ioffset;

	for(int i=0;i<n_cpu_lists;i++)
	{

		if(i < pdata->nspecies)
			ioffset = 0;
		else
			ioffset = i - pdata->nspecies;
		cpu_particles[i].init(problem,moments->cpu_moments+i,ioffset);

//		if(i < pdata->nspecies)
//			ioffset = 0;
//		else
//			ioffset += pdata->nptcls_device[pdata->thread_info[i]->device_type];
	}


	if(pdata->node_info->nGPU > 0)
	{
		printf("Initializing GPU particles\n");
//#pragma omp for
		for(int i=0;i<pdata->node_info->nGPU;i++)
		{
			CUDA_SAFE_CALL(cudaSetDevice(thread_info[pdata->node_info->nspecies+i]->gpu_info->igpu));
			gpu_particles[i].copy_from(cpu_particles+i+pdata->nspecies);
		}
	}

	if(pdata->node_info->nMIC > 0)
	{
//#pragma omp for
		for(int i=0;i<pdata->node_info->nMIC;i++)
		{

			mic_particles[i].copy_from(cpu_particles+1+i+pdata->node_info->nGPU);
		}
	}


}

int pthreads_return_val = 1;

void* pthreads_push_interface(void* data)
{

	pthreads_push_data* tdata = ((pthreads_push_data*)data);


	tdata->particles->pthreads_push(tdata->fields,tdata->moments,tdata->tid);



	pthread_exit((void*)&pthreads_return_val);

//	return (void*)&pthreads_return_val;
}

void NodeParticleList::pthreads_push(NodeFieldData* fields,NodeHOMoments* moments,int tid)
{
	int device_type = thread_info[tid]->device_type;
	int device_id = thread_info[tid]->device_id;

	printf("Pthreads pushing particle list %i\n",tid);


	if(device_type == 0)
	{
		cpu_push_time[device_id].reset();
		cpu_push_time[device_id].start();
		nsubcycles_pushed[tid] = (cpu_particles+device_id)->push(pdata,fields->cpu_fields,moments->cpu_moments);
		cpu_push_time[device_id].stop();

		printf("CPU Push took %e ms\n",cpu_push_time[device_id].get_cummulative());

	}
	else if(device_type == 1)
	{
		gpu_push_time[device_id].reset();
		gpu_push_time[device_id].start();
		CUDA_SAFE_CALL(cudaSetDevice(thread_info[tid]->gpu_info->igpu));

//		((fields->gpu_fields)+device_id)->copy_from(fields->cpu_fields);
		// Push the particles
		nsubcycles_pushed[tid] = (gpu_particles+device_id)->push(pdata,(fields->gpu_fields)+device_id,
				(moments->gpu_moments)+device_id);

		// Tramsfer the moments to the host
		(moments->gpu_moments+device_id)->copy_to(moments->cpu_moments+pdata->node_info->nCPU+device_id);
		gpu_push_time[device_id].stop();



		int gpu_device = thread_info[tid]->gpu_info->igpu;
//		printf("Pushing particles on GPU %i\n",gpu_device);

		printf("GPU %i Push took %e ms\n",thread_info[tid]->gpu_info->igpu,gpu_push_time[device_id].get_cummulative());
	}
	else if(device_type == 2)
	{
		nsubcycles_pushed[tid] = (mic_particles+device_id)->push(pdata,fields->mic_fields+device_id,moments->mic_moments);
		moments->cpu_moments[pdata->node_info->nCPU+
		pdata->node_info->nGPU+device_id].copy_from(moments->mic_moments+device_id);
	}
}

long long int NodeParticleList::push(NodeFieldData* fields,NodeHOMoments* moments)
{


	printf("Pushing Particles\n");
	int nLists = pdata->nspecies+pdata->node_info->nGPU+pdata->node_info->nMIC;
	pthread_t threads[nLists];
	pthreads_push_data* data = (pthreads_push_data*)malloc(nLists*sizeof(pthreads_push_data));

	int rc;
	void* status;

	pthread_attr_t attr;
   /* Initialize and set thread detached attribute */
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);


	omp_set_nested(0);
	omp_set_dynamic(1);
	for(int i=0;i<nLists;i++)
	{
		data[i].fields = fields;
		data[i].particles = this;
		data[i].moments = moments;
		data[i].tid = i;
		rc = pthread_create(&threads[i],&attr,pthreads_push_interface,&(data[i]));
	  if (rc) {
		 printf("ERROR; return code from pthread_create() is %d\n", rc);
		 exit(-1);
		 }

	}

	pthread_attr_destroy(&attr);

	int nDone = 0;
	while(nDone < nLists)
	{
		nDone = 0;
		for(int i=0;i<nLists;i++)
		{
			pthread_join(threads[i],&status);

			nDone += *((int*)status);
		}

	}


	printf("Finished pushing particles\n");
	free(data);
//	for(int i=0;i<nLists;i++)
//	{
//		pthreads_push(fields,moments,i);
//	}
//	pthread_exit(NULL);

	for(int i=1;i<nLists;i++)
		nsubcycles_pushed[0] += nsubcycles_pushed[i];



	return nsubcycles_pushed[0];
}


void NodeParticleList::copy_from(NodeParticleList* src)
{


	int nLists = pdata->nspecies+pdata->node_info->nGPU+pdata->node_info->nMIC;
	long long int result[nLists];

    omp_set_dynamic(0);
    omp_set_nested(1);
#pragma omp parallel num_threads(nLists)
    {
    	int tid = omp_get_thread_num();

    	int device_type = thread_info[tid]->device_type;
    	int device_id = thread_info[tid]->device_id;

    	if(device_type == 0)
    	{
    		cpu_particles[device_id].copy_from(src->cpu_particles+device_id);
    	}

    	if(device_type == 1)
    	{
    		CUDA_SAFE_CALL(cudaSetDevice(thread_info[tid]->gpu_info->igpu));
    		gpu_particles[device_id].copy_from(src->gpu_particles+device_id);
    	}

    	if(device_type == 2)
    	{
    		mic_particles[device_id].copy_from(src->mic_particles+device_id);
    	}



	}

}




















