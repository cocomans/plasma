#include "ParticleListGPUSorted.cuh"
#include "ParticleListCPU.h"
#include "CurrentTallyGPU.cuh"
#include "ChargeTallyGPU.h"
#include "StressTallyGPU.h"
#include "HOMomentsGPU.h"
#include "PlasmaData.h"
#include "ProblemInitializer.h"
#include "FieldDataGPU.h"
#include "math.h"
#include "omp.h"
#include "Util/GPU_Utils.h"
#include "ClusterInfo.cuh"
#include "PExitCheck.h"
#include "ParticleObjNT.h"
#include "FieldDataCPU.h"
#include "ParallelInfo.h"
#include "RunData.h"
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>



__constant__ PlasmaData pdata_c;
__constant__ HOMomentsGPU moments_c;

//FieldDataGPU* fields_gpu;
//__constant__ ParticleListGPU* particlesGPU_c;

const PlasmaData* pdata_c_address = &pdata_c;
const HOMomentsGPU* moments_c_address = &moments_c;
//const FieldDataGPU* fields_c_address = &fields_c;
//const ParticleListGPU* particlesGPU_c_address = &particlesGPU_c;

__host__ __device__ __attribute__((noinline))
ParticleListGPUSorted::ParticleListGPUSorted()
{
	device_type = 1;
}

__host__ __device__ __attribute__((noinline))
ParticleListGPUSorted::~ParticleListGPUSorted()
{
	/*
	// Free realkind arrays
	for(int i=0;i<ParticleList_nfloats;i++)
	{
		free(*get_float(i));
	}

	// Allocate int arrays
	for(int i=0;i<ParticleList_nints;i++)
	{
		free(*get_int(i));
	}

	// allocate short ints for cluster id's
	free(cluster_id);
	*/
}

void ParticleListGPUSorted::copy_from(const ParticleListGPUSorted* list_in)
{
	ispecies = list_in -> ispecies;

	enum cudaMemcpyKind kind;

	kind = cudaMemcpyDeviceToDevice;

	// Free realkind arrays
	for(int i=0;i<8;i++)
	{
		CUDA_SAFE_CALL(cudaMemcpyAsync(*get_float(i),*(list_in->get_float(i)),nptcls*sizeof(realkind),kind));
	}

	// Allocate int arrays
	for(int i=0;i<3;i++)
	{
		CUDA_SAFE_CALL(cudaMemcpyAsync(*get_int(i),*(list_in->get_int(i)),nptcls*sizeof(int),kind));
	}

	// allocate short ints for cluster id's
	CUDA_SAFE_CALL(cudaMemcpyAsync(cluster_id,(list_in->cluster_id),nptcls*sizeof(int),kind));
	CUDA_SAFE_CALL(cudaMemcpyAsync(num_subcycles,(list_in->num_subcycles),nptcls*sizeof(int),kind));
//	CUDA_SAFE_CALL(cudaMemcpyAsync(particleIDs,(list_in->particleIDs),nptcls*sizeof(int),kind));
//	CUDA_SAFE_CALL(cudaMemcpyAsync(particleIDs_original,(list_in->particleIDs_original),nptcls*sizeof(int),kind));

	CUDA_SAFE_CALL(cudaDeviceSynchronize());

}


void ParticleListGPUSorted::copy_from(const ParticleListCPU* list_in)
{
	ispecies = list_in -> ispecies;

	enum cudaMemcpyKind kind;

	kind = cudaMemcpyHostToDevice;


	// Free realkind arrays
	for(int i=0;i<7;i++)
	{
		CUDA_SAFE_CALL(cudaMemcpyAsync(*get_float(i),*(list_in->get_float(i)),nptcls*sizeof(realkind),kind));
	}

	// Allocate int arrays
	for(int i=0;i<3;i++)
	{
		CUDA_SAFE_CALL(cudaMemcpyAsync(*get_int(i),*(list_in->get_int(i)),nptcls*sizeof(int),kind));
	}

	// allocate short ints for cluster id's
//	CUDA_SAFE_CALL(cudaMemcpyAsync(cluster_id,(list_in->cluster_id),nptcls*sizeof(int),kind));
//	CUDA_SAFE_CALL(cudaMemcpyAsync(num_subcycles,(list_in->num_subcycles),nptcls*sizeof(int),kind));
//	CUDA_SAFE_CALL(cudaMemcpyAsync(particleIDs,(list_in->particleIDs),nptcls*sizeof(int),kind));
//	CUDA_SAFE_CALL(cudaMemcpyAsync(particleIDs_original,(list_in->particleIDs_original),nptcls*sizeof(int),kind));

	CUDA_SAFE_CALL(cudaDeviceSynchronize());

}

void ParticleListGPUSorted::copy_to(ParticleListCPU* list_in)
{
	list_in -> ispecies = ispecies;

	enum cudaMemcpyKind kind;

	kind = cudaMemcpyDeviceToHost;


	// Free realkind arrays
	for(int i=0;i<8;i++)
	{
		CUDA_SAFE_CALL(cudaMemcpyAsync(*(list_in->get_float(i)),*get_float(i),nptcls*sizeof(realkind),kind));
	}

	// Allocate int arrays
	for(int i=0;i<3;i++)
	{
		CUDA_SAFE_CALL(cudaMemcpyAsync(*(list_in->get_int(i)),*get_int(i),nptcls*sizeof(int),kind));
	}

	// allocate short ints for cluster id's
	CUDA_SAFE_CALL(cudaMemcpyAsync((list_in->cluster_id),cluster_id,nptcls*sizeof(int),kind));
	CUDA_SAFE_CALL(cudaMemcpyAsync((list_in->num_subcycles),num_subcycles,nptcls*sizeof(int),kind));
//	CUDA_SAFE_CALL(cudaMemcpyAsync(particleIDs,(list_in->particleIDs),nptcls*sizeof(int),kind));
//	CUDA_SAFE_CALL(cudaMemcpyAsync(particleIDs_original,(list_in->particleIDs_original),nptcls*sizeof(int),kind));

	CUDA_SAFE_CALL(cudaDeviceSynchronize());

}

void ParticleListGPUSorted::allocate(PlasmaData* pdata,int nptcls_in)
{
	printf("Allocating a GPU Particle List\n");

	ClusterSortDim = pdata->node_info->gpu_info->ClusterSortDim;
	device_type = 1;
	// Allocate memory for particles
	nptcls_allocated = nptcls_in;

	nptcls = nptcls_in;

	//plot = gnuplot_init();

	// Allocate realkind arrays
	for(int i=0;i<9;i++)
	{
		CUDA_SAFE_CALL(cudaMalloc((void**)get_float(i),nptcls_allocated*sizeof(realkind)));
	}

	// Allocate int arrays
	for(int i=0;i<5;i++)
	{
		CUDA_SAFE_CALL(cudaMalloc((void**)get_int(i),nptcls_allocated*sizeof(int)));
	}

	// allocate short ints for cluster id's
	CUDA_SAFE_CALL(cudaMalloc((void**)&cluster_id,nptcls_allocated*sizeof(int)));


	CUDA_SAFE_CALL(cudaMalloc((void**)&particleIDs,nptcls_allocated*sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&particleIDs_original,nptcls_allocated*sizeof(int)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&buffer64,nptcls_allocated*sizeof(realkind)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&buffer32,nptcls_allocated*sizeof(int)));




	gridsize = 6*16;
	blocksize = 128;

	gridsize = min(gridsize,(nptcls+blocksize-1)/blocksize);

	CUDA_SAFE_CALL(cudaMalloc((void**)&nsubcycles_thread,gridsize*blocksize*sizeof(int)));

	int nclusters_x = (pdata->nx+ClusterSortDim-1)/ClusterSortDim;
	int nclusters_y = (pdata->ny+ClusterSortDim-1)/ClusterSortDim;
	int nclusters_z = (pdata->nz+ClusterSortDim-1)/ClusterSortDim;

	nclusters = nclusters_x*nclusters_y*nclusters_z;
	printf("nclusters = %i\n",nclusters);

	CUDA_SAFE_CALL(cudaMalloc((void**)&clusters,nclusters*sizeof(ClusterInfo)));



	CUDA_SAFE_CALL(cudaGetSymbolAddress((void**)&pdata_d,pdata_c));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(pdata_c,pdata,sizeof(PlasmaData)));

}

void ParticleListGPUSorted::allocate(PlasmaData* pdata,ParticleListGPUSorted* _old)
{
	printf("Allocating a GPU Particle List\n");
	*this = *_old;
	//plot = gnuplot_init();

	// Allocate realkind arrays
	for(int i=0;i<7;i++)
	{
		CUDA_SAFE_CALL(cudaMalloc((void**)get_float(i),nptcls_allocated*sizeof(realkind)));
	}

	// Allocate int arrays
	for(int i=0;i<3;i++)
	{
		CUDA_SAFE_CALL(cudaMalloc((void**)get_int(i),nptcls_allocated*sizeof(int)));
	}

	// allocate short ints for cluster id's
	CUDA_SAFE_CALL(cudaMalloc((void**)&cluster_id,nptcls_allocated*sizeof(int)));


	CUDA_SAFE_CALL(cudaMalloc((void**)&particleIDs,nptcls_allocated*sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&particleIDs_original,nptcls_allocated*sizeof(int)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&buffer64,nptcls_allocated*sizeof(realkind)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&buffer32,nptcls_allocated*sizeof(int)));
	buffer = buffer64;

}

void ParticleListGPUSorted::init(ProblemInitializer* initializer, HOMomentsGPU* moments)
{

	CUDA_SAFE_CALL(cudaMemset(dt_finished,0,nptcls*sizeof(realkind)));




	if(moments->pdata->ndimensions == 2)
	{

	int blockSize = 512;
	int gridSize = 96;
	printf("Writing particle cluster IDs\n");
	CUDA_SAFE_KERNEL((write_cluster_ids<2,3><<<blockSize,gridSize>>>(
			pdata_d,*this,ifinished,nptcls)));


	//check_int_vals(cluster_id,nptcls);

	printf("Calling Thrust sort\n");
	SortByKey(particleIDs,cluster_id,nptcls);

	//check_int_vals(cluster_id,nptcls);
	ReorderData(particleIDs,nptcls);

	//check_int_vals(cluster_id,nptcls);

	}

}

realkind ParticleListGPUSorted::evaluate_energy(PlasmaData* pdata)
{
	double etotal = 0.0;
	for(int i=0;i<nptcls;i++)
	{
//		etotal += get_fvalue(i,3)* get_fvalue(i,3);
//		etotal += get_fvalue(i,4)* get_fvalue(i,4);
//		etotal += get_fvalue(i,5)* get_fvalue(i,5);
	}

	etotal = etotal * pdata->mspecies[ispecies] * 0.5/((double)pdata->nptcls);

	return etotal;
}

double4 ParticleListGPUSorted::subcycle_stats(PlasmaData* pdata)
{
	int* num_subcycles_temp = (int*)malloc(nptcls*sizeof(int));

	CUDA_SAFE_CALL(cudaMemcpy(num_subcycles_temp,num_subcycles,nptcls*sizeof(int),
			cudaMemcpyDeviceToHost));

	double scale = pdata->rstatus->npiccard_outer;
	double mean = 0;
	double mean2 = 0;
	double mins = num_subcycles_temp[0]/scale;
	double maxs = num_subcycles_temp[0]/scale;
	int imax = 0;
	int imin = 0;
	for(int i=0;i<nptcls;i++)
	{
		if(mins > num_subcycles_temp[i]/scale)
		{
			mins = num_subcycles_temp[i]/scale;
			imin = i;
		}
		if(maxs < num_subcycles_temp[i]/scale)
		{
			maxs = num_subcycles_temp[i]/scale;
			imax = i;
		}

		mean += num_subcycles_temp[i]/((double)nptcls*scale);
		mean2 += num_subcycles_temp[i]*num_subcycles_temp[i]/((double)nptcls*scale*scale);
	}

	double std_diff;

	std_diff = sqrt(fabs(mean*mean - mean2));
	printf("Particle Subcycle Stats:\n");
	printf("Avg Subcycles: %f +/- %f\n",mean,std_diff);
	printf("Min / Max: %f[%i] / %f[%i]\n",mins,imin,maxs,imax);

	free(num_subcycles_temp);

	return make_double4(mean,std_diff,mins,maxs);
}

double4 ParticleListGPUSorted::piccard_stats(PlasmaData* pdata)
{
	int* num_subcycles_temp = (int*)malloc(nptcls*sizeof(int));
	double* num_piccard_temp = (double*)malloc(nptcls*sizeof(double));
//	double* num_piccard2_temp = (double*)malloc(nptcls*sizeof(double));

	CUDA_SAFE_CALL(cudaMemcpy(num_subcycles_temp,num_subcycles,nptcls*sizeof(int),
			cudaMemcpyDeviceToHost));

	CUDA_SAFE_CALL(cudaMemcpy(num_piccard_temp,num_piccard,nptcls*sizeof(double),
			cudaMemcpyDeviceToHost));

//	CUDA_SAFE_CALL(cudaMemcpy(num_piccard2_temp,num_piccard2,nptcls*sizeof(double),
//			cudaMemcpyDeviceToHost));

	double scale = pdata->rstatus->npiccard_outer;
	double mean = 0;
	double mean2 = 0;
	double mins = num_piccard_temp[0]/num_subcycles_temp[0];
	double maxs = num_piccard_temp[0]/num_subcycles_temp[0];
	int imax = 0;
	int imin = 0;
	for(int i=0;i<nptcls;i++)
	{
		if(mins > num_piccard_temp[i]/num_subcycles_temp[i])
		{
			mins = num_piccard_temp[i]/num_subcycles_temp[i];
			imin = i;
		}
		if(maxs < num_piccard_temp[i]/num_subcycles_temp[i])
		{
			maxs = num_piccard_temp[i]/num_subcycles_temp[i];
			imax = i;
		}

		mean += num_piccard_temp[i]/((double)nptcls*num_subcycles_temp[i]);
		mean2 += num_piccard_temp[i]*num_piccard_temp[i]/((double)nptcls*num_subcycles_temp[i]*num_subcycles_temp[i]);
	}

	double std_diff;

	std_diff = sqrt(fabs(mean*mean - mean2));
	printf("Particle Piccard Stats(GPU):\n");
	printf("Avg Piccard: %f +/- %f\n",mean,std_diff);
	printf("Min / Max: %f[%i] / %f[%i]\n",mins,imin,maxs,imax);

	free(num_subcycles_temp);
	free(num_piccard_temp);


	return make_double4(mean,std_diff,mins,maxs);
}



void ParticleListGPUSorted::plot_particles(PlasmaData* pdata)
{
//	float* x_vals = (float*)malloc(nptcls*sizeof(float));
//	float* y_vals = (float*)malloc(nptcls*sizeof(float));
//
//	realkind* px_t = (realkind*)malloc(nptcls*sizeof(realkind));
//	int* ix_t = (int*)malloc(nptcls*sizeof(int));
//	realkind* vx_t = (realkind*)malloc(nptcls*sizeof(realkind));
//
//	CUDA_SAFE_CALL(cudaMemcpy(px_t,px,nptcls*sizeof(realkind),cudaMemcpyDeviceToHost));
//	CUDA_SAFE_CALL(cudaMemcpy(ix_t,ix,nptcls*sizeof(int),cudaMemcpyDeviceToHost));
//	CUDA_SAFE_CALL(cudaMemcpy(vx_t,vx,nptcls*sizeof(realkind),cudaMemcpyDeviceToHost));
//
//	for(int i=0;i<nptcls;i++)
//	{
//		x_vals[i] = (px_t[i]+ix_t[i])*pdata->dxdi + pdata->xmin;
//		y_vals[i] = vx_t[i];
//
////		printf("particle[%i]: %f %f\n",i,x_vals[i],y_vals[i]);
//	}
//
//	gnuplot_resetplot(plot);
//
//	gnuplot_plot_xy(plot,x_vals,y_vals,nptcls,NULL);
//
//
//	free(x_vals);
//	free(y_vals);
//	free(px_t);
//	free(vx_t);
//	free(ix_t);
}

void ParticleListGPUSorted::CPUfree()
{
//	// Allocate realkind arrays
//	for(int i=0;i<ParticleList_nfloats;i++)
//	{
//		CUDA_SAFE_CALL(cudaFree(*get_float(i)));
//	}
//
//	// Allocate int arrays
//	for(int i=0;i<ParticleList_nints;i++)
//	{
//		CUDA_SAFE_CALL(cudaFree(*get_int(i)));
//	}
//
//	// allocate short ints for cluster id's
//	CUDA_SAFE_CALL(cudaFree(cluster_id));
}

long long int ParticleListGPUSorted::push(PlasmaData* pdata,
		FieldDataGPU* fields, HOMomentsGPU* moments)
{
	long long int result;

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(pdata_c,pdata,sizeof(PlasmaData)));
	// Copy field data to the device

//	printf("Nulling device HOMoments\n");
	moments -> set_vals(0);

	// Change location of pdata in moments_d to device
	moments->pdata = pdata_d;
	fields->pdata = pdata_d;

//	FieldCheckGPU(*fields);
	// Push the particles
//	printf("Pushing particles on the GPU\n");
	//result = push_interface2(pdata,fields,moments);
	result = push_interface(pdata,fields,moments);
	printf("nsubcycles GPU = %i\n",result);
	// Tally Charge and S2 Moments


//	printf("More Pushing particles on the GPU\n");
	// Change the location of pdata in moments_d to host
	moments->pdata = pdata;
	fields->pdata = pdata;
	// Copy HO moments to the host
//	printf("Copying device HOMoments to Host\n");
//	moments->copy_from(moments_d);

	return result;
}

template<int nSpatial,int nVel,bool iEM> __global__
void SimpleGPUPush(PlasmaData* 				pdata,
					FieldDataGPU 			fields,
					HOMomentsGPU 			moments,
					ParticleListGPUSorted	particles,
					PExitCheckCPU			exit_checker,
					int* num_subcycles)
{
	int tidx = threadIdx.x;
	int bidx = blockIdx.x;
	int gidx = tidx+blockDim.x*bidx;
	int stride = blockDim.x*gridDim.x;

	int pid = gidx;

	long long int num_subcycles_thread = 0;

	__shared__ float currentx[1024];
	__shared__ float charge_s[1024];
	__shared__ float S2xx_s[1024];
	while(tidx < 1024)
	{
		currentx[tidx] = 0;
		charge_s[tidx] = 0;
		S2xx_s[tidx] = 0;
		tidx += blockDim.x;
	}
	tidx = threadIdx.x;



//	PlasmaData pdata = *pdata2;

	ParticleObjNT<1,nSpatial,nVel,iEM,
	ParticleListGPUSorted,FieldDataGPU,CurrentTallyGPU,
	PExitCheckCPU,PExitCheckCPU> particle(&pid,&exit_checker);
	typevecN<int,1> iter;

	CurrentTallyGPU currents((float*)currentx, (float*)&moments.get_val(0,0,0,particles.ispecies,HOMoments_currenty),
						  (float*)&moments.get_val(0,0,0,particles.ispecies,HOMoments_currentz),
						  pdata->nx,pdata->ny,pdata->nz,
						  0,0,0,
						  pdata->ndimensions);

	ChargeTallyGPU charge((float*)charge_s,
						  make_int3(currents.nx,currents.ny,currents.nz),
						  pdata->dxdi,pdata->dydi,pdata->dzdi,
						  1);

	StressTallyGPU stress((float*)S2xx_s,
			  make_int3(currents.nx,currents.ny,currents.nz),
			  pdata->dxdi,pdata->dydi,pdata->dzdi,
			  1,1);

	particle.species = particles.ispecies;


	while(pid < particles.nptcls)
	{
		iter(0) = 0;
		//printf("reading in particle %i\n",pid);
		particle.copy_in_gpu(particles,0);
		particle.dt_finished(0) = 0;
		particle.push(pdata,&fields,&currents,iter,pdata->nSubcycle_max);

//		printf("Writing Paerticles Back\n");
		particle.write_back(particles,0);

		particles.num_subcycles[pid] += iter(0);

		charge.tally1d(particles.px[pid],particles.ix[pid],1.0);
		stress.tally1d1v(particles.px[pid],particles.vx[pid],particles.ix[pid],1.0);

		num_subcycles_thread += iter(0);

		pid += stride;
	}

//	printf("exiting gpu shit\n");

	num_subcycles[gidx] = num_subcycles_thread;

	__syncthreads();

	pid = tidx;
	while(pid < pdata->nx)
	{
#ifdef DOUBLE_PRECISION
//		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currentx),((double*)currentx)[pid]);
//		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_charge),((double*)charge_s)[pid]);
//		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xx),((double*)S2xx_s)[pid]);

		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currentx),((float*)currentx)[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_charge),((float*)charge_s)[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xx),((float*)S2xx_s)[pid]);
#else
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currentx),currentx[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_charge),charge_s[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xx),S2xx_s[pid]);
#endif
		pid += blockDim.x;
	}

}

template<int nSpatial,int nVel,bool iEM> __global__
void SimpleGPUPush1D3V(PlasmaData* 				pdata,
					FieldDataGPU 			fields,
					HOMomentsGPU 			moments,
					ParticleListGPUSorted	particles,
					PExitCheckCPU			exit_checker,
					int* num_subcycles)
{
	int tidx = threadIdx.x;
	int bidx = blockIdx.x;
	int gidx = tidx+blockDim.x*bidx;
	int stride = blockDim.x*gridDim.x;

	int pid = gidx;

	long long int num_subcycles_thread = 0;

	__shared__ float currentx[256];
	__shared__ float currenty[256];
	__shared__ float currentz[256];

	__shared__ float charge_s[256];

	__shared__ float S2xx_s[256];
	__shared__ float S2xy_s[256];
	__shared__ float S2xz_s[256];
	__shared__ float S2yy_s[256];
	__shared__ float S2yz_s[256];
	__shared__ float S2zz_s[256];




	while(tidx < 256)
	{
		currentx[tidx] = 0;
		currenty[tidx] = 0;
		currentz[tidx] = 0;

		charge_s[tidx] = 0;
		S2xx_s[tidx] = 0;
		S2xy_s[tidx] = 0;
		S2xz_s[tidx] = 0;
		S2yy_s[tidx] = 0;
		S2yz_s[tidx] = 0;
		S2zz_s[tidx] = 0;

		tidx += blockDim.x;
	}
	tidx = threadIdx.x;



//	PlasmaData pdata = *pdata2;

	ParticleObjNT<1,nSpatial,nVel,iEM,
	ParticleListGPUSorted,FieldDataGPU,CurrentTallyGPU,
	PExitCheckCPU,PExitCheckCPU> particle(&pid,&exit_checker);
	typevecN<int,1> iter;

	CurrentTallyGPU currents((float*)currentx, (float*)currenty,
						  (float*)currentz,
						  pdata->nx,pdata->ny,pdata->nz,
						  0,0,0,
						  pdata->ndimensions);

	ChargeTallyGPU charge((float*)charge_s,
						  make_int3(currents.nx,currents.ny,currents.nz),
						  pdata->dxdi,pdata->dydi,pdata->dzdi,
						  1);

	StressTallyGPU stress((float*)S2xx_s,(float*)S2xy_s,(float*)S2xz_s,(float*)S2yy_s,(float*)S2yz_s,(float*)S2zz_s,
			  currents.nx,currents.ny,currents.nz,
			  nSpatial,nVel);

	particle.species = particles.ispecies;


	while(pid < particles.nptcls)
	{
		iter(0) = 0;
		//printf("reading in particle %i\n",pid);
		particle.copy_in_gpu(particles,0);
		particle.dt_finished(0) = 0;
		particle.push(pdata,&fields,&currents,iter,pdata->nSubcycle_max);

//		printf("Writing Paerticles Back\n");
		particle.write_back(particles,0);

		particles.num_subcycles[pid] += iter(0);

		charge.tally1d(particles.px[pid],particles.ix[pid],1.0);
		stress.tally1d3v(particles.px[pid],particles.vx[pid],particles.vy[pid],particles.vz[pid],particles.ix[pid],1.0);

		num_subcycles_thread += iter(0);

		pid += stride;
	}

//	printf("exiting gpu shit\n");

	num_subcycles[gidx] = num_subcycles_thread;

	__syncthreads();

	pid = tidx;
	while(pid < pdata->nx)
	{
#ifdef DOUBLE_PRECISION
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currentx),currentx[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currenty),currenty[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currentz),currentz[pid]);

		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_charge),charge_s[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xx),S2xx_s[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xy),S2xy_s[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xz),S2xz_s[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2yy),S2yy_s[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2yz),S2yz_s[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2zz),S2zz_s[pid]);


#else
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currentx),currentx[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currenty),currenty[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currentz),currentz[pid]);

		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_charge),charge_s[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xx),S2xx_s[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xy),S2xy_s[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xz),S2xz_s[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2yy),S2yy_s[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2yz),S2yz_s[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2zz),S2zz_s[pid]);
#endif
		pid += blockDim.x;
	}

}

template<int nSpatial,int nVel,bool iEM>
void SimpleGPUPushH(PlasmaData* 			pdata,
					FieldDataGPU* 			fields,
					HOMomentsGPU* 				moments,
					ParticleListGPUSorted*	particles,
					int* num_subcycles)
{
	int blocksize = particles->blocksize;
	int gridsize = particles->gridsize;
	PExitCheckCPU exit_checker(pdata->dt,pdata->nSubcycle_max);
	printf("GPU Pushing Particles\n");
	if((nSpatial == 1)&&(nVel == 1)){
	CUDA_SAFE_KERNEL((SimpleGPUPush<nSpatial,nVel,iEM><<<blocksize,gridsize>>>
			(particles->pdata_d,*fields,*moments,*particles,exit_checker,num_subcycles)));}
	else if((nSpatial == 1)&&(nVel == 3)){
	CUDA_SAFE_KERNEL((SimpleGPUPush1D3V<nSpatial,nVel,iEM><<<blocksize,gridsize>>>
			(particles->pdata_d,*fields,*moments,*particles,exit_checker,num_subcycles)));}
	//printf("GPU Finished Pushing Particles\n");
}


long long int ParticleListGPUSorted::push_interface(PlasmaData* pdata,
			FieldDataGPU* fields,
			HOMomentsGPU* moments)
{

//		FieldCheckGPU(*fields);

//		getchar();

	int iDevice;
	CUDA_SAFE_CALL(cudaGetDevice(&iDevice));
		printf("GPU %i Pushing particles of species %i\n",iDevice,ispecies);
		if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1)&&(pdata->iEM == 0))
			SimpleGPUPushH<1,1,0>(pdata,fields,moments,this,nsubcycles_thread);
		else if((pdata->ndimensions == 1)&&(pdata->nVelocity == 3)&&(pdata->iEM == 1))
			SimpleGPUPushH<1,3,1>(pdata,fields,moments,this,nsubcycles_thread);



		thrust::device_ptr<int> nsubcycles_t(nsubcycles_thread);


		// Reduce and return the total number of particle-subcycles taken
		long long int result = thrust::reduce(nsubcycles_t,nsubcycles_t+blocksize*gridsize);

//		int* nsubcycles_temp = (int*)malloc(gridsize*blocksize*sizeof(int));
//		CUDA_SAFE_CALL(cudaMemcpy(nsubcycles_temp,nsubcycles_thread,gridsize*blocksize*sizeof(int),cudaMemcpyDeviceToHost))
//
//		for(int i=0;i<gridsize*blocksize;i++)
//			result += nsubcycles_temp[i];

//		MomentCheckGPU(*moments);

		return result;


}


template<int nSpatial,int nVel> __global__
void write_cluster_ids(PlasmaData* pdata,
					ParticleListGPUSorted particles,
					int* ifinished,
					int nptcls_check)
{
//	int tidx = threadIdx.x;
//	int bidx = blockIdx.x;
//	int gidx = tidx+blockDim.x*bidx;
//	int stride = blockDim.x*gridDim.x;
//
//	int nclusters_x = (pdata->nx+pdata->ClusterSortDim-1)/pdata->ClusterSortDim;
//	int nclusters_y = (pdata->ny+pdata->ClusterSortDim-1)/pdata->ClusterSortDim;
//	int nclusters_z = (pdata->nz+pdata->ClusterSortDim-1)/pdata->ClusterSortDim;
//
//	while(gidx < nptcls_check)
//	{
//		int ix,iy,iz;
//		ix = particles.ix[gidx]/pdata->ClusterSortDim;
//		iy = particles.iy[gidx]/pdata->ClusterSortDim;
//		iz = particles.iz[gidx]/pdata->ClusterSortDim;
//
//		particles.cluster_id[gidx] = ix + nclusters_x*(iy + nclusters_y*iz);
//		particles.particleIDs[gidx] = gidx;
//		printf("ix,iy,iz[%i] = %i %i %i %i\n",gidx,ix,iy,iz,particles.cluster_id[gidx]);
//
//		gidx += stride;
//	}

}


template<int nSpatial,int nVel,bool iEM> __global__
void GPUBlockPush(PlasmaData* 				pdata,
					FieldDataGPU* 			fields,
					HOMomentsGPU* 				moments,
					ParticleListGPUSorted*	particles,
					int* num_subcycles,
					int tally_size,
					GPUInfo* gpu_info)
{
//	int tidx = threadIdx.x;
//	int bidx = blockIdx.x;
//	int gidx = tidx+blockDim.x*bidx;
//	int stride = blockDim.x*gridDim.x;
//
//
//	int pid;
//
//	long long int num_subcycles_thread = 0;
//
//	// Going to shoot for about 4 thread-blocks per cluster,
//	int ClustID = bidx/gpu_info->TBpCluster;
//	int ClustID2 = bidx%gpu_info->TBpCluster;
//
//	if(ClustID > particles->nclusters-1)
//		return;
//	// Check to see if any particles are in the cluster
//	if((particles->clusters[ClustID].ifirst < 0)||(particles->clusters[ClustID].ilast < 0))
//		return;
//
//	__syncthreads();
//	// Allocate shared memory for current tallies
//	__shared__ float currentx[33*33];
//	__shared__ float currenty[33*33];
//	__shared__ float currentz[33*33];
//
//
//
//
//	// Shared copy of the cluster info
//	__shared__ ClusterInfo cluster[1];
//
//
//	// Null out the tallies
//	while(tidx < tally_size)
//	{
//		currentx[tidx] = 0;
//		currenty[tidx] = 0;
//		currentz[tidx] = 0;
//		tidx += blockDim.x;
//	}
//
//	tidx = threadIdx.x;
//	__syncthreads();
//
//	if(tidx == 0)
//		*cluster = particles->clusters[ClustID];
//
//	__syncthreads();
//
//	int ix0c,iy0c,iz0c;
//	int ioffset;
//	int ix0,iy0,iz0;
//	int nx,ny;
//	// First figure out the current tally domain
//	ioffset = floor((gpu_info->ClusterStorDim - gpu_info->ClusterSortDim)/2.0f);
//	if(nSpatial == 1)
//	{
//		ix0c = cluster->clusterid*gpu_info->ClusterSortDim;
//		iy0c = 0;
//		iz0c = 0;
//
//		ix0 = ix0c - ioffset;
//		iy0 = 0;
//		iz0 = 0;
//	}
//	else if(nSpatial == 2)
//	{
//		int nclustx = (pdata->nx+gpu_info->ClusterSortDim-1)/gpu_info->ClusterSortDim;
//		iy0c = cluster->clusterid/nclustx;
//		ix0c = cluster->clusterid - iy0c*nclustx;
//		iy0c *= gpu_info->ClusterSortDim;
//		ix0c *= gpu_info->ClusterSortDim;
//
//		ix0 = ix0c - ioffset;
//		iy0 = iy0c - ioffset;
//		iz0 = 0;
//	}
//	else if(nSpatial == 3)
//	{
//		int nclusty = (pdata->ny+gpu_info->ClusterSortDim-1)/gpu_info->ClusterSortDim;
//		int nclustx = (pdata->nx+gpu_info->ClusterSortDim-1)/gpu_info->ClusterSortDim;
//		iz0c = cluster->clusterid/(nclusty*nclustx);
//		iy0c = (cluster->clusterid - iz0c*nclustx*nclusty)/nclustx;
//		ix0c = cluster->clusterid - (iy0c+iz0c*nclusty)*nclustx;
//		iz0c *= gpu_info->ClusterSortDim;
//		iy0c *= gpu_info->ClusterSortDim;
//		ix0c *= gpu_info->ClusterSortDim;
//
//		ix0 = ix0c - ioffset;
//		iy0 = iy0c - ioffset;
//		iz0 = iz0c - ioffset;
//	}
//	CurrentTallyGPU* currents = new CurrentTallyGPU(currentx,currenty,currentz,
//			gpu_info->ClusterStorDim,gpu_info->ClusterStorDim,gpu_info->ClusterStorDim,
//			ix0,iy0,iz0,pdata->ndimensions);
//	// Setup the PExitCheck object
//
//
//	// Load local CurrentTally Object to shared memory
////		*currents = currents_temp;
//
//	// Now time for the exit condition
//	// This guy needs to be 4 smaller than the tally region in order to avoid
//	// particles writing outside of the tally region
//	ix0 += 2;
//	iy0 += 2;
//	PExitCheck* exit_checker = new PExitCheckGPU(pdata->dt,pdata->nSubcycle_max,
//			ix0,iy0,iz0,
//			gpu_info->ClusterStorDim-4,gpu_info->ClusterStorDim-4,gpu_info->ClusterStorDim-4);
//
//
//
//
//	__syncthreads();
//
//	// Each thread figures out its own PID
//	pid = cluster->ifirst + ClustID2*blockDim.x + threadIdx.x;
////	printf("First particle in block %i is %i %i\n",bidx,cluster->ifirst+ClustID2*blockDim.x,cluster->ilast);
//
//	// Setup the ParticleObjectNT
//	ParticleObjNT<1,nSpatial,nVel,iEM> particle(&pid,exit_checker);
//
//
//	// Each thread pushes a bunch of particles
//	while(pid <= cluster->ilast)
//	{
//		typevecN<int,1> iter;
//		// Load a particle into a ParticleObjNT
//		particle.copy_in_gpu(particles,0);
//
//		// Set the subcycle count
//		iter(0) = particles->nsubcycles_current[pid];
//
//		// Push that particle
//		particle.push(pdata,fields,currents,iter,pdata->nSubcycle_max);
//
//		// ParticleObjNT has exited.
//		// Write ParticleObjNT data back to the Particle List
//		particle.write_back(particles,0);
//
//		// Save number of subcycles completed
//		// for the particle that was just pushed
//		particles.nsubcycles_current[pid] = iter(0);
//
//		// Mark particle finished or not for stream compaction array
//		if((particle.dt_finished(0) >= pdata->dt)||(iter(0) >= pdata->nSubcycle_max))
//			particles.ifinished[pid] = 1;
//
//		printf("Dt finished and iter done[%i] = %f %i\n",pid,particle.dt_finished(0),
//				iter(0));
//
//
//		// increment the number of subcycles this thread has processed
//		num_subcycles_thread += iter(0);
//
//		// increment pid and repeat
//		pid += 4*blockDim.x;
//	}
//
//	// synchronize threads
//	__syncthreads();
//
//
//	// Accumulate currents to global memory
//	pid = tidx;
//	while(pid < tally_size)
//	{
//		int ix,iy,iz;
//
//		if(nSpatial == 1)
//		{
//			ix = pid + currents->ix0;
//			iy = 0;
//			iz = 0;
//		}
//		else if(nSpatial == 2)
//		{
//			iy = pid/currents->nx;
//			ix = pid - iy*currents->nx;
//			iz = 0;
//
//			ix += currents->ix0;
//			iy += currents->iy0;
//		}
//		else if(nSpatial == 3)
//		{
//			iz = pid/(currents->nx*currents->ny);
//			iy = (pid - currents->nx*currents->ny*iz)/currents->nx;
//			ix = pid - currents->nx*(iy+currents->ny*iz);
//
//			ix += currents->ix0;
//			iy += currents->iy0;
//			iz += currents->iz0;
//		}
//
//#ifdef DOUBLE_PRECISION
//		atomicAddD(&moments->get_val(ix,iy,iz,particles.ispecies,HOMoments_currentx),currentx[pid]);
//		atomicAddD(&moments->get_val(ix,iy,iz,particles.ispecies,HOMoments_currenty),currenty[pid]);
//		atomicAddD(&moments->get_val(ix,iy,iz,particles.ispecies,HOMoments_currentz),currentz[pid]);
//#else
//		atomicAdd(&moments->get_val(ix,iy,iz,particles->ispecies,HOMoments_currentx),currentx[pid]);
//		atomicAdd(&moments->get_val(ix,iy,iz,particles->ispecies,HOMoments_currenty),currenty[pid]);
//		atomicAdd(&moments->get_val(ix,iy,iz,particles->ispecies,HOMoments_currentz),currentz[pid]);
//#endif
//		pid += blockDim.x;
//	}
}

__global__
void find_cluster_boundaries(ParticleListGPUSorted particles,ClusterInfo* bins,int nptcls_check)
{
//	int idx = threadIdx.x;
//	int gidx = idx+blockIdx.x*blockDim.x;
//
//	uint binindex;
//	uint binindex_left;
//	uint binindex_right;
//
//	while(gidx < nptcls_check)
//	{
//		if(gidx == 0)
//		{
//			binindex = particles.cluster_id[gidx];
//			bins[binindex].ifirst = gidx;
//			bins[binindex].clusterid = binindex;
//		}
//		else if(gidx == nptcls_check-1)
//		{
//			binindex = particles.cluster_id[gidx];
//			bins[binindex].ilast = gidx;
//			bins[binindex].clusterid = binindex;
//		}
//
//			binindex = particles.cluster_id[gidx];
//			binindex_left = particles.cluster_id[max(gidx-1,0)];
//			binindex_right = particles.cluster_id[min((gidx+1),(nptcls_check-1))];
//
//			if(binindex_left != binindex)
//			{
//				printf("ClusterID[%i] = %i\n",gidx,binindex);
//				bins[binindex].ifirst = gidx;
//				bins[binindex].clusterid = binindex;
//			}
//
//			if(binindex_right != binindex)
//			{
//				printf("ClusterID[%i] = %i\n",gidx,binindex);
//				bins[binindex].ilast = gidx;
//				bins[binindex].clusterid = binindex;
//			}
//
//
//
//
//
//
//		gidx += blockDim.x*gridDim.x;
//	}
}

__global__
void InitClusters(ClusterInfo* bins,
		int nbins)
{
//	int gidx = threadIdx.x+blockDim.x*blockIdx.x;
//
//	// Need to set ifirst, ilast to -1, so that once everything
//	// is populated, the bins with no particles have -1 as ifirst and ilast
//
//	while(gidx < nbins)
//	{
//		bins[gidx].ifirst = -1;
//		bins[gidx].ilast = -1;
//
//
//
//		gidx += blockDim.x*gridDim.x;
//	}
//

}

__global__
void FieldCheck_g(FieldDataGPU fields)
{
	int gidx = threadIdx.x+blockDim.x*blockIdx.x;
	int gidy = threadIdx.y+blockDim.y*blockIdx.y;
//	gidy = 0;




	while(gidy < 1)
	{
		while(gidx < fields.nx)
		{
			realkind Bx,By,Bz,Ex,Ey,Ez;

//			Bx = fields.getB(gidx,gidy,0,0);
//			By = fields.getB(gidx,gidy,0,1);
//			Bz = fields.getB(gidx,gidy,0,2);
			Ex = fields.getE(gidx,gidy,0,0);
//			Ey = fields.getE(gidx,gidy,0,1);
//			Ez = fields.getE(gidx,gidy,0,2);

//			Ex = fields.intrpET<0,FieldData_deriv_f>(0.0,0,0,gidx,0,0);

//			printf("GPU Field[%i %i] = %e, %e, %e, %e, %e, %e\n",gidx,gidy,Bx,By,Bz,Ex,Ey,Ez);

			printf("GPU Field[%i %i] = %e\n",gidx,gidy,Ex);

			gidx += blockDim.x*gridDim.x;
		}

		gidy += blockDim.y*gridDim.y;
	}

}




__global__
void MomentCheck_g(HOMomentsGPU moments)
{
	int gidx = threadIdx.x+blockDim.x*blockIdx.x;
	int gidy = threadIdx.y+blockDim.y*blockIdx.y;
//	gidy = 0;




	while(gidy < 1)
	{
		while(gidx < moments.nx)
		{
			realkind Bx,By,Bz,Jx,Ey,Ez;

//			Bx = fields.getB(gidx,gidy,0,0);
//			By = fields.getB(gidx,gidy,0,1);
//			Bz = fields.getB(gidx,gidy,0,2);
			Jx = moments.get_val(gidx,gidy,0,0,HOMoments_currentx);
//			Ey = fields.getE(gidx,gidy,0,1);
//			Ez = fields.getE(gidx,gidy,0,2);

//			Ex = fields.intrpET<0,FieldData_deriv_f>(0.0,0,0,gidx,0,0);

//			printf("GPU Field[%i %i] = %e, %e, %e, %e, %e, %e\n",gidx,gidy,Bx,By,Bz,Ex,Ey,Ez);

			printf("GPU Moments[%i %i] = %e\n",gidx,gidy,Jx);

			gidx += blockDim.x*gridDim.x;
		}

		gidy += blockDim.y*gridDim.y;
	}

}

void FieldCheckGPU(FieldDataGPU fields)
{
	dim3 blockSize(1,1,1);
	dim3 gridSize(1,1,1);

	CUDA_SAFE_KERNEL((FieldCheck_g<<<1,1>>>(fields)))

}

void MomentCheckGPU(HOMomentsGPU moments)
{
	dim3 blockSize(1,1,1);
	dim3 gridSize(1,1,1);

	CUDA_SAFE_KERNEL((MomentCheck_g<<<1,1>>>(moments)))

}

void ParticleListGPUSorted::SetupBlockingInfo(int nptcls_check)
{
//	int blockSize = 256;
//	int gridSize = (nclusters+blockSize-1)/blockSize;
//
//	// First we need to null out all of the blocking information
//	// set first and last id to 0, and nptcls_cluster = 0 for
//	// each cluster
//	CUDA_SAFE_KERNEL((InitClusters<<<blockSize,gridSize>>>
//			(clusters,nclusters)));
//
//
//	// Now we need to find the first and last particles in each cluster
//	gridSize = nclusters;
//	CUDA_SAFE_KERNEL((find_cluster_boundaries<<<blockSize,gridSize>>>
//			(*this,clusters,nptcls_check)));

}


void ParticleListGPUSorted::ReorderData(int* particleIDs,int nptcls_left_old)
{
	if(sizeof(realkind) == 8)
	{
		// use 64 bit swap

		for(int i=0;i<8;i++)
		{
			realkind* idata = *get_float(i);
			realkind* odata = buffer64;

			ReOrderData64_GPU(idata,odata,particleIDs,nptcls_left_old);

			*get_float(i) = odata;
			buffer64 = idata;
		}

	}
	else
	{
		for(int i=0;i<8;i++)
		{
			realkind* idata = *get_float(i);
			realkind* odata = (realkind*)buffer32;

			ReOrderData32_GPU(idata,odata,particleIDs,nptcls_left_old);

			*get_float(i) = odata;
			buffer32 = (int*)idata;
		}
	}

	// Do the ints
	for(int i=0;i<4;i++)
	{
		int* idata = *get_int(i);
		int* odata = buffer32;

		ReOrderData32_GPU(idata,odata,particleIDs,nptcls_left_old);

		*get_int(i) = odata;
		buffer32 = idata;
	}

	int* idata = nsubcycles_current;
	int* odata = buffer32;

	ReOrderData32_GPU(idata,odata,particleIDs,nptcls_left_old);

	nsubcycles_current = odata;
	buffer32 = idata;

	idata = particleIDs_original;
	odata = buffer32;

	ReOrderData32_GPU(idata,odata,particleIDs,nptcls_left_old);

	particleIDs_original = odata;
	buffer32 = idata;
}



template<int nSpatial,int nVel,bool iEM>
void ParticleListGPUSorted::SortedGPUPushH(PlasmaData* 			pdata,
					FieldDataGPU* 			fields,
					HOMomentsGPU* 				moments)
{

//
//	int blockSize = 512;
//	int gridSize = 96;
//	int tally_size = pdata->ClusterStorDim*pdata->ClusterStorDim;
//
//	FieldDataGPU fields_t = *fields;
//
//	// Make sure everything is in constant memory
//	CUDA_SAFE_CALL(cudaMemcpyAsync(fields_d_cp,fields,
//			sizeof(FieldDataGPU),cudaMemcpyHostToDevice));
//	CUDA_SAFE_CALL(cudaMemcpyAsync(moments_d_cp,moments,
//			sizeof(HOMoments),cudaMemcpyHostToDevice));
//	CUDA_SAFE_CALL(cudaMemcpyAsync(this_cp,this,
//			sizeof(ParticleListGPU),cudaMemcpyHostToDevice));
//
//	CUDA_SAFE_CALL(cudaDeviceSynchronize());
//
////	FieldCheck(fields_t);
//
//
//	int nptcls_left = nptcls;
//	// Particles should initially be sorted due to copy from particles_old
//
//	// Set dt_finished and nsubcycles_current = 0
//	CUDA_SAFE_CALL(cudaMemset(dt_finished,0,nptcls*sizeof(realkind)));
//	CUDA_SAFE_CALL(cudaMemset(nsubcycles_current,0,nptcls*sizeof(int)));
//
//	// Basically while any particles have subcycles left we will keep pushing -> sorting
//	while(nptcls_left > 0)
//	{
//		printf("Nptcls_left in GPU push = %i\n",nptcls_left);
//		ifinished = (int*)buffer32;
//		CUDA_SAFE_CALL(cudaMemset(ifinished,0,nptcls_left*sizeof(int)));
//		// Set up the blocking information
//		printf("Setting up blocking info\n");
//		SetupBlockingInfo(nptcls_left);
//		// Push the particles
////		void GPUBlockPush(PlasmaData* 				pdata,
////							FieldDataGPU* 			fields,
////							HOMoments* 				moments,
////							ParticleListGPU	particles,
////							int* num_subcycles,
////							int tally_size)
//		int pushBlockSize = 512;
//		int pushGridSize = pdata->TBpCluster*nclusters;
//		int shared_mem_alloc = tally_size*sizeof(realkind)+sizeof(CurrentTallyGPU)+sizeof(PExitCheckGPU)+sizeof(ClusterInfo);
//		printf("Pushing particles on gpu\n");
//		CUDA_SAFE_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferShared));
//		CUDA_SAFE_KERNEL((GPUBlockPush<nSpatial,nVel,iEM><<<pushBlockSize,pushGridSize>>>(
//				pdata_d,fields_t,moments_d_cp,
//				*this,nsubcycles_thread,tally_size)));
//
//		// Get partition indices, get the number of particles left
//		int nptcls_left_old = nptcls_left;
//		nptcls_left = GenPartitionIDs(ifinished,particleIDs,nptcls_left);
//
//		// if no particles are left -> break
//		if(nptcls_left <= 0)
//			break;
//
//		// Write the cluster indices
//		CUDA_SAFE_KERNEL((write_cluster_ids<nSpatial,nVel><<<blockSize,gridSize>>>(
//				pdata_d,*this,ifinished,nptcls_left_old)));
//
//
//		// Partition the cluster indices
//		ReOrderData32_GPU(cluster_id,buffer32,particleIDs,nptcls_left_old);
//		int* temp = cluster_id;
//		cluster_id = buffer32;
//		buffer32 = temp;
//
//		// Sort the partition indices using cluster indices as keys (nptcls_left)
//		SortByKey(particleIDs,cluster_id,nptcls_left);
//
//		// Reorder the remaining particle data (nptcls_left_old)
//		ReorderData(particleIDs,nptcls_left_old);
//
//
//		// repeat
//
//	}
//
//	// Sort the particles for charge and s2 tallies
//	// Write the cluster indices
//	CUDA_SAFE_KERNEL((write_cluster_ids<nSpatial,nVel><<<blockSize,gridSize>>>(
//			pdata_d,*this,ifinished,nptcls)));
//
//
//	SortByKey(particleIDs,cluster_id,nptcls);
//
//	ReorderData(particleIDs,nptcls);
//
//	// Setup blocking info
//	SetupBlockingInfo(nptcls);
//	// Do a charge and s2 tally
//
//



//	//printf("GPU Pushing Particles\n");
//	CUDA_SAFE_KERNEL((SimpleGPUPush<nSpatial,nVel,iEM><<<blocksize,gridsize>>>
//			(particles->pdata_d,*fields,*moments,*particles,num_subcycles)));
	//printf("GPU Finished Pushing Particles\n");
}


//long long int ParticleListGPUSorted::push_interface2(PlasmaData* pdata,
//		FieldDataGPU* fields,
//		HOMoments* moments)
//{
//
//
//	int num_threads = gridsize*blocksize;
//
//
//
//	// Template Selection
//	switch(pdata->ndimensions)
//	{
//	case 1:
//		switch(pdata->nVelocity)
//		{
//		case 1:
//			if(!pdata->iEM)
//				SortedGPUPushH<1,1,0>(pdata_d,fields,moments);
//			else
//				SortedGPUPushH<1,1,1>(pdata,fields,moments);
//			break;
//		case 2:
//			if(pdata->iEM == 0)
//				SortedGPUPushH<1,2,0>(pdata,fields,moments);
//			else
//				SortedGPUPushH<1,2,1>(pdata,fields,moments);
//			break;
//		case 3:
//			if(pdata->iEM == 0)
//				SortedGPUPushH<1,3,0>(pdata,fields,moments);
//			else
//				SortedGPUPushH<1,3,1>(pdata,fields,moments);
//			break;
//		default:
//			break;
//		}
//
//		break;
//	case 2:
//		switch(pdata->nVelocity)
//		{
//		case 2:
//			if(pdata->iEM == 0)
//				SortedGPUPushH<2,2,0>(pdata,fields,moments);
//			else
//				SortedGPUPushH<2,2,1>(pdata,fields,moments);
//			break;
//		case 3:
//			if(pdata->iEM == 0)
//				SortedGPUPushH<2,3,0>(pdata,fields,moments);
//			else
//				SortedGPUPushH<2,3,1>(pdata,fields,moments);
//			break;
//		default:
//			break;
//		}
//
//		break;
//	case 3:
//		switch(pdata->nVelocity)
//		{
//		case 3:
//			if(pdata->iEM == 0)
//				SortedGPUPushH<3,3,0>(pdata,fields,moments);
//			else
//				SortedGPUPushH<3,3,1>(pdata,fields,moments);
//			break;
//		default:
//			break;
//		}
//		break;
//
//	default:
//		break;
//
//	}
//
//
////	thrust::device_ptr<int> nsubcycles_t(nsubcycles);
//
//
//	// Reduce and return the total number of particle-subcycles taken
//	long long int result = 0;// = thrust::reduce(nsubcycles_t,nsubcycles_t+num_threads);
//
//	int* nsubcycles_temp = (int*)malloc(gridsize*blocksize*sizeof(int));
//	CUDA_SAFE_CALL(cudaMemcpy(nsubcycles_temp,nsubcycles_thread,gridsize*blocksize*sizeof(int),cudaMemcpyDeviceToHost))
//
//	for(int i=0;i<gridsize*blocksize;i++)
//		result += nsubcycles_temp[i];
//
//	free(nsubcycles_temp);
//
//	return result;
//
//
//
//
//
//}








