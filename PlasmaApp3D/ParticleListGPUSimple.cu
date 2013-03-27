#include "ParticleListGPUSimple.cuh"
#include "ParticleObjNT.h"
#include "CurrentTally.h"
#include "ChargeTally.h"
#include "StressTally.h"
#include "CurrentTallyGPU.cuh"
#include "CurrentTallyCPU.h"
#include "HOMoments.h"
#include "PlasmaData.h"
#include "ProblemInitializer.h"
#include "FieldDataGPU.cuh"
#include "math.h"
#include "omp.h"
#include <thrust/count.h>

extern __constant__ PlasmaData pdata_c;


//__host__ __device__
//ParticleListGPUSimple::ParticleListGPUSimple()
//{
//
//}
//__host__ __device__
//ParticleListGPUSimple::~ParticleListGPUSimple()
//{
//	/*
//	// Free realkind arrays
//	for(int i=0;i<ParticleList_nfloats;i++)
//	{
//		free(*get_float(i));
//	}
//
//	// Allocate int arrays
//	for(int i=0;i<ParticleList_nints;i++)
//	{
//		free(*get_int(i));
//	}
//
//	// allocate short ints for cluster id's
//	free(cluster_id);
//	*/
//}

void ParticleListGPUSimple::copy_from(const ParticleList* list_in)
{
	ispecies = list_in -> ispecies;

	enum cudaMemcpyKind kind;

	if(list_in->device_type != device_type)
	{
		if(list_in->device_type != 1)
			kind = cudaMemcpyHostToDevice;
		else
			kind = cudaMemcpyDeviceToHost;
	}
	else
	{
		if(device_type == 1)
			kind = cudaMemcpyDeviceToDevice;
		else
			kind = cudaMemcpyHostToHost;
	}
	// Free realkind arrays
	for(int i=0;i<ParticleList_nfloats;i++)
	{
		CUDA_SAFE_CALL(cudaMemcpyAsync(*get_float(i),*(list_in->get_float(i)),nptcls*sizeof(realkind),kind));
	}

	// Allocate int arrays
	for(int i=0;i<ParticleList_nints;i++)
	{
		CUDA_SAFE_CALL(cudaMemcpyAsync(*get_int(i),*(list_in->get_int(i)),nptcls*sizeof(int),kind));
	}

	// allocate short ints for cluster id's
	CUDA_SAFE_CALL(cudaMemcpyAsync(cluster_id,(list_in->cluster_id),nptcls*sizeof(short int),kind));
	CUDA_SAFE_CALL(cudaMemcpyAsync(num_subcycles,(list_in->num_subcycles),nptcls*sizeof(int),kind));


	CUDA_SAFE_CALL(cudaDeviceSynchronize());

}


void ParticleListGPUSimple::allocate(PlasmaData* pdata,int nptcls_in)
{
	device_type = 1;
	// Allocate memory for particles
	nptcls_allocated = nptcls_in;

	nptcls = nptcls_in;

//	plot = gnuplot_init();

	// Allocate realkind arrays
	for(int i=0;i<ParticleList_nfloats;i++)
	{
		CUDA_SAFE_CALL(cudaMalloc((void**)get_float(i),nptcls_allocated*sizeof(realkind)));
	}

	// Allocate int arrays
	for(int i=0;i<ParticleList_nints;i++)
	{
		CUDA_SAFE_CALL(cudaMalloc((void**)get_int(i),nptcls_allocated*sizeof(int)));
	}

	// allocate short ints for cluster id's
	CUDA_SAFE_CALL(cudaMalloc((void**)&cluster_id,nptcls_allocated*sizeof(short int)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&num_subcycles,nptcls_allocated*sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&num_piccard,nptcls_allocated*sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&num_piccard2,nptcls_allocated*sizeof(double)));

	CUDA_SAFE_CALL(cudaMemset(num_piccard,0,nptcls*sizeof(double)));
	CUDA_SAFE_CALL(cudaMemset(num_piccard2,0,nptcls*sizeof(double)));

	gridsize = 6*16;
	blocksize = 256;

	gridsize = min(gridsize,(nptcls+blocksize-1)/blocksize);

	CUDA_SAFE_CALL(cudaMalloc((void**)&nsubcycles_thread,gridsize*blocksize*sizeof(int)));

	CUDA_SAFE_CALL(cudaGetSymbolAddress((void**)&pdata_d,pdata_c));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(pdata_c,pdata,sizeof(PlasmaData)));

	moments_d = new HOMoments(pdata,1);
	fields_d = new FieldDataGPU();
	fields_d -> allocate(pdata);

}

void ParticleListGPUSimple::init(ProblemInitializer* initializer, HOMoments* moments)
{

	CurrentTallyCPU currents(&moments->get_val(0,0,0,ispecies,HOMoments_currentx),
						  &moments->get_val(0,0,0,ispecies,HOMoments_currenty),
						  &moments->get_val(0,0,0,ispecies,HOMoments_currentz),
						  make_int3(moments->pdata->nx,moments->pdata->ny,moments->pdata->nz),
						  moments->pdata->dxdi,moments->pdata->dydi,moments->pdata->dzdi,
						  moments->pdata->ndimensions);

	ChargeTally charge(&moments->get_val(0,0,0,ispecies,HOMoments_charge),
						  make_int3(moments->pdata->nx,moments->pdata->ny,moments->pdata->nz),
						  moments->pdata->dxdi,moments->pdata->dydi,moments->pdata->dzdi,
						  moments->pdata->ndimensions);

	StressTally stress(&moments->get_val(0,0,0,ispecies,HOMoments_S2xx),
						  make_int3(moments->pdata->nx,moments->pdata->ny,moments->pdata->nz),
						  moments->pdata->dxdi,moments->pdata->dydi,moments->pdata->dzdi,
						  moments->pdata->ndimensions,moments->pdata->nVelocity);


	moments -> set_vals(0);


	// do this on the CPU

	// allocate temporary space on the cpu
	realkind* pxt = (realkind*)malloc(nptcls*sizeof(realkind));
	realkind* pyt = (realkind*)malloc(nptcls*sizeof(realkind));
	realkind* pzt = (realkind*)malloc(nptcls*sizeof(realkind));
	realkind* vxt = (realkind*)malloc(nptcls*sizeof(realkind));
	realkind* vyt = (realkind*)malloc(nptcls*sizeof(realkind));
	realkind* vzt = (realkind*)malloc(nptcls*sizeof(realkind));

	int* ixt = (int*)malloc(nptcls*sizeof(int));
	int* iyt = (int*)malloc(nptcls*sizeof(int));
	int* izt = (int*)malloc(nptcls*sizeof(int));

	printf("Populating Particle data\n");

	for(int i=0;i<nptcls;i++)
	{
		realkind pxtt,pytt,pztt,vxtt,vytt,vztt;
		int ixtt,iytt,iztt;


		initializer->init_particle(pxtt,pytt,pztt,ixtt,iytt,iztt,vxtt,vytt,vztt,ispecies,i);




		currents.tally(pxtt,pytt,pztt,vxtt,vytt,vztt,ixtt,iytt,iztt,1.0f);

		charge.tally(pxtt,pytt,pztt,
				ixtt,iytt,iztt,
				1.0);

		stress.tally(pxtt,pytt,pztt,
				vxtt,vytt,vztt,
				ixtt,iytt,iztt,
				1.0);



		// Set Position Values, ifloat = 0-2
		pxt[i] = pxtt;
		pyt[i] = pytt;
		pzt[i] = pztt;

		// Set Position Index Values, iint = 0-2
		ixt[i] = ixtt;
		iyt[i] = iytt;
		izt[i] = iztt;

		// Set Velocity Values, ifloat = 3-5
		vxt[i] = vxtt;
		vyt[i] = vytt;
		vzt[i] = vztt;
	}

	printf("Copying particle data to the device\n");

	CUDA_SAFE_CALL(cudaMemset(dt_finished,0,nptcls*sizeof(realkind)));



	CUDA_SAFE_CALL(cudaMemcpyAsync(px,pxt,nptcls*sizeof(realkind),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyAsync(py,pyt,nptcls*sizeof(realkind),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyAsync(pz,pzt,nptcls*sizeof(realkind),cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMemcpyAsync(vx,vxt,nptcls*sizeof(realkind),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyAsync(vy,vyt,nptcls*sizeof(realkind),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyAsync(vz,vzt,nptcls*sizeof(realkind),cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMemcpyAsync(ix,ixt,nptcls*sizeof(int),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyAsync(iy,iyt,nptcls*sizeof(int),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyAsync(iz,izt,nptcls*sizeof(int),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemset(num_subcycles,0,nptcls*sizeof(int)));
	CUDA_SAFE_CALL(cudaDeviceSynchronize());


	free(pxt);
	free(pyt);
	free(pzt);
	free(vxt);
	free(vyt);
	free(vzt);
	free(ixt);
	free(iyt);
	free(izt);
}

realkind ParticleListGPUSimple::evaluate_energy(PlasmaData* pdata)
{
	double etotal = 0.0;
	for(int i=0;i<nptcls;i++)
	{
//		etotal += get_fvalue(i,3)* get_fvalue(i,3);
//		etotal += get_fvalue(i,4)* get_fvalue(i,4);
//		etotal += get_fvalue(i,5)* get_fvalue(i,5);
	}

	etotal = etotal * pdata->mspecies[ispecies] * 0.5/((double)pdata->nptcls_total);

	return etotal;
}

double4 ParticleListGPUSimple::subcycle_stats(PlasmaData* pdata)
{
	int* num_subcycles_temp = (int*)malloc(nptcls*sizeof(int));

	CUDA_SAFE_CALL(cudaMemcpy(num_subcycles_temp,num_subcycles,nptcls*sizeof(int),
			cudaMemcpyDeviceToHost));

	double scale = pdata->npiccard_outer;
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

double4 ParticleListGPUSimple::piccard_stats(PlasmaData* pdata)
{
	int* num_subcycles_temp = (int*)malloc(nptcls*sizeof(int));
	double* num_piccard_temp = (double*)malloc(nptcls*sizeof(double));
	double* num_piccard2_temp = (double*)malloc(nptcls*sizeof(double));

	CUDA_SAFE_CALL(cudaMemcpy(num_subcycles_temp,num_subcycles,nptcls*sizeof(int),
			cudaMemcpyDeviceToHost));

	CUDA_SAFE_CALL(cudaMemcpy(num_piccard_temp,num_piccard,nptcls*sizeof(double),
			cudaMemcpyDeviceToHost));

	CUDA_SAFE_CALL(cudaMemcpy(num_piccard2_temp,num_piccard2,nptcls*sizeof(double),
			cudaMemcpyDeviceToHost));

	double scale = pdata->npiccard_outer;
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
	free(num_piccard2_temp);

	return make_double4(mean,std_diff,mins,maxs);
}



void ParticleListGPUSimple::plot_particles(PlasmaData* pdata)
{
	float* x_vals = (float*)malloc(nptcls*sizeof(float));
	float* y_vals = (float*)malloc(nptcls*sizeof(float));

	realkind* px_t = (realkind*)malloc(nptcls*sizeof(realkind));
	int* ix_t = (int*)malloc(nptcls*sizeof(int));
	realkind* vx_t = (realkind*)malloc(nptcls*sizeof(realkind));

	CUDA_SAFE_CALL(cudaMemcpy(px_t,px,nptcls*sizeof(realkind),cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(ix_t,ix,nptcls*sizeof(int),cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(vx_t,vx,nptcls*sizeof(realkind),cudaMemcpyDeviceToHost));

	for(int i=0;i<nptcls;i++)
	{
		x_vals[i] = (px_t[i]+ix_t[i])*pdata->dxdi + pdata->xmin;
		y_vals[i] = vx_t[i];

//		printf("particle[%i]: %f %f\n",i,x_vals[i],y_vals[i]);
	}

	gnuplot_resetplot(plot);

	gnuplot_plot_xy(plot,x_vals,y_vals,nptcls,NULL);


	free(x_vals);
	free(y_vals);
	free(px_t);
	free(vx_t);
	free(ix_t);
}

void ParticleListGPUSimple::CPUfree()
{
	// Allocate realkind arrays
	for(int i=0;i<ParticleList_nfloats;i++)
	{
		CUDA_SAFE_CALL(cudaFree(*get_float(i)));
	}

	// Allocate int arrays
	for(int i=0;i<ParticleList_nints;i++)
	{
		CUDA_SAFE_CALL(cudaFree(*get_int(i)));
	}

	// allocate short ints for cluster id's
	CUDA_SAFE_CALL(cudaFree(cluster_id));
}

long long int ParticleListGPUSimple::push(PlasmaData* pdata, FieldData* fields, HOMoments* moments)
{
	long long int result;
	// Copy field data to the device
//	printf("Copying fields to device\n");
	fields_d->copy_from(fields);

//	printf("Nulling device HOMoments\n");
	moments_d -> set_vals(0);

	// Change location of pdata in moments_d to device
	moments_d->pdata = pdata_d;
	fields_d->pdata = pdata_d;

	// Push the particles
//	printf("Pushing particles on the GPU\n");
	result = push_interface3(pdata,fields_d,moments_d);
//	printf("More Pushing particles on the GPU\n");
	// Change the location of pdata in moments_d to host
	moments_d->pdata = pdata;
	fields_d->pdata = pdata;
	// Copy HO moments to the host
//	printf("Copying device HOMoments to Host\n");
	moments->copy_from(moments_d);

	return result;
}




template<int nSpatial,int nVel,bool iEM> __global__
void SimpleGPUPush(PlasmaData* 				pdata,
					FieldDataGPU 			fields,
					HOMoments 				moments2,
					ParticleListGPUSimple	particles,
					int* num_subcycles)
{
	int tidx = threadIdx.x;
	int bidx = blockIdx.x;
	int gidx = tidx+blockDim.x*bidx;
	int stride = blockDim.x*gridDim.x;

	int pid = gidx;

	long long int num_subcycles_thread = 0;

	__shared__ float currentx[256];
	__shared__ realkind charge_s[256];
	__shared__ realkind S2xx_s[256];
	while(tidx < 256)
	{
		currentx[tidx] = 0;
		charge_s[tidx] = 0;
		S2xx_s[tidx] = 0;
		tidx += blockDim.x;
	}
	tidx = threadIdx.x;


	ParticleListGPUSimple particles2 = particles;
	FieldDataGPU 			fields2 = fields;
	FieldData*				fields3 = &fields2;
	HOMoments 				moments = moments2;
//	PlasmaData pdata = *pdata2;

	ParticleObjNT<1,nSpatial,nVel,iEM> particle(&pid);
	typevecN<int,1> iter;

	CurrentTallyGPU currents((float*)currentx, (float*)&moments.get_val(0,0,0,particles.ispecies,HOMoments_currenty),
						  (float*)&moments.get_val(0,0,0,particles.ispecies,HOMoments_currentz),
						  pdata->nx,pdata->ny,pdata->nz,
						  0,0,0,
						  pdata->ndimensions);

	ChargeTally charge(charge_s,
						  make_int3(currents.nx,currents.ny,currents.nz),
						  pdata->dxdi,pdata->dydi,pdata->dzdi,
						  1);

	StressTally stress(S2xx_s,
			  make_int3(currents.nx,currents.ny,currents.nz),
			  pdata->dxdi,pdata->dydi,pdata->dzdi,
			  1,1);

	particle.species = particles2.ispecies;


	while(pid < particles.nptcls)
	{
		iter(0) = 0;
		printf("reading in particle %i\n",pid);
		particle.copy_in_gpu(particles2,0);
		particle.dt_finished(0) = 0;
		particle.push(pdata,fields3,&currents,iter,pdata->nSubcycle_max);

//		printf("Writing Paerticles Back\n");
		particle.write_back(particles2,0);

		particles2.num_subcycles[pid] += iter(0);

		charge.tally1d(particles2.px[pid],particles2.ix[pid],1.0);
		stress.tally1d1v(particles2.px[pid],particles2.vx[pid],particles2.ix[pid],1.0);

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
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_charge),charge_s[pid]);
		atomicAddD(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xx),S2xx_s[pid]);
#else
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_currentx),currentx[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_charge),charge_s[pid]);
		atomicAdd(&moments.get_val(pid,0,0,particles.ispecies,HOMoments_S2xx),S2xx_s[pid]);
#endif
		pid += blockDim.x;
	}

}

template<int nSpatial,int nVel,bool iEM>
void SimpleGPUPushH(PlasmaData* 			pdata,
					FieldDataGPU* 			fields,
					HOMoments* 				moments,
					ParticleListGPUSimple*	particles,
					int* num_subcycles)
{
	int blocksize = particles->blocksize;
	int gridsize = particles->gridsize;
	//printf("GPU Pushing Particles\n");
	CUDA_SAFE_KERNEL((SimpleGPUPush<nSpatial,nVel,iEM><<<blocksize,gridsize>>>
			(particles->pdata_d,*fields,*moments,*particles,num_subcycles)));
	//printf("GPU Finished Pushing Particles\n");
}


long long int ParticleListGPUSimple::push_interface3(PlasmaData* pdata,
		FieldDataGPU* fields,
		HOMoments* moments)
{


	int num_threads = gridsize*blocksize;



	// Template Selection
	switch(pdata->ndimensions)
	{
	case 1:
		switch(pdata->nVelocity)
		{
		case 1:
			if(!pdata->iEM)
				SimpleGPUPushH<1,1,0>(pdata_d,fields,moments,this,nsubcycles_thread);
			else
				SimpleGPUPushH<1,1,1>(pdata,fields,moments,this,nsubcycles_thread);
			break;
		case 2:
			if(pdata->iEM == 0)
				SimpleGPUPushH<1,2,0>(pdata,fields,moments,this,nsubcycles_thread);
			else
				SimpleGPUPushH<1,2,1>(pdata,fields,moments,this,nsubcycles_thread);
			break;
		case 3:
			if(pdata->iEM == 0)
				SimpleGPUPushH<1,3,0>(pdata,fields,moments,this,nsubcycles_thread);
			else
				SimpleGPUPushH<1,3,1>(pdata,fields,moments,this,nsubcycles_thread);
			break;
		default:
			break;
		}

		break;
	case 2:
		switch(pdata->nVelocity)
		{
		case 2:
			if(pdata->iEM == 0)
				SimpleGPUPushH<2,2,0>(pdata,fields,moments,this,nsubcycles_thread);
			else
				SimpleGPUPushH<2,2,1>(pdata,fields,moments,this,nsubcycles_thread);
			break;
		case 3:
			if(pdata->iEM == 0)
				SimpleGPUPushH<2,3,0>(pdata,fields,moments,this,nsubcycles_thread);
			else
				SimpleGPUPushH<2,3,1>(pdata,fields,moments,this,nsubcycles_thread);
			break;
		default:
			break;
		}

		break;
	case 3:
		switch(pdata->nVelocity)
		{
		case 3:
			if(pdata->iEM == 0)
				SimpleGPUPushH<3,3,0>(pdata,fields,moments,this,nsubcycles_thread);
			else
				SimpleGPUPushH<3,3,1>(pdata,fields,moments,this,nsubcycles_thread);
			break;
		default:
			break;
		}
		break;

	default:
		break;

	}


//	thrust::device_ptr<int> nsubcycles_t(nsubcycles);


	// Reduce and return the total number of particle-subcycles taken
	long long int result = 0;// = thrust::reduce(nsubcycles_t,nsubcycles_t+num_threads);

	int* nsubcycles_temp = (int*)malloc(gridsize*blocksize*sizeof(int));
	CUDA_SAFE_CALL(cudaMemcpy(nsubcycles_temp,nsubcycles_thread,gridsize*blocksize*sizeof(int),cudaMemcpyDeviceToHost))

	for(int i=0;i<gridsize*blocksize;i++)
		result += nsubcycles_temp[i];

	free(nsubcycles_temp);

	return result;





}















