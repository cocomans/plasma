/*-------------------------------------------------------------------------*/
/**
	@file		ParticleListCPU.cpp
	@author	J. Payne
	@date		1/09/2012
	@brief	Defines members of the ParticleListCPU class


*/
/*--------------------------------------------------------------------------*/
#include "ParticleListCPUSorted.h"
#include "ParticleObjNT.h"
#include "ParticleObjN.h"
#include "CurrentTallyCPU.h"
#include "ChargeTally.h"
#include "StressTally.h"
#include "HOMoments.h"
#include "PlasmaData.h"
#include "ProblemInitializer.h"
#include "math.h"
#include <sstream>
#include <omp.h>
#include <utmpx.h>
#include "CPUTimer.h"
#include <vector>



ParticleListCPUSorted::ParticleListCPUSorted()
{
	device_type = 0;
}

ParticleListCPUSorted::~ParticleListCPUSorted()
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

void ParticleListCPUSorted::copy_stats_from(const ParticleList* list_in)
{
	memcpy(num_subcycles,list_in->num_subcycles,nptcls*sizeof(int));

	memcpy(num_piccard,list_in->num_piccard,nptcls*sizeof(double));
}

void ParticleListCPUSorted::copy_from(const ParticleList* list_in)
{
	ispecies = list_in -> ispecies;
	// Free realkind arrays

	if(list_in -> device_type == 0){
		for(int i=0;i<ParticleList_nfloats;i++)
		{
			memcpy(*get_float(i),*(list_in->get_float(i)),nptcls*sizeof(realkind));
		}

		// Allocate int arrays
		for(int i=0;i<ParticleList_nints;i++)
		{
			memcpy(*get_int(i),*(list_in->get_int(i)),nptcls*sizeof(int));
		}

		// allocate short ints for cluster id's
		memcpy(cluster_id,list_in->cluster_id,nptcls*sizeof(int));

//		memcpy(num_subcycles,list_in->num_subcycles,nptcls*sizeof(int));
//
//		memcpy(num_piccard,list_in->num_piccard,nptcls*sizeof(double));
//		memcpy(num_piccard2,list_in->num_piccard2,nptcls*sizeof(double));
	}
	else if(list_in->device_type == 1)
	{

#ifndef NO_CUDA
		enum cudaMemcpyKind kind = cudaMemcpyDeviceToHost;
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
		CUDA_SAFE_CALL(cudaMemcpyAsync(cluster_id,(list_in->cluster_id),nptcls*sizeof(int),kind));

		CUDA_SAFE_CALL(cudaDeviceSynchronize());
#endif
	}
}

void ParticleListCPUSorted::allocate(PlasmaData* pdata,int nptcls_in)
{
	//printf("Allocating Particle List on the CPU\n");

//	if(pdata->plot_flag)
//	plot = gnuplot_init();

	//gnuplot_cmd(plot,"set pointsize 0.1");
	// Allocate memory for particles
	nptcls_allocated = nptcls_in;

	nptcls = nptcls_in;

	num_cores = pdata->num_cores;

	// Allocate realkind arrays
	for(int i=0;i<ParticleList_nfloats;i++)
	{
		*get_float(i) = (realkind*)malloc(nptcls_allocated*sizeof(realkind));
	}

	// Allocate int arrays
	for(int i=0;i<ParticleList_nints;i++)
	{
		*get_int(i) = (int*)malloc(nptcls_allocated*sizeof(realkind));
	}

	buffer = (realkind*)malloc(nptcls_allocated*sizeof(realkind));

	num_subcycles = (int*)malloc(nptcls_allocated*sizeof(realkind));
	memset(num_subcycles,0,nptcls_allocated*sizeof(int));

	num_piccard = (realkind*)malloc(nptcls_allocated*sizeof(double));
	memset(num_piccard,0,nptcls_allocated*sizeof(double));

	num_piccard2 = (realkind*)malloc(nptcls_allocated*sizeof(double));
	memset(num_piccard2,0,nptcls_allocated*sizeof(double));


	// allocate short ints for cluster id's
	cluster_id = (int*)malloc(nptcls_allocated*sizeof(int));
	ptcl_index = (int*)malloc(nptcls_allocated*sizeof(int));

	piccard_timer = (CPUTimer*)malloc(pdata->num_cores*sizeof(CPUTimer));
	accel_timer = (CPUTimer*)malloc(pdata->num_cores*sizeof(CPUTimer));
	tally_timer = (CPUTimer*)malloc(pdata->num_cores*sizeof(CPUTimer));
	crossing_timer = (CPUTimer*)malloc(pdata->num_cores*sizeof(CPUTimer));
	dtau_est_timer = (CPUTimer*)malloc(pdata->num_cores*sizeof(CPUTimer));
	tally_timer2 = (CPUTimer*)malloc(pdata->num_cores*sizeof(CPUTimer));
	load_store_timer = (CPUTimer*)malloc(pdata->num_cores*sizeof(CPUTimer));

	int tid;
	omp_set_num_threads(pdata->num_cores);
#pragma omp parallel private(tid) default(shared) num_threads(pdata->num_cores)
	{
		tid = omp_get_thread_num();
		piccard_timer[tid] = *(new CPUTimer());
		accel_timer[tid] = *(new CPUTimer());
		tally_timer[tid] = *(new CPUTimer());
		crossing_timer[tid] = *(new CPUTimer());
		dtau_est_timer[tid] = *(new CPUTimer());

		tally_timer2[tid] = *(new CPUTimer());
		load_store_timer[tid] = *(new CPUTimer());
	}

	push_timer = new CPUTimer();


}

void ParticleListCPUSorted::init(ProblemInitializer* initializer, HOMoments* moments)
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
			&moments->get_val(0,0,0,ispecies,HOMoments_S2xy),
			&moments->get_val(0,0,0,ispecies,HOMoments_S2xz),
			&moments->get_val(0,0,0,ispecies,HOMoments_S2yy),
			&moments->get_val(0,0,0,ispecies,HOMoments_S2yz),
			&moments->get_val(0,0,0,ispecies,HOMoments_S2zz),
						  moments->pdata->nx,moments->pdata->ny,moments->pdata->nz,
						  moments->pdata->ndimensions,moments->pdata->nVelocity);

	moments -> set_vals(0);
#pragma omp for
	for(int i=0;i<nptcls;i++)
	{
		realkind px,py,pz,vx,vy,vz;
		int ix,iy,iz;

		initializer->init_particle(px,py,pz,ix,iy,iz,vx,vy,vz,ispecies,i);





		dt_finished[i] = 0;


		// Set Position Values, ifloat = 0-2gmake
		this->get_fvalue(i,0) = px;
		this->get_fvalue(i,1) = py;
		this->get_fvalue(i,2) = pz;

		// Set Position Index Values, iint = 0-2
		this->get_ivalue(i,0) = ix;
		this->get_ivalue(i,1) = iy;
		this->get_ivalue(i,2) = iz;

		// Set Velocity Values, ifloat = 3-5
		this->get_fvalue(i,3) = vx;
		this->get_fvalue(i,4) = vy;
		this->get_fvalue(i,5) = vz;
	}

	for(int i=0;i<nptcls;i++)
	{
	currents.tally(px[i],py[i],pz[i],vx[i],vy[i],vz[i],ix[i],iy[i],iz[i],1.0);

	charge.tally(px[i],py[i],pz[i],
			ix[i],iy[i],iz[i],
			1.0);

	stress.tally(px[i],py[i],pz[i],
			vx[i],vy[i],vz[i],
			ix[i],iy[i],iz[i],
			1.0);
	}

	memset(num_subcycles,0,nptcls*sizeof(int));

}

realkind ParticleListCPUSorted::evaluate_energy(PlasmaData* pdata)
{
	double etotal = 0.0;
	for(int i=0;i<nptcls;i++)
	{
		etotal += get_fvalue(i,3)* get_fvalue(i,3);
	//	etotal += get_fvalue(i,4)* get_fvalue(i,4);
	//	etotal += get_fvalue(i,5)* get_fvalue(i,5);
	}

	etotal = etotal * pdata->mspecies[ispecies] * 0.5 * pdata->wspecies[ispecies];

	return etotal;
}

void ParticleListCPUSorted::plot_particles(PlasmaData* pdata)
{
	float* x_vals = (float*)malloc(nptcls*sizeof(realkind));
	float* y_vals = (float*)malloc(nptcls*sizeof(realkind));


	for(int i=0;i<nptcls;i++)
	{
		x_vals[i] = (px[i]+ix[i])*pdata->dxdi + pdata->xmin;
		//x_vals[i] = (py[i]+iy[i])*pdata->dydi + pdata->ymin;
		y_vals[i] = vx[i];
//		printf("particle[%i]: %f %f\n",i,x_vals[i],y_vals[i]);

	}

	gnuplot_resetplot(plot);

	gnuplot_plot_xy(plot,x_vals,y_vals,nptcls,NULL);


	free(x_vals);
	free(y_vals);
}

double4 ParticleListCPUSorted::subcycle_stats(PlasmaData* pdata)
{
	double scale = pdata->npiccard_outer;
	double mean = 0;
	double mean2 = 0;
	double mins = num_subcycles[0]/scale;
	double maxs = num_subcycles[0]/scale;
	int imax = 0;
	int imin = 0;
	double nsubcycles_total = 0;
	for(int i=0;i<nptcls;i++)
	{
		if(mins > num_subcycles[i]/scale)
		{
			mins = num_subcycles[i]/scale;
			imin = i;
		}
		if(maxs < num_subcycles[i]/scale)
		{
			maxs = num_subcycles[i]/scale;
			imax = i;
		}

		nsubcycles_total += num_subcycles[i];
		mean += num_subcycles[i]/((double)nptcls*scale);
		//mean2 += num_subcycles[i]*num_subcycles[i]/((double)nptcls*scale*scale*nptcls);
		mean2 = mean*mean;
	}

	double std_diff;

	std_diff = sqrt(fabs(mean*mean - mean2)/((double)nptcls*scale));
	printf("Particle Subcycle Stats:\n");
	printf("Avg Subcycles: %f +/- %f\n",mean,std_diff);
	printf("Min / Max: %f[%i] / %f[%i]\n",mins,imin,maxs,imax);
	printf("Total number of subcycles was %e\n",nsubcycles_total);

	return make_double4(mean,std_diff,mins,maxs);
}

double4 ParticleListCPUSorted::piccard_stats(PlasmaData* pdata)
{
	double scale = pdata->npiccard_outer;
	double mean = 0;
	double mean2 = 0;
	double mins = num_piccard[0]/num_subcycles[0];
	double maxs = num_piccard[0]/num_subcycles[0];
	int imax = 0;
	int imin = 0;

	for(int i=0;i<nptcls;i++)
	{
		if(mins > num_piccard[i]/num_subcycles[i])
		{
			mins = num_piccard[i]/num_subcycles[i];
			imin = i;
		}
		if(maxs < num_piccard[i]/num_subcycles[i])
		{
			maxs = num_piccard[i]/num_subcycles[i];
			imax = i;
		}

		mean += num_piccard[i]/((double)nptcls*num_subcycles[i]);
		mean2 += num_piccard[i]*num_piccard[i]/((double)nptcls*num_subcycles[i]*num_subcycles[i]);
	}

	double std_diff;

	std_diff = sqrt(fabs(mean*mean - mean2));
	printf("Particle Piccard Stats(CPU):\n");
	printf("Avg Piccard: %f +/- %f\n",mean,std_diff);
	printf("Min / Max: %f[%i] / %f[%i]\n",mins,imin,maxs,imax);

	return make_double4(mean,std_diff,mins,maxs);
}

/*-------------------------------------------------------------------------*/
/**
	@brief Return the cummulative time for the requested timer
*/
/*--------------------------------------------------------------------------*/
double ParticleListCPUSorted::get_cummulative_time(int itimer)
{
	CPUTimer* temp_timer;
	switch(itimer)
	{
	case 0:
		temp_timer = piccard_timer;
		break;
	case 1:
		temp_timer = accel_timer;
		break;
	case 2:
		temp_timer = tally_timer;
		break;
	case 3:
		temp_timer = crossing_timer;
		break;
	case 4:
		temp_timer = dtau_est_timer;
		break;
	case 5:
		temp_timer = tally_timer2;
		break;
	case 6:
		temp_timer = load_store_timer;
		break;
	case 7:
		// Case 7 = piccard other
		break;
	case 8:
		temp_timer = push_timer;
		break;
	case 9:
		// case 9 = push other
		break;
	default:
		temp_timer = piccard_timer;
		break;
	}



	double cummulative = 0;

	if(itimer == 8)
		cummulative = push_timer->get_cummulative() * num_cores;
	else if(itimer == 9)
		cummulative = (push_timer->get_cummulative()
		- (get_cummulative_time(6)
			+ get_cummulative_time(5)
			+ get_cummulative_time(0)))*num_cores;
	else
	for(int i=0;i<num_cores;i++)
	{
		if(itimer == 6)
			cummulative += load_store_timer[i].get_cummulative() - piccard_timer[i].get_cummulative();
		else if(itimer == 7)
			cummulative += piccard_timer[i].get_cummulative()
			- (tally_timer[i].get_cummulative() + crossing_timer[i].get_cummulative() + accel_timer[i].get_cummulative());
		else
			cummulative += temp_timer[i].get_cummulative();
	}



	printf("Cummulative time = %e\n",cummulative);

	return cummulative/((double)num_cores);
}

void ParticleListCPUSorted::CPUfree()
{
	// Allocate realkind arrays
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
}


template<typename T2>
void swap(int *a, T2 *vala, int *b, T2 *valb)
{
    int t;
    T2 tval;
    tval = *vala;
    t = *a;
    *a = *b;
    *vala = *valb;
    *b = t;
    *valb = tval;
}



template<typename T2>
void bitonic_sort_par(int start, int length, int *seq, T2 *vals, int m,int flag)
{
    int i;
    int split_length;

    if (length == 1)
        return;

    if (length % 2 !=0 )
    {
        printf("The length of a (sub)sequence is not divided by 2.\n");
        exit(0);
    }

    split_length = length / 2;

    // bitonic split
#pragma omp parallel for shared(seq,vals, flag, start, split_length) private(i)
    for (i = start; i < start + split_length; i++)
    {
        if (flag == 0)
        {
            if (seq[i] > seq[i + split_length])
                swap(&seq[i], &seq[i + split_length],vals+i,vals+i+split_length);
        }
        else
        {
            if (seq[i] < seq[i + split_length])
                swap(&seq[i], &seq[i + split_length],vals+i,vals+i+split_length);
        }
    }

    if (split_length > m)
    {
        // m is the size of sub part-> n/numThreads
        bitonic_sort_par(start, split_length, seq, vals,m, flag);
        bitonic_sort_par(start + split_length, split_length, seq, vals,m, flag);
    }

    return;
}
template<typename T2>
void bitonic_sort_seq(int start, int length, int *seq, T2 *vals, int m,int flag)
{
    int i;
    int split_length;

    if (length == 1)
        return;

    if (length % 2 !=0 )
    {
        printf("error\n");
        exit(0);
    }

    split_length = length / 2;

    // bitonic split
    for (i = start; i < start + split_length; i++)
    {
        if (flag == 0)
        {
            if (seq[i] > seq[i + split_length])
            	swap(&seq[i], &seq[i + split_length],vals+i,vals+i+split_length);
        }
        else
        {
            if (seq[i] < seq[i + split_length])
            	swap(&seq[i], &seq[i + split_length],vals+i,vals+i+split_length);
        }
    }

    bitonic_sort_seq(start, split_length, seq,vals,m, flag);
    bitonic_sort_seq(start + split_length, split_length, seq,vals,m, flag);
}

template<typename T2>
void bitonic_sort_by_key(int* keys, T2* vals,int n,int numThreads)
{

	omp_set_num_threads(numThreads);

	  int m,i,j,flag;
	    // making sure input is okay
	    if ( n < numThreads * 2 )
	    {
	        printf("The size of the sequence is less than 2 * the number of processes.\n");
	        exit(0);
	    }

	    // the size of sub part
	    m = n / numThreads;

	    // make the sequence bitonic - part 1
	    for (i = 2; i <= m; i = 2 * i)
	    {
	#pragma omp parallel for shared(i, keys,vals) private(j, flag)
	        for (j = 0; j < n; j += i)
	        {
	            if ((j / i) % 2 == 0)
	                flag = 0;
	            else
	                flag = 1;
	            bitonic_sort_seq(j, i, keys,vals,m, flag);
	        }
	    }

	    // make the sequence bitonic - part 2
	    for (i = 2; i <= numThreads; i = 2 * i)
	    {
	        for (j = 0; j < numThreads; j += i)
	        {
	            if ((j / i) % 2 == 0)
	                flag = 0;
	            else
	                flag = 1;
	            bitonic_sort_par(j*m, i*m, keys,vals,m, flag);
	        }
	#pragma omp parallel for shared(j)
	        for (j = 0; j < numThreads; j++)
	        {
	            if (j < i)
	                flag = 0;
	            else
	                flag = 1;
	            bitonic_sort_seq(j*m, m, keys,vals,m, flag);
	        }
	    }
}

int compare(const void* a, const void* b)
{
	return (((int2*)a)->x < ((int2*)b)->x);
}

void stupid_sort(int* keys_h, int* values_h, int nelements)
{

//	int* keys_h = (int*)malloc(nelements*sizeof(int));
//	int* values_h = (int*)malloc(nelements*sizeof(int));

	int2* dict = (int2*)malloc(nelements*sizeof(int2));

//	CUDA_SAFE_CALL(cudaMemcpy(keys_h,keys_d,nelements*sizeof(int),cudaMemcpyDeviceToHost));
//	CUDA_SAFE_CALL(cudaMemcpy(values_h,values_d,nelements*sizeof(int),cudaMemcpyDeviceToHost));

	for(int i=0;i<nelements;i++)
	{
		dict[i].x = keys_h[i];
		dict[i].y = values_h[i];
	}

	qsort(dict,nelements,sizeof(int2),compare);

	for(int i=0;i<nelements;i++)
	{
		keys_h[i] = dict[i].x;
		values_h[i] = dict[i].y;
	}

//	CUDA_SAFE_CALL(cudaMemcpy(keys_d,keys_h,nelements*sizeof(int),cudaMemcpyHostToDevice));
//	CUDA_SAFE_CALL(cudaMemcpy(values_d,values_h,nelements*sizeof(int),cudaMemcpyHostToDevice));

//	free(keys_h);
//	free(values_h);
	free(dict);
}

template<typename T>
void reorder_by_index(T* odata, T* idata,int* index_array,int n)
{
#pragma omp parallel for
		for(int i=0;i<n;i++)
		{
			int ogidx = index_array[i];
//			printf("ogidx[%i] = %i\n",i,ogidx);
			odata[i] = idata[ogidx];
		}

}

__inline__
int part1by1(int ix)
{
	ix = (ix | (ix << 16)) & 0x030000FF;
	ix = (ix | (ix <<  8)) & 0x0300F00F;
	ix = (ix | (ix <<  4)) & 0x030C30C3;
	ix = (ix | (ix <<  2)) & 0x09249249;

    return ix;
}

__inline__
int zorder2D(int ix,int iy)
{


	return part1by1(ix) | (part1by1(iy) << 1);

}

void ParticleListCPUSorted::sort(PlasmaData* pdata)
{
	printf("Sorting Particles\n");
	int nthreads = pdata->num_cores;
	omp_set_num_threads(nthreads);

	int subcycle_res = 5;
	int nx_bin = 2;
	int nsubcycle_bins = (pdata->nSubcycle_max+subcycle_res-1)/subcycle_res;
	int nxindex = (pdata->nx+nx_bin-1)/nx_bin;

	// Evaluate and fill cluster id's
#pragma omp parallel for
	for(int i=0;i<nptcls;i++)
	{
		int xindex = ix[i]/nx_bin;
		int yindex = iy[i]/nx_bin;
		int subindex = floor((double)num_subcycles[i]/(pdata->npiccard_outer*subcycle_res));
		cluster_id[i] = (xindex+nxindex*yindex);
		ptcl_index[i] = i;
	}

	// Sort the particle indicies using the cluster_id's
	printf("Sorting Particle ID's\n");
	stupid_sort(cluster_id,ptcl_index,nptcls);

//	printf("Reording Particles floats\n");
	for(int i=0;i<7;i++)
	{
		realkind* idata = *(get_float(i));
		realkind* odata = buffer;

		printf("Reording Particles floats\n");
		reorder_by_index<realkind>(odata,idata,ptcl_index,nptcls);

		*(get_float(i)) = odata;
		buffer = idata;
	}


	printf("Reording Particles ints\n");
	for(int i=0;i<3;i++)
	{
		int* idata = *(get_int(i));
		int* odata = (int*)buffer;

		reorder_by_index<int>(odata,idata,ptcl_index,nptcls);

		*(get_int(i)) = odata;
		buffer = (realkind*)idata;
	}

	reorder_by_index<realkind>(buffer,num_piccard,ptcl_index,nptcls);
	memcpy(num_piccard,buffer,nptcls*sizeof(realkind));

	reorder_by_index<int>((int*)buffer,num_subcycles,ptcl_index,nptcls);
	memcpy(num_subcycles,buffer,nptcls*sizeof(int));



}

/*-------------------------------------------------------------------------*/
/**
	@brief Interface to the particle physics object on the CPU
	@param[in] pdata Simulation information
	@param[in] fields Field Data
	@param[in,out] momments HO moments.

	@return total number of particle subcycles executed in push() call
*/
/*--------------------------------------------------------------------------*/
long long int ParticleListCPUSorted::push(PlasmaData* pdata, FieldData* fields, HOMoments* moments)
{

	sort(pdata);
	push_timer->start();
	printf("Actually Pushing Particles\n");
	long long int result = push_interface<VEC_LENGTH_MAX>(pdata,fields,moments);
	push_timer->stop();

	return result;
}

template<int VEC_LENGTH>
long long int ParticleListCPUSorted::push_interface(PlasmaData* pdata, FieldData* fields, HOMoments* moments)
{
	long long int result = 0;
	if(pdata->cpu_vec_length == VEC_LENGTH)
	{
		switch(pdata->ndimensions)
		{
		case 1:
			switch(pdata->nVelocity)
			{
			case 1:
				if(!pdata->iEM)
					result = pushT<VEC_LENGTH,1,1,0>(pdata,fields,moments);
				else
					result = pushT<VEC_LENGTH,1,1,1>(pdata,fields,moments);
				break;
			case 2:
				if(pdata->iEM == 0)
					result = pushT<VEC_LENGTH,1,2,0>(pdata,fields,moments);
				else
					result = pushT<VEC_LENGTH,1,2,1>(pdata,fields,moments);
				break;
			case 3:
				if(pdata->iEM == 0)
					result = pushT<VEC_LENGTH,1,3,0>(pdata,fields,moments);
				else
					result = pushT<VEC_LENGTH,1,3,1>(pdata,fields,moments);
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
					result = pushT<VEC_LENGTH,2,2,0>(pdata,fields,moments);
				else
					result = pushT<VEC_LENGTH,2,2,1>(pdata,fields,moments);
				break;
			case 3:
				if(pdata->iEM == 0)
					result = pushT<VEC_LENGTH,2,3,0>(pdata,fields,moments);
				else
					result = pushT<VEC_LENGTH,2,3,1>(pdata,fields,moments);
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
					result = pushT<VEC_LENGTH,3,3,0>(pdata,fields,moments);
				else
					result = pushT<VEC_LENGTH,3,3,1>(pdata,fields,moments);
				break;
			default:
				break;
			}
			break;

		default:
			break;

		}



	}
	else
	{
		result = push_interface<VEC_LENGTH-1>(pdata,fields,moments);
	}

	return result;
}


template<>
long long int ParticleListCPUSorted::push_interface<1>(PlasmaData* pdata, FieldData* fields, HOMoments* moments)
{
	long long int result = 0;
	if(pdata->cpu_vec_length == 1)
	{
		switch(pdata->ndimensions)
		{
		case 1:
			switch(pdata->nVelocity)
			{
			case 1:
				if(!pdata->iEM)
					result = pushT<1,1,1,0>(pdata,fields,moments);
				else
					result = pushT<1,1,1,1>(pdata,fields,moments);
				break;
			case 2:
				if(pdata->iEM == 0)
					result = pushT<1,1,2,0>(pdata,fields,moments);
				else
					result = pushT<1,1,2,1>(pdata,fields,moments);
				break;
			case 3:
				if(pdata->iEM == 0)
					result = pushT<1,1,3,0>(pdata,fields,moments);
				else
					result = pushT<1,1,3,1>(pdata,fields,moments);
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
					result = pushT<1,2,2,0>(pdata,fields,moments);
				else
					result = pushT<1,2,2,1>(pdata,fields,moments);
				break;
			case 3:
				if(!pdata->iEM)
					result = pushT<1,2,3,0>(pdata,fields,moments);
				else
					result = pushT<1,2,3,1>(pdata,fields,moments);
				break;
			default:
				break;
			}

			break;
		case 3:
			switch(pdata->nVelocity)
			{
			case 3:
				if(!pdata->iEM)
					result = pushT<1,3,3,0>(pdata,fields,moments);
				else
					result = pushT<1,3,3,1>(pdata,fields,moments);
				break;
			default:
				break;
			}
			break;

		default:
			break;

		}
	}

	return result;

}

template<int VEC_LENGTH,const int nSpatial,const int nVel,const bool iEM> __attribute__ ((noinline))
long long int ParticleListCPUSorted::pushT(PlasmaData* pdata, FieldData* fields, HOMoments* moments)
{

	int tid;
	int nthreads = pdata->num_cores;
	int stride = (nptcls+nthreads-1)/nthreads;

	long long int nSubSteps_proc[nthreads];

	omp_set_num_threads(nthreads);

//	for(int i=0;i<pdata->nx;i++)
//	{
//		realkind temp;
//		temp = fields->intrpE(0.5,0,0,i,0,0,0,FieldData_deriv_f);
//		printf("fields[%i] on cpu = %f\n",i,temp);
//	}

	//printf("particles ")

	//printf("nthreads = %i with vector length = %i\n",nthreads,VEC_LENGTH);


	// Start the parallel loop
#pragma omp parallel private(tid,stride) default(shared) num_threads(nthreads)
	{
		nthreads = omp_get_num_threads();
		//printf("nthreads = %i with vector length = %i\n",nthreads,VEC_LENGTH);
		//nthreads = 1;
		stride = (nptcls+nthreads-1)/nthreads;

		tid = omp_get_thread_num();
		//tid = 0;

//	    auto cpu = sched_getcpu();
//	    std::ostringstream os;
//	        os<<"\nThread "<<omp_get_thread_num()<<" on cpu "<<sched_getcpu()<<std::endl;
//	        std::cout<<os.str()<<std::flush;

		PlasmaData pdata_local = *pdata;

		// Each thread gets a separate copy of the accumulation arrays
		HOMoments* my_moment = moments+tid;

		// Initialize the moment values
		//printf("Initializing moment values\n");
		my_moment->set_vals(0);

		int nSubcycle_max = pdata->nSubcycle_max;

		int ptcl_start,ptcl_end;
		int nptcls_process;
		int nptcls_left;
		int ishrink = 0;
		int nptcl_replacements = 0;

		int nptcl_done;
		//int iptcl_max;
		int iptcl_new_v[VEC_LENGTH];
		int iptcl_v[VEC_LENGTH];
		int iter_array_v[VEC_LENGTH];

		int* iptcl_new = iptcl_new_v;
		int* iptcl = iptcl_v;
		int* iter_array = iter_array_v;

		long long int nSubSteps_done = 0;

		ptcl_start = stride*tid;
		ptcl_end = fmin(stride*(tid+1)-1,nptcls-1);

		nptcls_process = ptcl_end-ptcl_start+1;

		//printf("Thread %i starting at %i to %i with %i ptcls\n",
			//	tid,ptcl_start,ptcl_end,nptcls_process);


		ParticleObjNT<VEC_LENGTH,nSpatial,nVel,iEM> particle(iptcl);

		// Populate the timers
		particle.piccard_timer = piccard_timer+tid;
		particle.accel_timer = accel_timer+tid;
		particle.tally_timer = tally_timer+tid;
		particle.crossing_timer = crossing_timer+tid;
		particle.dtau_est_timer = dtau_est_timer+tid;

//		ParticleObjN<VEC_LENGTH> particle(iptcl);

		typevecN<int,VEC_LENGTH> iter;


		iter = 0;
		for(int i=0;i<VEC_LENGTH;i++)
			iter_array[i] = 0;

		CurrentTallyCPU currents(&my_moment->get_val(0,0,0,ispecies,HOMoments_currentx),
							  &my_moment->get_val(0,0,0,ispecies,HOMoments_currenty),
							  &my_moment->get_val(0,0,0,ispecies,HOMoments_currentz),
							  make_int3(moments->pdata->nx,moments->pdata->ny,moments->pdata->nz),
							  moments->pdata->dxdi,moments->pdata->dydi,moments->pdata->dzdi,
							  moments->pdata->ndimensions);

		ChargeTally charge(&my_moment->get_val(0,0,0,ispecies,HOMoments_charge),
							  make_int3(moments->pdata->nx,moments->pdata->ny,moments->pdata->nz),
							  moments->pdata->dxdi,moments->pdata->dydi,moments->pdata->dzdi,
							  moments->pdata->ndimensions);

		StressTally stress(&moments->get_val(0,0,0,ispecies,HOMoments_S2xx),
				&moments->get_val(0,0,0,ispecies,HOMoments_S2xy),
				&moments->get_val(0,0,0,ispecies,HOMoments_S2xz),
				&moments->get_val(0,0,0,ispecies,HOMoments_S2yy),
				&moments->get_val(0,0,0,ispecies,HOMoments_S2yz),
				&moments->get_val(0,0,0,ispecies,HOMoments_S2zz),
							  moments->pdata->nx,moments->pdata->ny,moments->pdata->nz,
							  moments->pdata->ndimensions,moments->pdata->nVelocity);

		for(int i=0;i<VEC_LENGTH;i++)
			iptcl[i] = ptcl_start+i;

		nptcl_done = 0;


		load_store_timer[tid].start();
		particle = *this;




		//for(int i=0;i<VEC_LENGTH;i++)
		//	particle.dt_finished(i) = 0;

		// Each thread loops over its own particles
		// In order to avoid SIMD divergence we loop until
		// all particles in the threads work que have been
		// pushed. Anytime a particle finishes a subcycle
		// it is written back to the main list and a new particle
		// takes its slot
		while(nptcl_done < nptcls_process)
		{
			nptcls_left = nptcls_process-nptcl_done;

			//printf("nptcls_left = %i, ntpcl_done = %i\n",nptcls_left,nptcl_done);

			if((nptcls_left <= VEC_LENGTH)&&(VEC_LENGTH > 1))
			{
				if(ishrink == 0)
				{
					for(int j=0;j<VEC_LENGTH;j++)
					{
						//printf("iptcl[%i] = %i\n",j,iptcl[0][j]);
						particle.write_back(*this,j);
					}

					int k = 0;
					for(int l=0;l<VEC_LENGTH;l++)
					{

						bool idone = 0;

						//printf("iter2(%i) = %f\n",j,particles2.dt_finished(j));
						if(particle.dt_finished(l) >= pdata->dt)
						{
							idone = 1;
						}
						else if(iter(l) >= pdata->nSubcycle_max)
						{
							idone = 1;
//							printf("warning particle finished before time step was finished dt_left[%i] = %e\n",iptcl[l],pdata->dt-particle.dt_finished(l));
						}
						else if(iptcl[l] > ptcl_end)
							idone = 1;
						else
							idone = 0;


						if(idone)
						{
							nSubSteps_done += iter(l);
							num_subcycles[iptcl[l]] += iter(l);
							iter(l) = 0;

							// Accumulate Charge and S2 moment

						}
						else
						{
							iptcl[k] = iptcl[l];
							iter_array[k] = iter(l);


							k++;
						}
					}

					nptcl_done = nptcls_process - k ;
					nptcls_left = k;

					ishrink = 1;
				}

// Hack to compile all versions of ParticleObjN template
				shrink_pushT<VEC_LENGTH,nSpatial,nVel,iEM>(pdata,fields,&currents,this,
									&iter_array,&iptcl,&iptcl_new,
									nptcls_left,nptcl_done,nptcls_process,nSubSteps_done);
//				shrink_push<VEC_LENGTH>(pdata,fields,&currents,this,
//									&iter_array,&iptcl,&iptcl_new,
//									nptcls_left,nptcl_done,nptcls_process,nSubSteps_done);


			}
			else
			{
//				for(int j=0;j<VEC_LENGTH;j++)
//					printf("particle %i done = %f, %f, %f, %i, %i, %i, %f, %f, %f\n",
//							iptcl[j],particle.px(j),particle.py(j),particle.pz(j),
//							particle.ix(j),particle.iy(j),particle.iz(j),
//							particle.vx(j),particle.vy(j),particle.vz(j));

				// Here our particle vector size is the same
				// size as our system vector size, and won't
				// change from step to step

				particle.push(pdata,fields,&currents,iter,nSubcycle_max);



				// Replace the particle (or particles) that
				// have finished their subcycle steps
				//int k = 0;
				for(int j=0;j<VEC_LENGTH;j++)
				{
					bool idone = 0;

					if(particle.dt_finished(j) >= pdata->dt)
					{
						idone = 1;
					}
					else if(iter(j) >= pdata->nSubcycle_max)
					{
						idone = 1;

//						printf("warning particle finished before time step was finished dt_left[%i] = %e\n",iptcl[j],pdata->dt-particle.dt_finished(j));
					}

					if(idone)
					{
						// Accumulate Charge and S2 moment

//						printf("particle %i done = %f, %f, %f, %i, %i, %i, %f, %f, %f\n",
//								iptcl[j],particle.px(j),particle.py(j),particle.pz(j),
//								particle.ix(j),particle.iy(j),particle.iz(j),
//								particle.vx(j),particle.vy(j),particle.vz(j));

						// Write results, and get a new particle from the list
						particle.write_back(*this,j);

						num_subcycles[iptcl[j]] += iter(j);

						iptcl[j] = ptcl_start + nptcl_done + VEC_LENGTH;
						nptcl_done++;

						if(nptcls_process-nptcl_done > 0)
						{
							particle.copy_in(*this,j);
						}

						nSubSteps_done += iter(j);

						iter(j) = 0;
						particle.dt_finished(j) = 0.0f;


					}
				} /* for(int j=0;j<nptcls_left;j++) */
				//printf("nptcls_left = %i, ntpcl_done = %i\n",nptcls_left,nptcl_done);

			} /* else */

			nptcl_replacements++;

		} /* while(nptcl_done < nptcls_process) */

		load_store_timer[tid].stop();

		tally_timer2[tid].start();
		// accumulate charge and s2 moment
		for(int i=ptcl_start;i<=ptcl_end;i++)
		{
			charge.tally(px[i],py[i],pz[i],
					ix[i],iy[i],iz[i],
					1.0);

			stress.tally1d1v(px[i],
					vx[i],
					ix[i],
					1.0f);


			//if(fabs(dt_finished[i] - pdata->dt) > 1.0e-5)
			//	printf("particle %i dt_finished = %e\n",i,dt_finished[i]);

			dt_finished[i] = 0.0f;

		}
		tally_timer2[tid].stop();

		//nSubSteps_proc[0] = nSubSteps_done;

		nSubSteps_proc[tid] = nSubSteps_done;

//		printf("average particles processed per replacement: %f\n",nptcls_process/((double)nptcl_replacements));



	} /* pragma omp parallel */

	for(int i=1;i<nthreads;i++)
		nSubSteps_proc[0] += nSubSteps_proc[i];

	//printf("nsteps avg = %i\n",nSubSteps_proc[0]);

	return nSubSteps_proc[0];


}
