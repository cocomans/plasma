/*-------------------------------------------------------------------------*/
/**
  @file		HOMoments.cu
*/
/*--------------------------------------------------------------------------*/
#include "HOMomentsGPUs.h"
#include "PlasmaData.h"
#include "ParallelInfo.h"
#include "RunData.h"
#include "math.h"
#include <omp.h>
#include "HOMomentsCPU.h"

HOMomentsGPUs::HOMomentsGPUs(PlasmaData* pdata_in)
{


	pdata = pdata_in;

	device_type = 0;

	nx = pdata->nx;
	ny = pdata->ny;
	nz = pdata->nz;
	nspecies = pdata->nspecies;

	int ntotal = nx*ny*nz*pdata->nspecies;
	long long int alloc_size = nx*ny*nz*nspecies*10;


	CUDA_SAFE_CALL(cudaMalloc((void**)&all_data,alloc_size*sizeof(realkind)));
//		charge 		= (realkind*)malloc(ntotal*sizeof(realkind));
//
//		currentx 	= (realkind*)malloc(ntotal*sizeof(realkind));
//		currenty 	= (realkind*)malloc(ntotal*sizeof(realkind));
//		currentz 	= (realkind*)malloc(ntotal*sizeof(realkind));
//
//		S2xx = (realkind*)malloc(ntotal*sizeof(realkind));
//		S2xy = (realkind*)malloc(ntotal*sizeof(realkind));
//		S2xz = (realkind*)malloc(ntotal*sizeof(realkind));
//		S2yy = (realkind*)malloc(ntotal*sizeof(realkind));
//		S2yz = (realkind*)malloc(ntotal*sizeof(realkind));
//		S2zz = (realkind*)malloc(ntotal*sizeof(realkind));

	charge = all_data;
	currentx = charge+ntotal;
	currenty = currentx+ntotal;
	currentz = currenty+ntotal;
	S2xx = currentz+ntotal;
	S2xy = S2xx+ntotal;
	S2xz = S2xy+ntotal;
	S2yy = S2xz+ntotal;
	S2yz = S2yy+ntotal;
	S2zz = S2yz+ntotal;

	set_vals(0.0);
}

__host__
void HOMomentsGPUs::set_vals(realkind val_in)
{
	nx = pdata->nx;
	ny = pdata->ny;
	nz = pdata->nz;

	int ntotal = nx*ny*nz*pdata->nspecies;



	CUDA_SAFE_CALL(cudaMemset(charge,val_in,10*ntotal*sizeof(realkind)));



}
//
//__host__
//void HOMomentsGPUs::set_vals(float val_in)
//{
//	int nx = pdata->nx;
//	int ny = pdata->ny;
//	int nz = pdata->nz;
//
//	int ntotal = nx*ny*nz*pdata->nspecies;
//
//
//
//	CUDA_SAFE_CALL(cudaMemset(charge,val_in,10*ntotal*sizeof(realkind)));
//
//
//
//}

void HOMomentsGPUs::apply_weights(void)
{


}

realkind HOMomentsGPUs::evaluate_energy(void)
{

}


void HOMomentsGPUs::copy_from(HOMomentsCPU* src)
{
	// Copy all the moment values from src to this
	int nalloc = pdata->nx*pdata->ny*pdata->nz*pdata->nspecies;

		enum cudaMemcpyKind kind;

		kind = cudaMemcpyHostToDevice;


		CUDA_SAFE_CALL(cudaMemcpy(all_data,src->all_data,10*nalloc*sizeof(realkind),kind));



}


void HOMomentsGPUs::copy_to(HOMomentsCPU* dst)
{
	// Copy all the moment values from src to this
	int nalloc = pdata->nx*pdata->ny*pdata->nz*pdata->nspecies;

		enum cudaMemcpyKind kind;

		kind = cudaMemcpyDeviceToHost;


		CUDA_SAFE_CALL(cudaMemcpy(dst->all_data,all_data,10*nalloc*sizeof(realkind),kind));



}

__attribute__((noinline))
realkind& HOMomentsGPUs::get_val(const int ix, const int iy, const int iz,
		const int ispecies,enum HOMoments_moment moment)
{
	realkind* result;

	int ix2,iy2,iz2;


	ix2 = ((ix%nx)+nx)%nx;
	iy2 = ((iy%ny)+ny)%ny;
	iz2 = ((iz%nz)+nz)%nz;


	int iout = ix2 + nx * (iy2 + ny * (iz2 + nz * ispecies));
	switch(moment)
	{
	case HOMoments_charge:
		result = charge + iout;
		break;
	case HOMoments_currentx:
		result = currentx + iout;
		break;
	case HOMoments_currenty:
		result = currenty + iout;
		break;
	case HOMoments_currentz:
		result = currentz + iout;
		break;
	case HOMoments_S2xx:
		result = S2xx + iout;
		break;
	case HOMoments_S2xy:
		result = S2xy + iout;
		break;
	case HOMoments_S2xz:
		result = S2xz + iout;
		break;
	case HOMoments_S2yy:
		result = S2yy + iout;
		break;
	case HOMoments_S2yz:
		result = S2yz + iout;
		break;
	case HOMoments_S2zz:
		result = S2zz + iout;
		break;
	default:
		break;
	}

	return *result;
}

__attribute__((noinline))
const realkind& HOMomentsGPUs::get_val(const int ix, const int iy, const int iz,
		const int ispecies,enum HOMoments_moment moment)
const
{
	realkind* result;

	int ix2,iy2,iz2;

	ix2 = ((ix%nx)+nx)%nx;
	iy2 = ((iy%ny)+ny)%ny;
	iz2 = ((iz%nz)+nz)%nz;


	int iout = ix2 + nx * (iy2 + ny * (iz2 + nz * ispecies));
	switch(moment)
	{
	case HOMoments_charge:
		result = charge + iout;
		break;
	case HOMoments_currentx:
		result = currentx + iout;
		break;
	case HOMoments_currenty:
		result = currenty + iout;
		break;
	case HOMoments_currentz:
		result = currentz + iout;
		break;
	case HOMoments_S2xx:
		result = S2xx + iout;
		break;
	case HOMoments_S2xy:
		result = S2xy + iout;
		break;
	case HOMoments_S2xz:
		result = S2xz + iout;
		break;
	case HOMoments_S2yy:
		result = S2yy + iout;
		break;
	case HOMoments_S2yz:
		result = S2yz + iout;
		break;
	case HOMoments_S2zz:
		result = S2zz + iout;
		break;
	default:
		break;
	}

	return *result;
}
