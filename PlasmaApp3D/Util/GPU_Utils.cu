/*-------------------------------------------------------------------------*/
/**
	@file		GPU_Utils.cu
	@author	J. Payne
	@date		3/01/2013
	@brief	Defines functions for a thrust interface and some general GPU utility kernels

	something funky was happening with thrust during a device code link, so this gets compiled independently.
	This file also contains some general GPU utility kernals for sorting and data reorder

*/
/*--------------------------------------------------------------------------*/

#include "GPU_Utils.h"
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/transform.h>
#include <thrust/count.h>
#include <thrust/sort.h>



template<typename T> __global__
void reorder_data(T* idata,T* odata,int* id_old, int n)
{
	int idx = threadIdx.x;
	int gidx = idx + blockDim.x*blockIdx.x;

	while(gidx < n)
	{
		int ogidx = id_old[gidx];

		odata[gidx] = idata[ogidx];

		gidx += blockDim.x*gridDim.x;
	}

}

/*-------------------------------------------------------------------------*/
/**
  	@brief Reorder an array of 32-bit words according to keys, all arrays on GPU

  	Example keys[0] = 5 means that the information stored at inData[5]
  	will be copied to outData[0]

  	@param[in] inData device pointer to array of 32-bit words to be reordered according to keys
  	@param[in] outData device pointer to array where re-ordered data will be stored
  	@param[in] keys pointer to array of integers containing old data indexes
  	@paran[in] n length of array to be reordered


*/
/*--------------------------------------------------------------------------*/
void ReOrderData32_GPU(void* inData,void* outData,int* keys,int n)
{
	int blockSize = 512;
	int gridSize = (n + blockSize -1)/blockSize;

	CUDA_SAFE_KERNEL((reorder_data<<<blockSize,gridSize>>>(
			(int*)inData,(int*)outData,keys,n)))
}


/*-------------------------------------------------------------------------*/
/**
  	@brief Reorder an array of 64-bit words according to keys, all arrays on GPU

  	Example keys[0] = 5 means that the information stored at inData[5]
  	will be copied to outData[0]

  	@param[in] inData device pointer to array of 64-bit words to be reordered according to keys
  	@param[in] outData device pointer to array where re-ordered data will be stored
  	@param[in] keys device pointer to array of integers containing old data indexes
  	@paran[in] n length of array to be reordered
*/
/*--------------------------------------------------------------------------*/
void ReOrderData64_GPU(void* inData,void* outData,int* keys,int n)
{
	int blockSize = 512;
	int gridSize = (n + blockSize -1)/blockSize;

	printf("gridsize = %i\n",gridSize);

	CUDA_SAFE_KERNEL((reorder_data<double><<<blockSize,gridSize>>>(
			(double*)inData,(double*)outData,keys,n)))
}

/*-------------------------------------------------------------------------*/
/**
  	@brief Interface to the thrust::sort_by_key() function

  	@param[in] values device pointer to array of values to be sorted according to keys
  	@param[in] keys device pointer to array of keys to sort
  	@paran[in] n length of array to be reordered


*/
/*--------------------------------------------------------------------------*/
void SortByKey(int* values,int* keys,int n)
{
	thrust::device_ptr<int> thrust_keys(keys);
	thrust::device_ptr<int> thrust_values(values);

	thrust::stable_sort_by_key(thrust_keys,thrust_keys+n,thrust_values);
	cudaDeviceSynchronize();
}

/*-------------------------------------------------------------------------*/
/**
  	@brief Interface to the thrust::scan() function

	@param[in] vals pointer to values to be scanned
  	@paran[in] n length of array


*/
/*--------------------------------------------------------------------------*/
void GPUscan_int(int* vals,int n)
{

	thrust::device_ptr<int> thrust_values(vals);

	thrust::inclusive_scan(thrust_values,thrust_values+n,thrust_values);

	cudaDeviceSynchronize();
}


__global__
void GenPartitionIDs_g(int* data_out,int* scan_data,int n_elements)
{
	// Copy data from data_in to a new location in data_out if condition is true
	int idx = threadIdx.x;
	int gidx = blockIdx.x*blockDim.x+idx;

	int scan_total = scan_data[n_elements-1];


	while(gidx < n_elements)
	{
		int oidxm = 0;
		int oidxt = scan_data[gidx];
		int oidx;

		if(gidx > 0)
			oidxm = scan_data[gidx-1];

		if(oidxt != oidxm)
		{
			oidx = oidxt - 1;
		}
		else
		{
			oidx = gidx - oidx + scan_total;
		}

		data_out[oidx] = gidx;

		gidx += blockDim.x*gridDim.x;
	}
}

/*-------------------------------------------------------------------------*/
/**
  	@brief Generate an array of reorder indexes partitioned according to the
  	true / false value in condition

	@param[in] condition true / false condition by which id's will be partitioned
	@param[in,out] index_out array of partitioned indicies
  	@paran[in] n length of array to be partitioned


*/
/*--------------------------------------------------------------------------*/
int GenPartitionIDs(int* condition,int* index_out,int n)
{

	int blockSize = 512;
	int gridSize = (n + 4*blockSize -1)/(4*blockSize);

//	// First we need to do an inclusive scan of the conditions
//	thrust::device_ptr<int> thrust_values(condition);
//
//	thrust::inclusive_scan(thrust_values,thrust_values+n,thrust_values);

	int* temp = (int*)malloc(n*sizeof(int));
	CUDA_SAFE_CALL(cudaMemcpy(temp,condition,n*sizeof(int),cudaMemcpyDeviceToHost))



	for(int i=1;i<n;i++)
		temp[i] += temp[i-1];

	CUDA_SAFE_CALL(cudaMemcpy(condition,temp,n*sizeof(int),cudaMemcpyHostToDevice))

	// Now we launch a kernel that populates indicies according to the condition
	CUDA_SAFE_KERNEL((GenPartitionIDs_g<<<blockSize,gridSize>>>(
			index_out,condition,n)));

	int ntrue;
	CUDA_SAFE_CALL(cudaMemcpy(&ntrue,condition+n-1,sizeof(int),cudaMemcpyDeviceToHost));

	free(temp);

	return n-ntrue;

}


/*-------------------------------------------------------------------------*/
/**
  	@brief Interface to the thrust::reduce() function

	@param[in] vals values to be reduced
  	@paran[in] n length of array to be reduced


*/
/*--------------------------------------------------------------------------*/
void GPUreduce_int(int* vals,int n)
{
	thrust::device_ptr<int> thrust_values(vals);

	thrust::reduce(thrust_values,thrust_values+n);

	cudaDeviceSynchronize();
}

__global__
void check_int_vals_g(int* array,int ncheck)
{
	int gidx = threadIdx.x+blockIdx.x*blockDim.x;


	while(gidx < ncheck)
	{
		printf("check_vals[%i] = %i\n",gidx,array[gidx]);

		gidx += blockDim.x*gridDim.x;
	}
}

void check_int_vals(int* array,int ncheck)
{
	int cudaBlockSize = 256;
	int cudaGridSize = 96;

	CUDA_SAFE_KERNEL((check_int_vals_g<<<cudaBlockSize,cudaGridSize>>>(array,ncheck)));
}

