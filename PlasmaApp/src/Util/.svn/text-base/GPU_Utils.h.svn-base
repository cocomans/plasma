/*-------------------------------------------------------------------------*/
/**
	@file		GPU_Utils.h
	@author	J. Payne
	@date		3/01/2013
	@brief	Declares functions for a thrust interface and some general GPU utility kernels

	something funky was
	happening with thrust during a device code link, so this gets compiled independently.
	This file also contains some general GPU utility kernals for sorting and data reorder

*/
/*--------------------------------------------------------------------------*/
#ifndef GPU_UTILS_H
#define GPU_UTILS_H

#  define CUDA_SAFE_KERNEL(call) {                                         \
	call;																					\
	cudaDeviceSynchronize();														\
	cudaError err = cudaGetLastError();										\
    if ( cudaSuccess != err) {                                               \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
                exit(EXIT_FAILURE);                                                  \
    } }

#  define CUDA_SAFE_CALL(call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }


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
void ReOrderData32_GPU(void* inData,void* outData,int* keys,int n);


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
void ReOrderData64_GPU(void* inData,void* outData,int* keys,int n);

/*-------------------------------------------------------------------------*/
/**
  	@brief Interface to the thrust::sort_by_key() function

  	@param[in] values device pointer to array of values to be sorted according to keys
  	@param[in] keys device pointer to array of keys to sort
  	@paran[in] n length of array to be reordered


*/
/*--------------------------------------------------------------------------*/
void SortByKey(int* values,int* keys,int n);

/*-------------------------------------------------------------------------*/
/**
  	@brief Generate an array of reorder indexes partitioned according to the
  	true / false value in condition

	@param[in] condition true / false condition by which id's will be partitioned
	@param[in,out] index_out array of partitioned indicies
  	@paran[in] n length of array to be partitioned


*/
/*--------------------------------------------------------------------------*/
int GenPartitionIDs(int* condition,int* index_out,int n);

/*-------------------------------------------------------------------------*/
/**
  	@brief Interface to the thrust::scan() function

	@param[in] vals pointer to values to be scanned
  	@paran[in] n length of array


*/
/*--------------------------------------------------------------------------*/
void GPUscan_int(int* vals,int n);

/*-------------------------------------------------------------------------*/
/**
  	@brief Interface to the thrust::reduce() function

	@param[in] vals values to be reduced
  	@paran[in] n length of array to be reduced


*/
/*--------------------------------------------------------------------------*/
void GPUreduce_int(int* vals,int n);


void check_int_vals(int* array,int ncheck);


#endif /* GPU_UTILS_H */
