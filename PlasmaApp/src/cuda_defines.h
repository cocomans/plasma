/*
 * cuda_defines.h
 *
 *  Created on: May 20, 2013
 *      Author: payne
 */

#ifndef CUDA_DEFINES_H_
#define CUDA_DEFINES_H_

#ifndef NO_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#else

#define __device__

#define __host__

typedef struct{
	int x,y,z;
} int3;

typedef struct{
	double x,y;
} double2;

static __inline__ __host__ __device__ int3 make_int3(int x, int y, int z)
{
	int3 result;
	result.x = x;
	result.y = y;
	result.z = z;
	return result;
}

class float3
{
public:
	float x,y,z;

	operator int3()
	{

		return make_int3(x,y,z);
	}
};



typedef struct{
	double x,y,z;
} double3;

typedef struct{
	double w,x,y,z;
} double4;

static __inline__ __host__ __device__ float3 make_float3(float x, float y, float z)
{
	float3 result;
	result.x = x;
	result.y = y;
	result.z = z;
	return result;
}

static __inline__ __host__ __device__ double4 make_double4(double x, double y, double z, double w)
{
  double4 t; t.x = x; t.y = y; t.z = z; t.w = w; return t;
}







#endif



#endif /* CUDA_DEFINES_H_ */
