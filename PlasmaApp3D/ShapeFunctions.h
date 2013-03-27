#ifndef SHAPE_FUNCTIONS_H
#define SHAPE_FUNCTINOS_H

#include "PlasmaData.h"
#include <stdlib.h>
#include <math.h>

#ifndef CUDA_CODE
#include <immintrin.h>
#endif

#define SIGN_MASK 0x7fffffff
union combo_hack{
	        unsigned int in;
	        float fl;
	};

__device__ __host__ static __inline__ float b0(float x){
  combo_hack flip; flip.fl = x;
  flip.in = flip.in & SIGN_MASK;
  if(flip.fl > 0.5f) return 0.0f;
  return 1.0f;
//  if( x >= -0.5 && x <= 0.5) return 1.0;
//  else                       return 0.0;
}
__device__ __host__ static __inline__  float b1(float x){
  combo_hack flip; flip.fl = x;
  flip.in = flip.in & SIGN_MASK;
  if(flip.fl > 1.0f) return 0.0f;
  return (1.0f - flip.fl);
//  if(x >= -1.0 && x < 0.0)      return 1.0 + x;
//  else if(x >= 0.0 && x <= 1.0) return 1.0 - x;
//  else                          return 0.0;
}
__device__ __host__ static __inline__  float b2(float x){
  combo_hack flip; flip.fl = x;
  flip.in = flip.in & SIGN_MASK;
/*  int offset = (flip.fl > 0.5f) + (flip.fl >= 1.5f);
  switch(offset) {
    case 0:
      return __fmaf_rn(-x, x, 0.75f);
    case 1:
      return __fmaf_rn( 0.5f*flip.fl, (flip.fl - 3.0f), 1.125f);
    case 2:
      return 0.f;
  }
  return 0.f;
*/
  if(flip.fl <= 1.5f) {
    if(flip.fl >0.5f) {
      return ( (0.5f*flip.fl)*(flip.fl - 3.0f)+ 1.125f);
      //return (1.125f + 0.5f*flip.fl*flip.fl + 0.5f*flip.fl*flip.fl(-3.0f)));
    }
    return (-x*x+0.75f);
    //return 0.75f - x*x;
  }
  return 0.0f;
//    if(x >= -1.5 && x < -0.5)     return 0.125*(3.0 + 2.0*x)*(3.0 + 2.0*x);
//    else if(x >= -0.5 && x < 0.5) return 0.75 - x*x;
//    else if(x >= 0.5 && x <= 1.5) return 0.125*(3.0 - 2.0*x)*(3.0 - 2.0*x);
//    else                          return 0.0;
}

template <typename T>
__device__ __host__ static __inline__
int sgn(T val) {
    return (val > T(0)) - (val < T(0));

}

__device__ __host__ static __inline__
realkind S1_shape(realkind x)
{
	// x = x/dx

	realkind result = 0;

	x = fabsf(x);

	result = (1.0f - x)*(x <= 1.0f);

	return result;

}

//static const __m256i sign_mask = _mm256_set1_epi32(0x7fffffff);
//
//static const __m256 float_1 = _mm256_set1_ps(1.0);
//
//static const __m256 float_15 = _mm256_set1_ps(1.5);
//
//static const __m256 float_075 = _mm256_set1_ps(0.75);
//
//static const __m256 float_05 = _mm256_set1_ps(0.5);
//
//static __inline__
//__m256 S1_shape(__m256 x)
//{
//	// abs
//	x = _mm256_and_ps(x,sign_mask);
//	__m256 xa = _mm256_cmp_ps(x,float_1,_CMP_LE_OQ);
//
//	x = _mm256_sub_ps(float_1,x);
//
//	return _mm256_and_ps(x,xa);
//
//}

__device__ __host__ static __inline__
realkind S2_shape(realkind x)
{
	// x = x/dx

	realkind result = 0;

	x = fabsf(x);

	if(x <= 0.5f)
	{
		result = 0.75f - (x*x);
	}
	else if(x <= 1.5f)
	{
		result = 0.5f * (1.5f - x) * (1.5f - x);
	}

	return result;
}

//static __inline__
//__m256 S2_shape(__m256 x)
//{
//	// abs
//	x = _mm256_and_ps(x,sign_mask);
//	__m256 comp1 = _mm256_cmp_ps(x,float_05,_CMP_LE_OQ);
//	__m256 comp2 = _mm256_cmp_ps(x,float_15,_CMP_LE_OQ);
//
//	__m256 x05,x15;
//	// x05 = 0.75 - (x*x)
//	x05 = _mm256_sub_ps(float_075,_mm256_mul_ps(x,x));
//
//	// check if (x <= 0.5)
//	x05 = _mm256_and_ps(x05,comp1);
//
//	// x15 = 1.5-x
//	x15 = _mm256_sub_ps(float_15,x);
//	// x15 = 0.5*(1.5-x)*(1.5-x)
//	x15 = _mm256_mul_ps(float_05,_mm256_mul_ps(x15,x15));
//
//	// check if (x <= 1.5 and x > 0.5)
//	x15 = _mm256_andnot_ps(_mm256_and_ps(x15,comp2),comp1);
//
//	return _mm256_add_ps(x05,x15);
//
//}

__device__ __host__ static __inline__
realkind dS1_shape(realkind x)
{
	// x = x/dx

	float result = 0;

	result = -1.0f*sgn(x)*(fabs(x) <= 1.0f);

	return result;
}

__device__ __host__ static __inline__
float dS2_shape(float x)
{
	// x = x/dx

	float result = 0;


	if(fabs(x) <= 0.5f)
	{
		result = -2.0f*x;
	}
	else if(fabs(x) <= 1.5f)
	{
		result = (1.5f - x);
	}

	return result;
}


#endif /* SHAPE_FUNCTIONS_H */
