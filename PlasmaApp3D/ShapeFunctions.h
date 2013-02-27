#ifndef SHAPE_FUNCTIONS_H
#define SHAPE_FUNCTINOS_H

#include "PlasmaData.h"
#include <stdlib.h>
#include <math.h>

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
