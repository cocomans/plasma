#ifndef TYPEVEC_N_H
#define TYPEVEC_N_H


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

#ifdef GPU_CODE
#define CUDA_CODE
#endif

#include "vec_funcs.h"
#ifndef CUDA_CODE
#include <immintrin.h>
#endif

// get AVX intrinsics

#define VEC_LENGTH_MAX 32

#define SIGN_MASK 0x7FFFFFFF

template<typename T, const int N>
class typevecN
{
public:
	__attribute__ ((__aligned__(32))) T values[N];

	/* Member Functions and Overloaded Operators */
	const int getsize(void){return N;}

	// Access operator
	__host__ __device__ __inline__
	T& operator()(const int i){return values[i];}

	__host__ __device__ __inline__
	const T& operator()(const int i)const{return values[i];}

	// Populating a typevecN<T,N> from an array
	__host__ __device__ __inline__
	typevecN<T,N>& operator=(const T* array)
	{

		for(int i=0;i<N;i++)
		{
			this->values[i] = array[i];
		}

		return *this;
	}

	// Populating a typevecN<T,N> from a single value
	__host__ __device__ __inline__
	typevecN<T,N>& operator=(const T value)
	{


		for(int i=0;i<N;i++)
		{
			this->values[i] = value;
		}


		return *this;
	}

	template<typename T2> __host__ __device__ __inline__
	typevecN<T,N>& operator=(const typevecN<T2,N>& a)
	{

		for(int i=0;i<N;i++)
		{
			this->values[i] = (T)(a(i));
		}


		return *this;
	}

	// Operations between two typevecN's
	// Addition
	template<typename T2> __host__ __device__ __inline__
	typevecN<T,N> operator+(const typevecN<T2,N>& b)
	const
	{
		typevecN<T,N> c;

		for(int i=0;i<N;i++)
		{
			c.values[i] = values[i] + b.values[i];
		}


		return c;
	}

	__host__ __device__ __inline__
	typevecN<T,N> operator+(const T b)
	const
	{
		typevecN<T,N> c;

		for(int i=0;i<N;i++)
		{
			c.values[i] = values[i] + b;
		}


		return c;
	}

	// Subtraction
	template<typename T2> __host__ __device__ __inline__
	typevecN<T,N> operator-(const typevecN<T2,N>& b)
	const
	{
		typevecN<T,N> c;


		for(int i=0;i<N;i++)
		{
			c.values[i] = values[i] - b.values[i];
		}

		return c;
	}

	__host__ __device__ __inline__
	typevecN<T,N> operator-(const typevecN<T,N>& b)
	const
	{
		typevecN<T,N> c;


		for(int i=0;i<N;i++)
		{
			c.values[i] = values[i] - b.values[i];
		}

		return c;
	}

	// Subtraction
	template<typename T2> __host__ __device__ __inline__
	typevecN<T,N> operator-(const T2 b)
	const
	{
		typevecN<T,N> c;


		for(int i=0;i<N;i++)
		{
			c.values[i] = values[i] - b;
		}

		return c;
	}

	// Multiplication
	template<typename T2> __host__ __device__ __inline__
	typevecN<T,N> operator*(const typevecN<T2,N>& b)
	{
		typevecN<T,N> c;


		for(int i=0;i<N;i++)
		{
			c.values[i] = values[i] * b.values[i];
		}

		return c;
	}

	// Multiplication
	template<typename T2>__host__ __device__ __inline__
	typevecN<T,N> operator*(const T2& b)
	{
		typevecN<T,N> c;


		for(int i=0;i<N;i++)
		{
			c.values[i] = values[i] * b;
		}

		return c;
	}

	// Division
	template<typename T2> __host__ __device__ __inline__
	typevecN<T,N> operator/(const typevecN<T2,N>& b)
	{
		typevecN<T,N> c;


		for(int i=0;i<N;i++)
		{
			c.values[i] = values[i] / b.values[i];
		}

		return c;
	}

	// Division
	template<typename T2> __host__ __device__ __inline__
	typevecN<T,N> operator/(const T2 b)
	{
		typevecN<T,N> c;


		for(int i=0;i<N;i++)
		{
			c.values[i] = values[i] / b;
		}

		return c;
	}



	// Increment Operators
	__host__ __device__ __inline__
	typevecN<T,N>& operator+=(const typevecN<T,N>& b)
	{
		for(int i=0;i<N;i++)
		{
			values[i] +=  b.values[i];
		}

		return *this;
	}

	__host__ __device__ __inline__
	typevecN<T,N>& operator-=(const typevecN<T,N>& b)
	{
		for(int i=0;i<N;i++)
		{
			values[i] -=  b.values[i];
		}

		return *this;
	}




	__host__ __device__ __inline__
	typevecN<T,N>& operator++(int)
	{
		for(int i=0;i<N;i++)
		{
			values[i]++;
		}

		return *this;
	}

	__host__ __device__ __inline__
	typevecN<T,N>& operator--(int)
	{
		for(int i=0;i<N;i++)
		{
			values[i]--;
		}

		return *this;
	}







};

// Subtraction
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> operator-(const float& a,const typevecN<T,N>& b)
{
	typevecN<T,N> c;


	for(int i=0;i<N;i++)
	{
		c.values[i] = a - b.values[i];
	}

	return c;
}

template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> operator-(const double& a,const typevecN<T,N>& b)
{
	typevecN<T,N> c;


	for(int i=0;i<N;i++)
	{
		c.values[i] = a - b.values[i];
	}

	return c;
}

// Division
template<typename T1,typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> operator/(const T1& a,const typevecN<T,N>& b)
{
	typevecN<T,N> c;


	for(int i=0;i<N;i++)
	{
		c.values[i] = a / b.values[i];
	}

	return c;
}





// Fused Multiply Add
template<typename T, const int N> static __host__ __device__ __inline__
typevecN<T,N> __fmad(const typevecN<T,N>& a, const typevecN<T,N>& b, const typevecN<T,N>& c)
{
	// Return a*b + c
	typevecN<T,N> d;


	for(int i=0;i<N;i++)
	{
		d.values[i] = a.values[i] * b.values[i] + c.values[i];
	}

	return d;
}

// Fused Multiply Add
template<typename T, typename T2,const int N> static __host__ __device__ __inline__
typevecN<T,N> __fmad(const typevecN<T2,N>& a, const typevecN<T,N>& b, const typevecN<T,N>& c)
{
	// Return a*b + c
	typevecN<T,N> d;


	for(int i=0;i<N;i++)
	{
		d.values[i] = a.values[i] * b.values[i] + c.values[i];
	}

	return d;
}

// make_float3
template<const int N> static __host__ __device__ __inline__
typevecN<float3,N> make_float3N(const typevecN<float,N>& x,
								const typevecN<float,N>& y,
								const typevecN<float,N>& z)
{
	typevecN<float3,N> d;


	for(int i=0;i<N;i++)
	{
		d.values[i].x = x.values[i];
		d.values[i].y = y.values[i];
		d.values[i].z = z.values[i];
	}

	return d;
}

// Cross Product
template<const int N> static __host__ __device__ __inline__
typevecN<float3,N> cross_productN(const typevecN<float3,N>& a,
								  const typevecN<float3,N>& b)
{
	// Return c = a x b
	typevecN<float3,N> c;


	for(int i=0;i<N;i++)
	{
		c.values[i] = cross_product(a.values[i],b.values[i]);
	}

	return c;
}

template<const int N> static __host__ __device__ __inline__
typevecN<typevecN<float,N>,3> cross_productN3(const typevecN<float,N>& ax,
											 const typevecN<float,N>& ay,
											 const typevecN<float,N>& az,
								  const typevecN<typevecN<float,N>,3>& b)
{
	// Return c = a x b
	typevecN<typevecN<float,N>,3> c;



	for(int i=0;i<N;i++)
	{
		(c.values[0]).values[i] = ay.values[i]*(b.values[2]).values[i]
		                        - (b.values[1]).values[i] * az.values[i];
	}

	for(int i=0;i<N;i++)
	{
		(c.values[1]).values[i] = az.values[i]*(b.values[0]).values[i]
		                        - (b.values[2]).values[i] * ax.values[i];
	}

	for(int i=0;i<N;i++)
	{
		(c.values[2]).values[i] = ax.values[i]*(b.values[1]).values[i]
		                        - (b.values[0]).values[i] * ay.values[i];
	}


	return c;
}

template<const int N> static __host__ __device__ __inline__
typevecN<typevecN<double,N>,3> cross_productN3(const typevecN<double,N>& ax,
											 const typevecN<double,N>& ay,
											 const typevecN<double,N>& az,
								  const typevecN<typevecN<double,N>,3>& b)
{
	// Return c = a x b
	typevecN<typevecN<double,N>,3> c;



	for(int i=0;i<N;i++)
	{
		(c.values[0]).values[i] = ay.values[i]*(b.values[2]).values[i]
		                        - (b.values[1]).values[i] * az.values[i];
	}

	for(int i=0;i<N;i++)
	{
		(c.values[1]).values[i] = az.values[i]*(b.values[0]).values[i]
		                        - (b.values[2]).values[i] * ax.values[i];
	}

	for(int i=0;i<N;i++)
	{
		(c.values[2]).values[i] = ax.values[i]*(b.values[1]).values[i]
		                        - (b.values[0]).values[i] * ay.values[i];
	}


	return c;
}

template<const int N,const int N2> static __host__ __device__ __inline__
typevecN<typevecN<double,N>,N2> cross_productNV(const typevecN<double,N>& ax,
											 const typevecN<double,N>& ay,
											 const typevecN<double,N>& az,
								  const typevecN<typevecN<double,N>,N2>& b)
{
	// Return c = a x b
	typevecN<typevecN<double,N>,N2> c;


	if(N2 == 3){

	for(int i=0;i<N;i++)
	{
		(c.values[0]).values[i] = ay.values[i]*(b.values[2]).values[i]
		                        - (b.values[1]).values[i] * az.values[i];
	}

	for(int i=0;i<N;i++)
	{
		(c.values[1]).values[i] = az.values[i]*(b.values[0]).values[i]
		                        - (b.values[2]).values[i] * ax.values[i];
	}

	for(int i=0;i<N;i++)
	{
		(c.values[2]).values[i] = ax.values[i]*(b.values[1]).values[i]
		                        - (b.values[0]).values[i] * ay.values[i];
	}
	}

	return c;
}

template<const int N,const int N2> static __host__ __device__ __inline__
typevecN<typevecN<float,N>,N2> cross_productNV(const typevecN<float,N>& ax,
											 const typevecN<float,N>& ay,
											 const typevecN<float,N>& az,
								  const typevecN<typevecN<float,N>,N2>& b)
{
	// Return c = a x b
	typevecN<typevecN<float,N>,N2> c;


	if(N2 == 3){

	for(int i=0;i<N;i++)
	{
		(c.values[0]).values[i] = ay.values[i]*(b.values[2]).values[i]
		                        - (b.values[1]).values[i] * az.values[i];
	}

	for(int i=0;i<N;i++)
	{
		(c.values[1]).values[i] = az.values[i]*(b.values[0]).values[i]
		                        - (b.values[2]).values[i] * ax.values[i];
	}

	for(int i=0;i<N;i++)
	{
		(c.values[2]).values[i] = ax.values[i]*(b.values[1]).values[i]
		                        - (b.values[0]).values[i] * ay.values[i];
	}
	}

	return c;
}

// l1norm
template<const int N> static __host__ __device__ __inline__
typevecN<float,N> l1normN(const typevecN<float3,N>& a)
{
	// Return c = a x b
	typevecN<float,N> c;


	for(int i=0;i<N;i++)
	{
		c.values[i] = fabsf(a.values[i].x) + fabsf(a.values[i].y) + fabsf(a.values[i].z);
	}

	return c;
}

// abs
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> abs(const typevecN<T,N>& a)
{
	typevecN<T,N> c;

#if !(defined CUDA_CODE || defined NO_HAND_VEC)
	if(N > 3)
	{
		 static const __m256d sign_mask = _mm256_set1_pd(-0.0);
		for(int i=0;i<N;i+=4)
		{
			__m256d& av = *(__m256d*)(a.values+i);
			__m256d& cv = *(__m256d*)(c.values+i);
			cv = _mm256_andnot_pd(sign_mask,av);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c.values[i] = fabs(a.values[i]);
#else

	for(int i=0;i<N;i++)
		c.values[i] = fabs(a.values[i]);
#endif

	return c;
}

// sqrt
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> sqrt(const typevecN<T,N>& a)
{
	typevecN<T,N> c;

	for(int i=0;i<N;i++)
		c.values[i] = sqrtf((float)(a.values[i]));

	return c;
}

// max
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> max(const typevecN<T,N>& a,const typevecN<T,N>& b)
{
	typevecN<T,N> c;

#if !(defined CUDA_CODE || defined NO_HAND_VEC)
	if(N > 3)
	{
		for(int i=0;i<N;i+=4)
		{
			__m256d& av = *(__m256d*)(a.values+i);
			__m256d& bv = *(__m256d*)(b.values+i);
			__m256d& cv = *(__m256d*)(c.values+i);
			cv = _mm256_max_pd(av,bv);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c.values[i] = fmax(a.values[i],b.values[i]);
#else

	for(int i=0;i<N;i++)
		c.values[i] = fmax(a.values[i],b.values[i]);
#endif
	return c;
}

// min
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> vmin(const typevecN<T,N>& a,const typevecN<T,N>& b)
{
	typevecN<T,N> c;

#if !(defined CUDA_CODE || defined NO_HAND_VEC)
	if(N > 3)
	{
		for(int i=0;i<N;i+=4)
		{
			__m256d& av = *(__m256d*)(a.values+i);
			__m256d& bv = *(__m256d*)(b.values+i);
			__m256d& cv = *(__m256d*)(c.values+i);
			cv = _mm256_min_pd(av,bv);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c.values[i] = fmin(a.values[i],b.values[i]);
#else

	for(int i=0;i<N;i++)
		c.values[i] = fminf(a.values[i],b.values[i]);
#endif
	return c;
}


// quadratic_equation
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<typevecN<T,N>,2> quadratic(const typevecN<T,N>& a,const typevecN<T,N>& b,const typevecN<T,N>& c)
{
	typevecN<typevecN<T,N>,2> d;

	for(int i=0;i<N;i++)
	{
		T radical = sqrtf(b(i)*b(i) - 4.0f*a(i)*c(i));
		d(0).values[i] = (-b(i)-radical)/(2.0f*a(i));
		d(1).values[i] = (-b(i)+radical)/(2.0f*a(i));
	}


	return d;
}

// floor
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> floor(const typevecN<T,N>& a)
{
	typevecN<T,N> c;

	for(int i=0;i<N;i++)
		c.values[i] = floor((double)(a.values[i]));

	return c;
}

// floor
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<int,N> ifloor(const typevecN<T,N>& a)
{
	typevecN<int,N> c;


#if !(defined CUDA_CODE || defined NO_HAND_VEC)
	if(N > 3)
	{
		for(int i=0;i<N;i+=4)
		{
			__m256d& av = *(__m256d*)(a.values+i);
			__m128i& cv = *(__m128i*)(c.values+i);
			cv = _mm256_cvtpd_epi32(_mm256_floor_pd(av));
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c.values[i] = floor(a.values[i]);
#else

	for(int i=0;i<N;i++)
		c.values[i] = floor(a.values[i]);
#endif

	return c;
}

template<typename T,const int N,const int N2> static __host__ __device__ __inline__
typevecN<typevecN<int,N>,N2> ifloor(const typevecN<typevecN<T,N>,N2>& a)
{
	typevecN<typevecN<int,N>,N2> c;

	for(int j=0;j<N2;j++)
	for(int i=0;i<N;i++)
		c(j)(i) = floor(a(j)(i));

	return c;
}

template<const int N> static __host__ __device__ __inline__
typevecN<int,N> imod(const typevecN<int,N>& a,int i_in)
{
	typevecN<int,N> c;

	i_in-=1;

#if !(defined CUDA_CODE || defined NO_HAND_VEC)
	if(N > 7)
	{
		static const __m256i iv = _mm256_set1_epi32(i_in);
		static const __m256& iv2 = *((__m256*)&iv);
		for(int i=0;i<N;i+=8)
		{
			__m256& av = *(__m256*)(a.values+i);
			__m256& cv = *(__m256*)(c.values+i);
			cv = _mm256_and_ps(av,iv2);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c.values[i] = (a.values[i])&(i_in);
#else

	for(int i=0;i<N;i++)
	{
		c.values[i] = (a.values[i])&(i_in);
	}
#endif
	return c;
}

template <typename T> __host__ __device__ int fastsgn(T val) {
    return (T(0) < val) - (val < T(0));
}


template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> sgn(const typevecN<T,N>& a)
{
	typevecN<T,N> c;

	for(int i=0;i<N;i++)
	{
		c.values[i] = fastsgn(a(i));
	}

	return c;
}

template<typename T,const int i>
class tloop {
    tloop<T,i-1> x;
    typevecN<T,i> vec;
};

template<typename T>
class tloop<T,1> {
	typevecN<T,1> vec;
};



#endif /* TYPEVEC_N_H */
