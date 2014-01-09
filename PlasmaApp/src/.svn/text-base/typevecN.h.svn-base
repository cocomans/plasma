#ifndef TYPEVEC_N_H
#define TYPEVEC_N_H


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "PlasmaData.h"

#ifdef GPU_CODE
#define CUDA_CODE
#endif


#include "vec_funcs.h"
#ifndef CUDA_CODE
#include <immintrin.h>
#endif


// get AVX intrinsics

#define VEC_LENGTH_MAX 32


template<typename T, const int N>
class typevecN
{
public:
	__attribute__ ((__aligned__(32))) T values[N];

	/* Member Functions and Overloaded Operators */
	__host__ __device__ __inline__
	const int getsize(void)const{return N;}

	// Access operator
	__host__ __device__ __inline__
	T& operator()(const int i){return values[i];}

	__host__ __device__ __inline__
	const T& operator()(const int i)const{return values[i];}

	__host__ __device__ __inline__
	const T* getPtr(void)const{return &(values[0]);}

	__host__ __device__ __inline__
	T* getPtr(void){return &(values[0]);}



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
			c(i) = values[i] + b(i);
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
			c(i) = values[i] + b;
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
			c(i) = values[i] - b(i);
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
			c(i) = values[i] - b(i);
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
			c(i) = values[i] - b;
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
			c(i) = values[i] * b(i);
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
			c(i) = values[i] * b;
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
			c(i) = values[i] / b(i);
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
			c(i) = values[i] / b;
		}

		return c;
	}



	// Increment Operators
	__host__ __device__ __inline__
	typevecN<T,N>& operator+=(const typevecN<T,N>& b)
	{
		for(int i=0;i<N;i++)
		{
			values[i] +=  b(i);
		}

		return *this;
	}

	__host__ __device__ __inline__
	typevecN<T,N>& operator-=(const typevecN<T,N>& b)
	{
		for(int i=0;i<N;i++)
		{
			values[i] -=  b(i);
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

//template<typename T>
//class typevecN<T,1>
//{
//public:
//    T values;
//
//	/* Member Functions and Overloaded Operators */
//	__host__ __device__ __inline__
//	const int getsize(void)const{return 1;}
//
//	// Access operator
//	__host__ __device__ __inline__
//	T& operator()(const int i){return values;}
//
//	__host__ __device__ __inline__
//	const T& operator()(const int i)const{return values;}
//
//	__host__ __device__ __inline__
//	const T* getPtr(void)const{return &values;}
//
//	__host__ __device__ __inline__
//	T* getPtr(void){return &values;}
//
//	// Populating a typevecN<T,N> from an array
//	__host__ __device__ __inline__
//	typevecN<T,1>& operator=(const T* array)
//	{
//
//
//			this->values = array[0];
//
//
//		return *this;
//	}
//
//	// Populating a typevecN<T,N> from a single value
//	__host__ __device__ __inline__
//	typevecN<T,1>& operator=(const T value)
//	{
//
//
//	this->values = value;
//
//
//
//		return *this;
//	}
//
//	template<typename T2> __host__ __device__ __inline__
//	typevecN<T,1>& operator=(const typevecN<T2,1>& a)
//	{
//
//
//			this->values = (T)(a(0));
//
//
//
//		return *this;
//	}
//
//	// Operations between two typevecN's
//	// Addition
//	template<typename T2> __host__ __device__ __inline__
//	typevecN<T,1> operator+(const typevecN<T2,1>& b)
//	const
//	{
//		typevecN<T,1> c;
//
//
//			c.values = values + b.values;
//
//
//
//		return c;
//	}
//
//	__host__ __device__ __inline__
//	typevecN<T,1> operator+(const T b)
//	const
//	{
//		typevecN<T,1> c;
//
//			c.values = values + b;
//
//
//
//		return c;
//	}
//
//	// Subtraction
//	template<typename T2> __host__ __device__ __inline__
//	typevecN<T,1> operator-(const typevecN<T2,1>& b)
//	const
//	{
//		typevecN<T,1> c;
//
//
//
//			c.values = values - b.values;
//
//
//		return c;
//	}
//
//	__host__ __device__ __inline__
//	typevecN<T,1> operator-(const typevecN<T,1>& b)
//	const
//	{
//		typevecN<T,1> c;
//
//
//			c.values = values - b.values;
//
//
//		return c;
//	}
//
//	// Subtraction
//	template<typename T2> __host__ __device__ __inline__
//	typevecN<T,1> operator-(const T2 b)
//	const
//	{
//		typevecN<T,1> c;
//
//
//
//			c.values = values - b;
//
//
//		return c;
//	}
//
//	// Multiplication
//	template<typename T2> __host__ __device__ __inline__
//	typevecN<T,1> operator*(const typevecN<T2,1>& b)
//	{
//		typevecN<T,1> c;
//
//
//		c.values = values * b.values;
//
//
//		return c;
//	}
//
//	// Multiplication
//	template<typename T2>__host__ __device__ __inline__
//	typevecN<T,1> operator*(const T2& b)
//	{
//		typevecN<T,1> c;
//
//
//		c.values = values * b;
//
//
//		return c;
//	}
//
//	// Division
//	template<typename T2> __host__ __device__ __inline__
//	typevecN<T,1> operator/(const typevecN<T2,1>& b)
//	{
//		typevecN<T,1> c;
//
//
//
//			c.values = values / b.values;
//
//
//		return c;
//	}
//
//	// Division
//	template<typename T2> __host__ __device__ __inline__
//	typevecN<T,1> operator/(const T2 b)
//	{
//		typevecN<T,1> c;
//
//
//		c.values = values / b;
//
//
//		return c;
//	}
//
//
//
//	// Increment Operators
//	__host__ __device__ __inline__
//	typevecN<T,1>& operator+=(const typevecN<T,1>& b)
//	{
//
//			values +=  b.values;
//
//
//		return *this;
//	}
//
//	__host__ __device__ __inline__
//	typevecN<T,1>& operator-=(const typevecN<T,1>& b)
//	{
//
//			values -=  b.values;
//
//
//		return *this;
//	}
//
//
//
//
//	__host__ __device__ __inline__
//	typevecN<T,1>& operator++(int)
//	{
//
//			values++;
//
//
//		return *this;
//	}
//
//	__host__ __device__ __inline__
//	typevecN<T,1>& operator--(int)
//	{
//
//			values--;
//
//
//		return *this;
//	}
//
//};

// Subtraction
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> operator-(const float& a,const typevecN<T,N>& b)
{
	typevecN<T,N> c;


	for(int i=0;i<N;i++)
	{
		c(i) = a - b(i);
	}

	return c;
}

template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> operator-(const double& a,const typevecN<T,N>& b)
{
	typevecN<T,N> c;


	for(int i=0;i<N;i++)
	{
		c(i) = a - b(i);
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
		c(i) = a / b(i);
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
		d(i) = a(i) * b(i);
	}

	for(int i=0;i<N;i++)
	{
		d(i) += c(i);
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
		d(i) = a(i) * b(i);
	}

	for(int i=0;i<N;i++)
	{
		d(i) += c(i);
	}

	return d;
}

template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> __fmad(const typevecN<T,N>& a, const T& b, const typevecN<T,N>& c)
{
	// Return a*b + c
	typevecN<T,N> d;


	for(int i=0;i<N;i++)
	{
		d(i) = a(i) * b;
	}

	for(int i=0;i<N;i++)
	{
		d(i) += c(i);
	}

	return d;
}

// Fused Multiply Add
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> __fmad(const typevecN<T,N>& a, const T& b, const T& c)
{
	// Return a*b + c
	typevecN<T,N> d;


	for(int i=0;i<N;i++)
	{
		d(i) = a(i) * b;
	}

	for(int i=0;i<N;i++)
	{
		d(i) += c;
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
		d(i).x = x(i);
		d(i).y = y(i);
		d(i).z = z(i);
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
		c(i) = cross_product(a(i),b(i));
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
		(c.values[0])(i) = ay(i)*(b.values[2])(i)
		                        - (b.values[1])(i) * az(i);
	}

	for(int i=0;i<N;i++)
	{
		(c.values[1])(i) = az(i)*(b.values[0])(i)
		                        - (b.values[2])(i) * ax(i);
	}

	for(int i=0;i<N;i++)
	{
		(c.values[2])(i) = ax(i)*(b.values[1])(i)
		                        - (b.values[0])(i) * ay(i);
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
		(c(0))(i) = ay(i)*(b(2))(i)
		                        - (b(1))(i) * az(i);
	}

	for(int i=0;i<N;i++)
	{
		(c(1))(i) = az(i)*(b(0))(i)
		                        - (b(2))(i) * ax(i);
	}

	for(int i=0;i<N;i++)
	{
		(c(2))(i) = ax(i)*(b(1))(i)
		                        - (b(0))(i) * ay(i);
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
		(c(0))(i) = ay(i)*b(2)(i)
		                        - b(1)(i) * az(i);
	}

	for(int i=0;i<N;i++)
	{
		(c(1))(i) = az(i)*b(0)(i)
		                        - b(2)(i) * ax(i);
	}

	for(int i=0;i<N;i++)
	{
		(c(2))(i) = ax(i)*(b(1))(i)
		                        - (b(0))(i) * ay(i);
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
		(c(0))(i) = ay(i)*b(2)(i)
		                        - b(1)(i) * az(i);
	}

	for(int i=0;i<N;i++)
	{
		(c(1))(i) = az(i)*b(0)(i)
		                        - b(2)(i) * ax(i);
	}

	for(int i=0;i<N;i++)
	{
		(c(2))(i) = ax(i)*(b(1))(i)
		                        - (b(0))(i) * ay(i);
	}
	}

	return c;
}

// l1norm
template<const int N> static __host__ __device__ __inline__
typevecN<float,N> l1normN(const typevecN<float3,N>& a)
{
	// Return c = a x b
	register typevecN<float,N> c;


	for(int i=0;i<N;i++)
	{
		c(i) = fabsf(a(i).x) + fabsf(a(i).y) + fabsf(a(i).z);
	}

	return c;
}

// abs
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> abs(const typevecN<T,N>& a)
{
	typevecN<T,N> c;

#if (!(defined CUDA_CODE || defined NO_HAND_VEC))
#ifdef DOUBLE_PRECISION
	if(N > 3)
	{
		 static const __m256d sign_mask = _mm256_set1_pd(-0.0);
		for(int i=0;i<N;i+=4)
		{
			__m256d& av = *(__m256d*)(a.getPtr()+i);
			__m256d& cv = *(__m256d*)(c.getPtr()+i);
			cv = _mm256_andnot_pd(sign_mask,av);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c(i) = fabs(a(i));
#else
	if(N > 7)
	{
		 static const __m256 sign_mask = _mm256_set1_ps(-0.0);
		for(int i=0;i<N;i+=7)
		{
			__m256& av = *(__m256*)(a.getPtr()+i);
			__m256& cv = *(__m256*)(c.getPtr()+i);
			cv = _mm256_andnot_ps(sign_mask,av);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c(i) = fabsf(a(i));
#endif
#else
#ifdef DOUBLE_PRECISION
	for(int i=0;i<N;i++)
		c(i) = fabs(a(i));
#else
	for(int i=0;i<N;i++)
		c(i) = fabsf(a(i));
#endif
#endif

	return c;
}

// sqrt
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> sqrtv(const typevecN<T,N>& a)
{
	typevecN<T,N> c;
#ifdef DOUBLE_PRECISION
	for(int i=0;i<N;i++)
		c(i) = sqrt((a(i)));
#else
	for(int i=0;i<N;i++)
		c(i) = sqrtf((a(i)));
#endif

	return c;
}

// max
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> max(const typevecN<T,N>& a,const typevecN<T,N>& b)
{
	typevecN<T,N> c;

#if (!(defined CUDA_CODE || defined NO_HAND_VEC))
#ifdef DOUBLE_PRECISION
	if(N > 3)
	{
		for(int i=0;i<N;i+=4)
		{
			__m256d& av = *(__m256d*)(a.getPtr()+i);
			__m256d& bv = *(__m256d*)(b.getPtr()+i);
			__m256d& cv = *(__m256d*)(c.getPtr()+i);
			cv = _mm256_max_pd(av,bv);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c(i) = fmax(a(i),b(i));
#else
	if(N > 7)
	{
		for(int i=0;i<N;i+=8)
		{
			__m256& av = *(__m256*)(a.getPtr()+i);
			__m256& bv = *(__m256*)(b.getPtr()+i);
			__m256& cv = *(__m256*)(c.getPtr()+i);
			cv = _mm256_max_ps(av,bv);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c(i) = fmaxf(a(i),b(i));
#endif
#else
#ifdef DOUBLE_PRECISION
	for(int i=0;i<N;i++)
		c(i) = fmax(a(i),b(i));
#else
	for(int i=0;i<N;i++)
		c(i) = fmaxf(a(i),b(i));
#endif
#endif
	return c;
}

// min
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> vmin(const typevecN<T,N>& a,const typevecN<T,N>& b)
{
	typevecN<T,N> c;

#if (!(defined CUDA_CODE || defined NO_HAND_VEC))
#ifdef DOUBLE_PRECISION
	if(N > 3)
	{
		for(int i=0;i<N;i+=4)
		{
			__m256d& av = *(__m256d*)(a.getPtr()+i);
			__m256d& bv = *(__m256d*)(b.getPtr()+i);
			__m256d& cv = *(__m256d*)(c.getPtr()+i);
			cv = _mm256_min_pd(av,bv);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c(i) = fmin(a(i),b(i));
#else
	if(N > 7)
	{
		for(int i=0;i<N;i+=8)
		{
			__m256& av = *(__m256*)(a.getPtr()+i);
			__m256& bv = *(__m256*)(b.getPtr()+i);
			__m256& cv = *(__m256*)(c.getPtr()+i);
			cv = _mm256_min_ps(av,bv);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c(i) = fminf(a(i),b(i));
#endif

#else
#ifdef DOUBLE_PRECISION
	for(int i=0;i<N;i++)
		c(i) = fmin(a(i),b(i));
#else
	for(int i=0;i<N;i++)
		c(i) = fminf(a(i),b(i));
#endif
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
		d(0)(i) = (-b(i)-radical)/(2.0f*a(i));
		d(1)(i) = (-b(i)+radical)/(2.0f*a(i));
	}


	return d;
}

// floor
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<T,N> floor(const typevecN<T,N>& a)
{
	typevecN<T,N> c;

	for(int i=0;i<N;i++)
		c(i) = floor((double)(a(i)));

	return c;
}

// floor
template<typename T,const int N> static __host__ __device__ __inline__
typevecN<int,N> ifloor(const typevecN<T,N>& a)
{
	typevecN<int,N> c;


#if !(defined CUDA_CODE || defined NO_HAND_VEC)
#ifdef DOUBLE_PRECISION
	if(N > 3)
	{
		for(int i=0;i<N;i+=4)
		{
			__m256d& av = *(__m256d*)(a.getPtr()+i);
			__m128i& cv = *(__m128i*)(c.getPtr()+i);
			cv = _mm256_cvtpd_epi32(_mm256_floor_pd(av));
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c(i) = floor(a(i));
#else
	if(N > 7)
	{
		for(int i=0;i<N;i+=8)
		{
			__m256& av = *(__m256*)(a.getPtr()+i);
			__m256i& cv = *(__m256i*)(c.getPtr()+i);
			cv = _mm256_cvtps_epi32(_mm256_floor_ps(av));
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c(i) = floorf(a(i));
#endif
#else
#ifdef DOUBLE_PRECISION
	for(int i=0;i<N;i++)
		c(i) = floor(a(i));
#else
	for(int i=0;i<N;i++)
		c(i) = floorf(a(i));
#endif
#endif

	return c;
}

template<typename T,const int N,const int N2> static __host__ __device__ __inline__
typevecN<typevecN<int,N>,N2> ifloor(const typevecN<typevecN<T,N>,N2>& a)
{
	typevecN<typevecN<int,N>,N2> c;

	for(int j=0;j<N2;j++)
		c(j) = ifloor(a(j));

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
			__m256& av = *(__m256*)(a.getPtr() + i);
			__m256& cv = *(__m256*)(c.getPtr() + i);
			cv = _mm256_and_ps(av,iv2);
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
	for(int i=0;i<N;i++)
		c(i) = (a(i))&(i_in);
#else

	for(int i=0;i<N;i++)
	{
		c(i) = (a(i))&(i_in);
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
		c(i) = fastsgn(a(i));
	}

	return c;
}

template<const int N> static __host__ __device__ __inline__
typevecN<realkind,N> sgn(const typevecN<realkind,N>& a)
{
	typevecN<realkind,N> c;




#if !(defined CUDA_CODE || defined NO_HAND_VEC)
#ifdef DOUBLE_PRECISION
	if(N > 3)
	{
		const __m256d v0 = _mm256_setzero_pd();
		const __m256i one_mask = _mm256_setr_epi64x(0x3FF0000000000000,
				0x3FF0000000000000,
				0x3FF0000000000000,
				0x3FF0000000000000);
		const __m256d one_mask_d = (__m256d)one_mask;

		for(int i=0;i<N;i+=4)
		{
			__m256d& av = *(__m256d*)(a.getPtr()+i);
			__m256d& cv = *(__m256d*)(c.getPtr()+i);
			cv = _mm256_sub_pd(_mm256_and_pd(_mm256_cmp_pd(v0,av,_CMP_LT_OQ),one_mask_d),
						_mm256_and_pd(_mm256_cmp_pd(av,v0,_CMP_LT_OQ),one_mask_d));
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
		for(int i=0;i<N;i++)
		{
			c(i) = fastsgn(a(i));
		}
#else
	if(N > 7)
	{
		const __m256 v0 = _mm256_setzero_ps();
		const __m256i one_mask =  _mm256_setr_epi32(0x3F800000,
				0x3F800000,
				0x3F800000,
				0x3F800000,0x3F800000,0x3F800000,0x3F800000,0x3F800000);
		const __m256 one_mask_d = (__m256)one_mask;

		for(int i=0;i<N;i+=8)
		{
			__m256& av = *(__m256*)(a.getPtr()+i);
			__m256& cv = *(__m256*)(c.getPtr()+i);
			cv = _mm256_sub_ps(_mm256_and_ps(_mm256_cmp_ps(v0,av,_CMP_LT_OQ),one_mask_d),
						_mm256_and_ps(_mm256_cmp_ps(av,v0,_CMP_LT_OQ),one_mask_d));
			//_mm256_store_pd(c.values+i,cv);
		}
	}
	else
		for(int i=0;i<N;i++)
		{
			c(i) = fastsgn(a(i));
		}
#endif
#else

	for(int i=0;i<N;i++)
	{
		c(i) = fastsgn(a(i));
	}

#endif

	return c;
}


//template<typename T,const int i>
//class tloop {
//    tloop<T,i-1> x;
//    typevecN<T,i> vec;
//};
//
//template<typename T>
//class tloop<T,1> {
//	typevecN<T,1> vec;
//};



#endif /* TYPEVEC_N_H */
