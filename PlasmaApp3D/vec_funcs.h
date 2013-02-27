#ifndef VEC_FUNCS_H
#define VEC_FUNCS_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>

static __inline__
float3 cross_product(const float3& A, const float3& B);

static __inline__
float dot_product(const float3& A, const float3& B);

static __inline__
float3 operator+(const float3& A, const float3& B);

static __inline__
float3 operator-(const float3& A, const float3& B);

static __inline__
float3 operator*(const float& A, const float3& B);

static __inline__
float3 operator*(const float3& A, const float& B);

static __inline__
float3 operator*(const float3& A, const float3& B);

static __inline__
float3& operator+=(float3& A, const float3& B);

static __inline__
float3 __fmadf3(const float& A, const float3& B, const float3& C);

static __inline__
float3 __fmadf3(const float3& A, const float3& B, const float3& C);

static __inline__
float l1norm(const float3& A);


static __inline__
float3 cross_product(const float3& A, const float3& B)
{
	// return C = A x B
	float3 C;

	C.x = A.y*B.z - B.y*A.z;
	C.y = A.z*B.x - A.x*B.z;
	C.z = A.x*B.y - A.y*B.x;

	return C;
}

static __inline__
float dot_product(const float3& A, const float3& B)
{
	// return C = A dot B
	float C;

	C = A.x*B.x + A.y*B.y + A.z*B.z;

	return C;
}


static __inline__
float3 operator+(const float3& A, const float3& B)
{
	float3 result;

	result.x = A.x + B.x;
	result.y = A.y + B.y;
	result.z = A.z + B.z;

	return result;
}

static __inline__
float3 operator-(const float3& A, const float3& B)
{
	float3 result;

	result.x = A.x - B.x;
	result.y = A.y - B.y;
	result.z = A.z - B.z;

	return result;
}

static __inline__
float3 operator*(const float& A, const float3& B)
{
	float3 result;

	result.x = A * B.x;
	result.y = A * B.y;
	result.z = A * B.z;

	return result;
}

static __inline__
float3 operator*(const float3& A, const float& B)
{
	float3 result;

	result.x = A.x * B;
	result.y = A.y * B;
	result.z = A.z * B;

	return result;
}

static __inline__
float3 operator*(const float3& A, const float3& B)
{
	float3 result;

	result.x = A.x * B.x;
	result.y = A.y * B.y;
	result.z = A.z * B.z;

	return result;
}

static __inline__
float3 operator/(const float3& A, const float& B)
{
	float3 result;

	result.x = A.x / B;
	result.y = A.y / B;
	result.z = A.z / B;

	return result;
}

static __inline__
float3& operator+=(float3& A, const float3& B)
{
	A.x +=  B.x;
	A.y +=  B.y;
	A.z +=  B.z;

	return A;
}

static __inline__
float3 __fmadf3(const float& A, const float3& B, const float3& C)
{
	float3 result;

	result.x = A * B.x + C.x;
	result.y = A * B.y + C.y;
	result.z = A * B.z + C.z;

	return result;
}

static __inline__
float3 __fmadf3(const float3& A, const float3& B, const float3& C)
{
	float3 result;

	result.x = A.x * B.x + C.x;
	result.y = A.y * B.y + C.y;
	result.z = A.z * B.z + C.z;

	return result;
}


static __inline__
float l1norm(const float3& A)
{
	return fabs(A.x) + fabs(A.y) + fabs(A.z);
}





#endif /* VEC_FUNCS_H */
