/*-------------------------------------------------------------------------*/
/**
  @file		FieldDataGPU.cpp
*/
/*--------------------------------------------------------------------------*/
#include "FieldDataGPUSimple.cuh"
#include "FieldDataCPU.h"
#include "ShapeFunctions.h"
#include "vec_funcs.h"
#include "PlasmaData.h"
#include "ParallelInfo.h"
#include "mpi.h"

__host__ __device__ __inline__
int zorder(int ix,int iy,int iz)
{
	// Spread the bits of each index
	// so that there are 2 empty bits between each bit
	ix = (ix | (ix << 16)) & 0x030000FF;
	ix = (ix | (ix <<  8)) & 0x0300F00F;
	ix = (ix | (ix <<  4)) & 0x030C30C3;
	ix = (ix | (ix <<  2)) & 0x09249249;

	iy = (iy | (iy << 16)) & 0x030000FF;
	iy = (iy | (iy <<  8)) & 0x0300F00F;
	iy = (iy | (iy <<  4)) & 0x030C30C3;
	iy = (iy | (iy <<  2)) & 0x09249249;

	iz = (iz | (iz << 16)) & 0x030000FF;
	iz = (iz | (iz <<  8)) & 0x0300F00F;
	iz = (iz | (iz <<  4)) & 0x030C30C3;
	iz = (iz | (iz <<  2)) & 0x09249249;

	return ix | (iy << 1) | (iz << 2);

}

__host__
void FieldDataGPUSimple::allocate(int nx_in, int ny_in, int nz_in, int nspecies_in)
{
	nx = nx_in;
	ny = ny_in;
	nz = nz_in;
	nspecies = nspecies_in;

	alloc_size = (nx*ny*nz);

	CUDA_SAFE_CALL(cudaMalloc((void**)&all_data,6*alloc_size*sizeof(realkind)));

	Ex = all_data;
	Ey = Ex+alloc_size;
	Ez = Ey+alloc_size;
	Bx = Ez+alloc_size;
	By = Bx+alloc_size;
	Bz = By+alloc_size;

	CUDA_SAFE_CALL(cudaMalloc((void**)&data,alloc_size*sizeof(FieldValues)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&q2m,nspecies*sizeof(realkind)));

}

__host__
void FieldDataGPUSimple::allocate(PlasmaData* pdata)
{
	allocate(pdata->nx,pdata->ny,pdata->nz,pdata->nspecies);
	ndimensions = pdata->ndimensions;
}


template<int icomponent>__device__
realkind& FieldDataGPUSimple::getET(int ix, int iy, int iz)
{
	realkind* result;
/*
	if(ix >= nx)
		ix -= nx;
	else if(ix < 0)
		ix += nx;

	if(iy >= ny)
		iy -= ny;
	else if(iy < 0)
		iy += ny;

	if(iz >= nz)
		iz -= nz;
	else if(iz < 0)
		iz += nz;
*/
	ix = (((ix&(nx-1)) + nx)&(nx-1));
	iy = (((iy&(ny-1)) + ny)&(ny-1));
	iz = (((iz&(nz-1)) + nz)&(nz-1));

	if(icomponent == 0)
	{
		result = Ex + ix + (nx)*(iy + ny * iz);
	}
	else if(icomponent == 1)
	{
		result = Ey + ix + (nx)*(iy + (ny) * iz);
	}
	else if(icomponent == 2)
	{

		result = Ez + ix + (nx)*(iy + ny * iz);
	}

	return *result;
}

template<int icomponent>__device__
realkind& FieldDataGPUSimple::getBT(int ix, int iy, int iz)
{
	realkind* result;

	ix = (((ix&(nx-1)) + nx)&(nx-1));
	iy = (((iy&(ny-1)) + ny)&(ny-1));
	iz = (((iz&(nz-1)) + nz)&(nz-1));

	if(icomponent == 0)
	{
		result = Bx + ix + (nx)*(iy + ny * iz);
	}
	else if(icomponent == 1)
	{
		result = By + ix + (nx)*(iy + (ny) * iz);
	}
	else if(icomponent == 2)
	{
		result = Bz + ix + (nx)*(iy + ny * iz);
	}

	return *result;
}

template<int icomponent>__device__
realkind& FieldDataGPUSimple::getETz(int ix, int iy, int iz)
{
	realkind* result;

	ix = (((ix&(nx-1)) + nx)&(nx-1));
	iy = (((iy&(ny-1)) + ny)&(ny-1));
	iz = (((iz&(nz-1)) + nz)&(nz-1));

	int iout = ix + (nx)*(iy + ny * 0);

	if(icomponent == 0)
	{
		result = Ex + iout;
	}
	else if(icomponent == 1)
	{
		result = Ey + iout;
	}
	else if(icomponent == 2)
	{

		result = Ez + iout;
	}

	return *result;
}

template<int icomponent>__device__
realkind& FieldDataGPUSimple::getBTz(int ix, int iy, int iz)
{
	realkind* result;

	ix = (((ix&(nx-1)) + nx)&(nx-1));
	iy = (((iy&(ny-1)) + ny)&(ny-1));
	iz = (((iz&(nz-1)) + nz)&(nz-1));

	int iout = ix + (nx)*(iy + ny * 0);

	if(icomponent == 0)
	{
		result = Bx + iout;
	}
	else if(icomponent == 1)
	{
		result = By + iout;
	}
	else if(icomponent == 2)
	{
		result = Bz + iout;
	}

	return *result;
}

__device__
realkind& FieldDataGPUSimple::getE(int ix,int iy,int iz,int icomponent)
{
	realkind* result;

	switch(icomponent)
	{
	case 0:
		result = &getET<0>(ix,iy,iz);
		break;
	case 1:
		result = &getET<1>(ix,iy,iz);
		break;
	case 2:
		result = &getET<2>(ix,iy,iz);
		break;
	default:
		result = NULL;
		break;
	}
	return *result;
}

__device__
realkind& FieldDataGPUSimple::getB(int ix,int iy,int iz,int icomponent)
{
	realkind* result;

	switch(icomponent)
	{
	case 0:
		result = &getBT<0>(ix,iy,iz);
		break;
	case 1:
		result = &getBT<1>(ix,iy,iz);
		break;
	case 2:
		result = &getBT<2>(ix,iy,iz);
		break;
	default:
		result = NULL;
		break;
	}
	return *result;
}

template<int icomponent, enum FieldData_deriv ideriv>__device__ __attribute__((noinline))
realkind FieldDataGPUSimple::intrpET(realkind x, realkind y, realkind z, int icellx, int icelly, int icellz)
{
	realkind result = 0;
	if(ndimensions == 1)
	{
		switch(icomponent)
		{
		case 0:
			// x component

			for(int i=0;i<2;i++)
			{
				realkind xp, yp, zp;

				realkind Etemp = getET<0>(icellx+i,icelly,icellz);

				xp = i-x;

				//printf("xp = %f, %f, %f\n",S1_shape(xp,dx),S2_shape(yp,dy),S2_shape(zp,dz));

				if(ideriv == FieldData_deriv_f)
				{
					result += Etemp*S1_shape(xp);
				}
				else if(ideriv == FieldData_deriv_dfdx)
				{
					result += Etemp*dS1_shape(xp) / dx;
				}
				else if(ideriv == FieldData_deriv_dfdy)
				{
					result += 0;
				}
				else if(ideriv == FieldData_deriv_dfdz)
				{
					result += 0;
				}

			}

			break;
		case 1:
			// y component
			result = 0;
			break;
		case 2:

			// z component
			result = 0;
			break;
		default:
			break;
		}
	} /* if(ndimensions == 1) */
	else if(ndimensions == 2)
	{
		switch(icomponent)
		{
		case 0:
			// x component
			for(int j=-1;j<2;j++)
			{
				for(int i=0;i<2;i++)
				{
					realkind xp, yp, zp;

					realkind Etemp = getET<0>(icellx+i,icelly+j,icellz);

					xp = i - x;
					yp = j + 0.5 - y;


					//printf("xp = %f, %f, %f\n",S1_shape(xp,dx),S2_shape(yp,dy),S2_shape(zp,dz));

					if(ideriv == FieldData_deriv_f)
					{
						result += Etemp*S1_shape(xp) * S2_shape(yp);
					}
					else if(ideriv == FieldData_deriv_dfdx)
					{
						result += Etemp*dS1_shape(xp) * S2_shape(yp) / dx;
					}
					else if(ideriv == FieldData_deriv_dfdy)
					{
						result += Etemp*S1_shape(xp) * dS2_shape(yp) / dy;
					}
					else if(ideriv == FieldData_deriv_dfdz)
					{
						result += Etemp*S1_shape(xp) * S2_shape(yp) / dz;
					}

				}
			}




			break;
		case 1:
			// y component
			for(int j=0;j<2;j++)
			{
				for(int i=-1;i<2;i++)
				{
					realkind xp, yp, zp;

					realkind Etemp = getET<1>(icellx+i,icelly+j,icellz);

					xp = i + 0.5 - x;
					yp = j - y;

					if(ideriv == FieldData_deriv_f)
					{
						result += Etemp*S2_shape(xp) * S1_shape(yp);
					}
					else if(ideriv == FieldData_deriv_dfdx)
					{
						result += Etemp*dS2_shape(xp) * S1_shape(yp) / dx;
					}
					else if(ideriv == FieldData_deriv_dfdy)
					{
						result += Etemp*S2_shape(xp) * dS1_shape(yp) / dy;
					}
					else if(ideriv == FieldData_deriv_dfdz)
					{
						result += Etemp*S2_shape(xp) * S1_shape(yp) / dz;
					}

				}
			}

			break;
		case 2:

			// z component
			for(int j=-1;j<2;j++)
			{
				for(int i=-1;i<2;i++)
				{
					realkind xp, yp, zp;

					realkind Etemp = getET<2>(icellx+i,icelly+j,icellz);

					xp = i + 0.5 - x;
					yp = j + 0.5 - y;

					if(ideriv == FieldData_deriv_f)
					{
						result += Etemp*S2_shape(xp) * S2_shape(yp);
					}
					else if(ideriv == FieldData_deriv_dfdx)
					{
						result += Etemp*dS2_shape(xp) * S2_shape(yp) / dx;
					}
					else if(ideriv == FieldData_deriv_dfdy)
					{
						result += Etemp*S2_shape(xp) * dS2_shape(yp) / dy;
					}
					else if(ideriv == FieldData_deriv_dfdz)
					{
						result += Etemp*S2_shape(xp) * S2_shape(yp) / dz;
					}
				}
			}

			break;
		default:
			break;
		}
	} /* if(ndimensions == 1) */
	else
	{
		switch(icomponent)
		{
		case 0:
			// x component
			for(int k=-1;k<2;k++)
			{
				for(int j=-1;j<2;j++)
				{
					for(int i=0;i<2;i++)
					{
						realkind xp, yp, zp;

						realkind Etemp = getET<0>(icellx+i,icelly+j,icellz+k);

						xp = i - x;
						yp = j + 0.5 - y;
						zp = k + 0.5 - z;


						//printf("xp = %f, %f, %f\n",S1_shape(xp,dx),S2_shape(yp,dy),S2_shape(zp,dz));

						if(ideriv == FieldData_deriv_f)
						{
							result += Etemp*S1_shape(xp) * S2_shape(yp) * S2_shape(zp);
						}
						else if(ideriv == FieldData_deriv_dfdx)
						{
							result += Etemp*dS1_shape(xp) * S2_shape(yp) * S2_shape(zp) / dx;
						}
						else if(ideriv == FieldData_deriv_dfdy)
						{
							result += Etemp*S1_shape(xp) * dS2_shape(yp) * S2_shape(zp) / dy;
						}
						else if(ideriv == FieldData_deriv_dfdz)
						{
							result += Etemp*S1_shape(xp) * S2_shape(yp) * dS2_shape(zp) / dz;
						}

					}
				}
			}



			break;
		case 1:
			// y component
			for(int k=-1;k<2;k++)
			{
				for(int j=0;j<2;j++)
				{
					for(int i=-1;i<2;i++)
					{
						realkind xp, yp, zp;

						realkind Etemp = getET<1>(icellx+i,icelly+j,icellz+k);

						xp = i + 0.5 - x;
						yp = j - y;
						zp = k + 0.5 - z;

						if(ideriv == FieldData_deriv_f)
						{
							result += Etemp*S2_shape(xp) * S1_shape(yp) * S2_shape(zp);
						}
						else if(ideriv == FieldData_deriv_dfdx)
						{
							result += Etemp*dS2_shape(xp) * S1_shape(yp) * S2_shape(zp) / dx;
						}
						else if(ideriv == FieldData_deriv_dfdy)
						{
							result += Etemp*S2_shape(xp) * dS1_shape(yp) * S2_shape(zp) / dy;
						}
						else if(ideriv == FieldData_deriv_dfdz)
						{
							result += Etemp*S2_shape(xp) * S1_shape(yp) * dS2_shape(zp) / dz;
						}

					}
				}
			}
			break;
		case 2:

			// z component
			for(int k=0;k<2;k++)
			{
				for(int j=-1;j<2;j++)
				{
					for(int i=-1;i<2;i++)
					{
						realkind xp, yp, zp;

						realkind Etemp = getET<2>(icellx+i,icelly+j,icellz+k);

						xp = i + 0.5 - x;
						yp = j + 0.5 - y;
						zp = k - z;

						if(ideriv == FieldData_deriv_f)
						{
							result += Etemp*S2_shape(xp) * S2_shape(yp) * S1_shape(zp);
						}
						else if(ideriv == FieldData_deriv_dfdx)
						{
							result += Etemp*dS2_shape(xp) * S2_shape(yp) * S1_shape(zp) / dx;
						}
						else if(ideriv == FieldData_deriv_dfdy)
						{
							result += Etemp*S2_shape(xp) * dS2_shape(yp) * S1_shape(zp) / dy;
						}
						else if(ideriv == FieldData_deriv_dfdz)
						{
							result += Etemp*S2_shape(xp) * S2_shape(yp) * dS1_shape(zp) / dz;
						}
					}
				}
			}
			break;
		default:
			break;
		}
	} /* if(ndimensions != 1) */

	return result;
}


template<int icomponent, enum FieldData_deriv ideriv>__device__ __attribute__((noinline))
realkind FieldDataGPUSimple::intrpBT(realkind x, realkind y, realkind z, int icellx, int icelly, int icellz)
{

	realkind result = 0;

	// Check Boundary Conditions

	if(ndimensions == 1)
	{
		switch(icomponent)
		{
		case 0:
			// x component

			for(int i=0;i<2;i++)
			{
				realkind xp, yp, zp;

				realkind Etemp = getB(icellx+i,icelly,icellz,0);

				xp = i-x;

				//printf("xp = %f, %f, %f\n",S1_shape(xp,dx),S2_shape(yp,dy),S2_shape(zp,dz));

				if(ideriv == FieldData_deriv_f)
				{
					result += Etemp*S1_shape(xp);
				}
				else if(ideriv == FieldData_deriv_dfdx)
				{
					result += Etemp*dS1_shape(xp) / dx;
				}
				else if(ideriv == FieldData_deriv_dfdy)
				{
					result += 0;
				}
				else if(ideriv == FieldData_deriv_dfdz)
				{
					result += 0;
				}

			}

			break;
		case 1:
			// y component
			for(int i=-1;i<2;i++)
			{
				realkind xp, yp, zp;

				realkind Etemp = getB(icellx+i,icelly,icellz,1);

				xp = i + 0.5 -x;

				//printf("xp = %f, %f, %f\n",S1_shape(xp,dx),S2_shape(yp,dy),S2_shape(zp,dz));

				if(ideriv == FieldData_deriv_f)
				{
					result += Etemp*S1_shape(xp);
				}
				else if(ideriv == FieldData_deriv_dfdx)
				{
					result += Etemp*dS1_shape(xp) / dx;
				}
				else if(ideriv == FieldData_deriv_dfdy)
				{
					result += 0;
				}
				else if(ideriv == FieldData_deriv_dfdz)
				{
					result += 0;
				}

			}
			break;
		case 2:

			// z component
			for(int i=-1;i<2;i++)
			{
				realkind xp, yp, zp;

				realkind Etemp = getB(icellx+i,icelly,icellz,2);

				xp = i + 0.5 -x;

				//printf("xp = %f, %f, %f\n",S1_shape(xp,dx),S2_shape(yp,dy),S2_shape(zp,dz));

				if(ideriv == FieldData_deriv_f)
				{
					result += Etemp*S1_shape(xp);
				}
				else if(ideriv == FieldData_deriv_dfdx)
				{
					result += Etemp*dS1_shape(xp) / dx;
				}
				else if(ideriv == FieldData_deriv_dfdy)
				{
					result += 0;
				}
				else if(ideriv == FieldData_deriv_dfdz)
				{
					result += 0;
				}

			}
			break;
		default:
			break;
		}
	} /* if(ndimensions == 1) */
	else if(ndimensions == 2)
	{

		// Ndimensions == 2 always
//		const realkind s11[4] = {S1_shape(x),S1_shape(1-x),S1_shape(y),S1_shape(1-y)};
//		const realkind s12[4] = {s11[2],s11[3],s11[1],s11[0]};
//
//		const realkind Barr[4] ={getBT<icomponent>(icellx,icelly,icellz),
//			getBT<icomponent>(icellx+1,icelly+1,icellz),
//			getBT<icomponent>(icellx+1,icelly,icellz),
//			getBT<icomponent>(icellx,icelly+1,icellz)};
////			icache = 1;
////		}
//
//
//		__m256d s11v = _mm256_load_pd(s11);
//		__m256d s12v = _mm256_load_pd(s12);
//		__m256d Bv = _mm256_load_pd(Barr);
//
//		__m256d temp1 = _mm256_mul_pd(_mm256_mul_pd(s12v,s11v),Bv);
//
//		temp1 = _mm256_hadd_pd(temp1,temp1);
//
//
//		__m128d hi128 = _mm256_extractf128_pd(temp1,1);
//		__m128d dotproduct = _mm_add_pd(*(__m128d*)&temp1,hi128);
//
//		result = *(double*)&dotproduct;
		switch(icomponent)
		{
		case 0:
			// x component
			for(int j=0;j<2;j++)
			{
				for(int i=0;i<2;i++)
				{
					realkind xp, yp, zp;

					realkind Etemp = getB(icellx+i,icelly+j,icellz,0);

					xp = i - x;
					yp = j - y;


					//printf("xp = %f, %f, %f\n",S1_shape(xp,dx),S2_shape(yp,dy),S2_shape(zp,dz));

					if(ideriv == FieldData_deriv_f)
					{
						result += Etemp*S1_shape(xp) * S1_shape(yp);
					}
					else if(ideriv == FieldData_deriv_dfdx)
					{
						result += Etemp*dS1_shape(xp) * S1_shape(yp) / dx;
					}
					else if(ideriv == FieldData_deriv_dfdy)
					{
						result += Etemp*S1_shape(xp) * dS1_shape(yp) / dy;
					}
					else if(ideriv == FieldData_deriv_dfdz)
					{
						result += Etemp*S1_shape(xp) * S1_shape(yp) / dz;
					}

				}
			}




			break;
		case 1:
			// y component
			for(int j=0;j<2;j++)
			{
				for(int i=0;i<2;i++)
				{
					realkind xp, yp, zp;

					realkind Etemp = getB(icellx+i,icelly+j,icellz,1);

					xp = i - x;
					yp = j - y;

					if(ideriv == FieldData_deriv_f)
					{
						result += Etemp*S1_shape(xp) * S1_shape(yp);
					}
					else if(ideriv == FieldData_deriv_dfdx)
					{
						result += Etemp*dS1_shape(xp) * S1_shape(yp) / dx;
					}
					else if(ideriv == FieldData_deriv_dfdy)
					{
						result += Etemp*S1_shape(xp) * dS1_shape(yp) / dy;
					}
					else if(ideriv == FieldData_deriv_dfdz)
					{
						result += Etemp*S1_shape(xp) * S1_shape(yp) / dz;
					}

				}
			}

			break;
		case 2:

			// z component
			for(int j=0;j<2;j++)
			{
				for(int i=0;i<2;i++)
				{
					realkind xp, yp, zp;

					realkind Etemp = getB(icellx+i,icelly+j,icellz,2);

					xp = i - x;
					yp = j - y;

					if(ideriv == FieldData_deriv_f)
					{
						result += Etemp*S1_shape(xp) * S1_shape(yp);
					}
					else if(ideriv == FieldData_deriv_dfdx)
					{
						result += Etemp*dS1_shape(xp) * S1_shape(yp) / dx;
					}
					else if(ideriv == FieldData_deriv_dfdy)
					{
						result += Etemp*S1_shape(xp) * dS1_shape(yp) / dy;
					}
					else if(ideriv == FieldData_deriv_dfdz)
					{
						result += Etemp*S1_shape(xp) * S1_shape(yp) / dz;
					}
				}
			}

			break;
		default:
			break;
		}
	} /* if(ndimensions == 1) */
	else
	{

		switch(icomponent)
		{
		case 0:
			// x component
			for(int k=-1;k<2;k++)
			{
				for(int j=-1;j<2;j++)
				{
					for(int i=0;i<2;i++)
					{
						realkind xp, yp, zp;

						realkind Etemp = getBT<0>(icellx+i,icelly+j,icellz+k);

						xp = i - x;
						yp = j  - y;
						zp = k  - z;


						//printf("xp = %f, %f, %f\n",S1_shape(xp,dx),S2_shape(yp,dy),S2_shape(zp,dz));

						if(ideriv == FieldData_deriv_f)
						{
							result += Etemp*S1_shape(xp) * S2_shape(yp) * S2_shape(zp);
						}
						else if(ideriv == FieldData_deriv_dfdx)
						{
							result += Etemp*dS1_shape(xp) * S2_shape(yp) * S2_shape(zp) / dx;
						}
						else if(ideriv == FieldData_deriv_dfdy)
						{
							result += Etemp*S1_shape(xp) * dS2_shape(yp) * S2_shape(zp) / dy;
						}
						else if(ideriv == FieldData_deriv_dfdz)
						{
							result += Etemp*S1_shape(xp) * S2_shape(yp) * dS2_shape(zp) / dz;
						}

					}
				}
			}



			break;
		case 1:
			// y component
			for(int k=-1;k<2;k++)
			{
				for(int j=0;j<2;j++)
				{
					for(int i=-1;i<2;i++)
					{
						realkind xp, yp, zp;

						realkind Etemp = getBT<1>(icellx+i,icelly+j,icellz+k);

						xp = i + 0.5 - x;
						yp = j - y;
						zp = k + 0.5 - z;

						if(ideriv == FieldData_deriv_f)
						{
							result += Etemp*S2_shape(xp) * S1_shape(yp) * S2_shape(zp);
						}
						else if(ideriv == FieldData_deriv_dfdx)
						{
							result += Etemp*dS2_shape(xp) * S1_shape(yp) * S2_shape(zp) / dx;
						}
						else if(ideriv == FieldData_deriv_dfdy)
						{
							result += Etemp*S2_shape(xp) * dS1_shape(yp) * S2_shape(zp) / dy;
						}
						else if(ideriv == FieldData_deriv_dfdz)
						{
							result += Etemp*S2_shape(xp) * S1_shape(yp) * dS2_shape(zp) / dz;
						}

					}
				}
			}
			break;
		case 2:

			// z component
			for(int k=0;k<2;k++)
			{
				for(int j=-1;j<2;j++)
				{
					for(int i=-1;i<2;i++)
					{
						realkind xp, yp, zp;

						realkind Etemp = getBT<2>(icellx+i,icelly+j,icellz+k);

						xp = i + 0.5 - x;
						yp = j + 0.5 - y;
						zp = k - z;

						if(ideriv == FieldData_deriv_f)
						{
							result += Etemp*S2_shape(xp) * S2_shape(yp) * S1_shape(zp);
						}
						else if(ideriv == FieldData_deriv_dfdx)
						{
							result += Etemp*dS2_shape(xp) * S2_shape(yp) * S1_shape(zp) / dx;
						}
						else if(ideriv == FieldData_deriv_dfdy)
						{
							result += Etemp*S2_shape(xp) * dS2_shape(yp) * S1_shape(zp) / dy;
						}
						else if(ideriv == FieldData_deriv_dfdz)
						{
							result += Etemp*S2_shape(xp) * S2_shape(yp) * dS1_shape(zp) / dz;
						}
					}
				}
			}
			break;
		default:
			break;
		}
	} /* if(ndimensions != 1) */

	return result;
}



__device__
realkind FieldDataGPUSimple::intrpE(realkind x, realkind y, realkind z,
			int icellx, int icelly, int icellz,
			const int icomponent, const enum FieldData_deriv ideriv)
{

	realkind result;

	switch(ideriv)
	{
	case FieldData_deriv_f:
		switch(icomponent)
		{
		case 0:
			result = intrpET<0,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
			break;
		case 1:
			result = intrpET<1,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
			break;
		case 2:
			result = intrpET<2,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
			break;
		default:
			result = 0;
			break;
		}
		break;

	case FieldData_deriv_dfdx:
		switch(icomponent)
		{
		case 0:
			result = intrpET<0,FieldData_deriv_dfdx>(x,y,z,icellx,icelly,icellz);
			break;
		case 1:
			result = intrpET<1,FieldData_deriv_dfdx>(x,y,z,icellx,icelly,icellz);
			break;
		case 2:
			result = intrpET<2,FieldData_deriv_dfdx>(x,y,z,icellx,icelly,icellz);
			break;
		default:
			result = 0;
			break;
		}
		break;

	case FieldData_deriv_dfdy:
		switch(icomponent)
		{
		case 0:
			result = intrpET<0,FieldData_deriv_dfdy>(x,y,z,icellx,icelly,icellz);
			break;
		case 1:
			result = intrpET<1,FieldData_deriv_dfdy>(x,y,z,icellx,icelly,icellz);
			break;
		case 2:
			result = intrpET<2,FieldData_deriv_dfdy>(x,y,z,icellx,icelly,icellz);
			break;
		default:
			result = 0;
			break;
		}
		break;

	case FieldData_deriv_dfdz:
		switch(icomponent)
		{
		case 0:
			result = intrpET<0,FieldData_deriv_dfdz>(x,y,z,icellx,icelly,icellz);
			break;
		case 1:
			result = intrpET<1,FieldData_deriv_dfdz>(x,y,z,icellx,icelly,icellz);
			break;
		case 2:
			result = intrpET<2,FieldData_deriv_dfdz>(x,y,z,icellx,icelly,icellz);
			break;
		default:
			result = 0;
			break;
		}
		break;

	default:
		result = 0;
		break;

	}

	return result;

}

__device__
realkind FieldDataGPUSimple::intrpB(realkind x, realkind y, realkind z,
			int icellx, int icelly, int icellz,
			const int icomponent, const enum FieldData_deriv ideriv)
{

	realkind result;

	switch(ideriv)
	{
	case FieldData_deriv_f:
		switch(icomponent)
		{
		case 0:
			result = intrpBT<0,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
			break;
		case 1:
			result = intrpBT<1,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
			break;
		case 2:
			result = intrpBT<2,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
			break;
		default:
			result = 0;
			break;
		}
		break;

	case FieldData_deriv_dfdx:
		switch(icomponent)
		{
		case 0:
			result = intrpBT<0,FieldData_deriv_dfdx>(x,y,z,icellx,icelly,icellz);
			break;
		case 1:
			result = intrpBT<1,FieldData_deriv_dfdx>(x,y,z,icellx,icelly,icellz);
			break;
		case 2:
			result = intrpBT<2,FieldData_deriv_dfdx>(x,y,z,icellx,icelly,icellz);
			break;
		default:
			result = 0;
			break;
		}
		break;

	case FieldData_deriv_dfdy:
		switch(icomponent)
		{
		case 0:
			result = intrpBT<0,FieldData_deriv_dfdy>(x,y,z,icellx,icelly,icellz);
			break;
		case 1:
			result = intrpBT<1,FieldData_deriv_dfdy>(x,y,z,icellx,icelly,icellz);
			break;
		case 2:
			result = intrpBT<2,FieldData_deriv_dfdy>(x,y,z,icellx,icelly,icellz);
			break;
		default:
			result = 0;
			break;
		}
		break;

	case FieldData_deriv_dfdz:
		switch(icomponent)
		{
		case 0:
			result = intrpBT<0,FieldData_deriv_dfdz>(x,y,z,icellx,icelly,icellz);
			break;
		case 1:
			result = intrpBT<1,FieldData_deriv_dfdz>(x,y,z,icellx,icelly,icellz);
			break;
		case 2:
			result = intrpBT<2,FieldData_deriv_dfdz>(x,y,z,icellx,icelly,icellz);
			break;
		default:
			result = 0;
			break;
		}
		break;

	default:
		result = 0;
		break;

	}

	return result;

}


__device__
void FieldDataGPUSimple::intrpAccel1D1V(realkind x,
				realkind vx,
				int icellx,
/* Outputs */		realkind &accelx)
{
	accelx = 0;

	for(int i=0;i<2;i++)
	{
		realkind xp, yp, zp;

		realkind Etemp = getET<0>(icellx+i,0,0);

		xp = i-x;


		accelx += Etemp*S1_shape(xp);


	}



}

__device__
void FieldDataGPUSimple::intrpAccel(realkind x, realkind y, realkind z,
				realkind vx,realkind vy, realkind vz,
				int icellx, int icelly, int icellz,
/* Outputs */		realkind &accelx,realkind& accely, realkind& accelz)
{
	realkind accelout;

	realkind Bxs,Bys,Bzs;
	realkind Exs,Eys,Ezs;

	Exs = intrpET<0,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
	Eys = intrpET<1,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
	Ezs = intrpET<2,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);

	Bxs = intrpBT<0,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
	Bys = intrpBT<1,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);
	Bzs = intrpBT<2,FieldData_deriv_f>(x,y,z,icellx,icelly,icellz);

	accely = vz*Bxs - vx*Bzs + Eys;
	accelz = vx*Bys - vy*Bxs + Ezs;
	accelx = vy*Bzs - vz*Bys + Exs;



}

void FieldDataGPUSimple::intrpFields(realkind x, realkind y, realkind z,
					int icellx, int icelly, int icellz,
/* Outputs */		realkind &Ex,realkind& Ey, realkind& Ez,
					realkind &Bx,realkind& By, realkind& Bz)
{

		int ixm = ((((icellx-1)&(nx-1)) + nx)&(nx-1));
		int iym = ((((icelly-1)&(ny-1)) + ny)&(ny-1));
		int ix = ((((icellx)&(nx-1)) + nx)&(nx-1));
		int iy = ((((icelly)&(ny-1)) + ny)&(ny-1));
		int ixp = ((((icellx+1)&(nx-1)) + nx)&(nx-1));
		int iyp = ((((icelly+1)&(ny-1)) + ny)&(ny-1));
		int iz = 0;

		realkind result = 0;



		// Ndimensions == 2 always
		realkind s11[4] = {S1_shape(x),S1_shape(1.0f-x),S1_shape(y),S1_shape(1.0f-y)};


		realkind s11ex[6] = {S2_shape(0.5f+y),s11[1],S2_shape(0.5f-y),s11[1],S2_shape(1.5f-y),s11[1]};

		realkind s11ey[6] = {s11[2],s11[2],s11[2],S2_shape(0.5f+x),S2_shape(0.5f-x),S2_shape(1.5f-x)};



		realkind s11ez[9] = {s11ex[0],s11ex[0],s11ex[0],
							s11ex[2],s11ex[2],s11ex[2],
							s11ex[4],s11ex[4],s11ex[4]};

		realkind s12ex[6] = {s11[0],s11ex[0],
							 s11[0],s11ex[2],
							 s11[0],s11ex[4]};

		realkind s12ey[6] = {s11ey[3],s11ey[4],s11ey[5],
							 s11[2],s11[2],s11[2]};

		realkind s12ez[9] = {s11ey[2],s11ey[3],s11ey[4],
							 s11ey[2],s11ey[3],s11ey[4],
							 s11ey[2],s11ey[3],s11ey[4]};

		realkind s12[4] = {s11[2],s11[3],s11[1],s11[0]};


		realkind Bxt[4];

		realkind Byt[4];

		realkind Bzt[4];

		realkind Ext[6];

		realkind Eyt[6];

		realkind Ezt[9];

		// ixm, iy
		Eyt[0] = getET<1>(ixm,iy,iz);
		Ezt[3] = getET<2>(ixm,iy,iz);

		// ix, iy
		Ext[2] = getET<0>(ix,iy,iz);
		Eyt[1] = getET<1>(ix,iy,iz);
		Ezt[4] = getET<2>(ix,iy,iz);
		Bxt[0] = getBT<0>(ix,iy,iz);
		Byt[0] = getBT<1>(ix,iy,iz);
		Bzt[0] = getBT<2>(ix,iy,iz);

		// ixp, iy
		Ext[3] = getET<0>(ixp,iy,iz);
		Eyt[2] = getET<1>(ixp,iy,iz);
		Ezt[5] = getET<2>(ixp,iy,iz);
		Bxt[2] = getBT<0>(ixp,iy,iz);
		Byt[2] = getBT<1>(ixp,iy,iz);
		Bzt[2] = getBT<2>(ixp,iy,iz);

		// ixm, iyp
		Eyt[3] = getET<1>(ixm,iyp,iz);
		Ezt[6] = getET<2>(ixm,iyp,iz);

		// ix, iyp
		Ext[4] = getET<0>(ix,iyp,iz);
		Eyt[4] = getET<1>(ix,iyp,iz);
		Ezt[7] = getET<2>(ix,iyp,iz);
		Bxt[3] = getBT<0>(ix,iyp,iz);
		Byt[3] = getBT<1>(ix,iyp,iz);
		Bzt[3] = getBT<2>(ix,iyp,iz);

		// ixp, iyp
		Ext[5] = getET<0>(ixp,iyp,iz);
		Eyt[5] = getET<1>(ixp,iyp,iz);
		Ezt[8] = getET<2>(ixp,iyp,iz);
		Bxt[1] = getBT<0>(ixp,iyp,iz);
		Byt[1] = getBT<1>(ixp,iyp,iz);
		Bzt[1] = getBT<2>(ixp,iyp,iz);

		// ixm, iym
		Ezt[0] = getET<2>(ixm,iym,iz);

		// ix, iym
		Ext[0] = getET<0>(ix,iym,iz);
		Ezt[1] = getET<2>(ix,iym,iz);

		// ixp, iym
		Ext[1] = getETz<0>(ixp,iym,iz);
		Ezt[2] = getETz<2>(ixp,iym,iz);


		realkind Bxs = 0;
		realkind Bys = 0;
		realkind Bzs = 0;
		realkind Exs = 0;
		realkind Eys = 0;
		realkind Ezs = 0;

		for(int i=0;i<4;i++)
			Bxs += Bxt[i]*s11[i]*s12[i];
		for(int i=0;i<4;i++)
			Bys += Byt[i]*s11[i]*s12[i];
		for(int i=0;i<4;i++)
			Bzs += Bzt[i]*s11[i]*s12[i];


		for(int i=0;i<6;i++)
			Exs += Ext[i]*s11ex[i]*s12ex[i];
		for(int i=0;i<6;i++)
			Eys += Eyt[i]*s11ey[i]*s12ey[i];


		for(int i=0;i<9;i++)
			Ezs += Ezt[i]*s11ez[i]*s12ez[i];



		Ex = Exs;
		Ey = Eys;
		Ez = Ezs;
		Bx = Bxs;
		By = Bys;
		Bz = Bzs;




}

__device__
void FieldDataGPUSimple::intrpFields1D3V(realkind x, realkind y, realkind z,
					int icellx, int icelly, int icellz,
/* Outputs */		realkind &Ex,realkind& Ey, realkind& Ez,
					realkind &Bx,realkind& By, realkind& Bz)
{



	int ixm = ((((icellx-1)&(nx-1)) + nx)&(nx-1));
	int ix = ((((icellx)&(nx-1)) + nx)&(nx-1));
	int ixp = ((((icellx+1)&(nx-1)) + nx)&(nx-1));

		// Ndimensions == 2 always
		realkind s1x[2] = {S1_shape(x),S1_shape(1.0-x)};

//		realkind s12x[3] = {S1_shape(0.5+x),S1_shape(0.5-x),S1_shape(1.5-x)};

		realkind s2x[3] = {S2_shape(0.5+x),S2_shape(0.5-x),S2_shape(1.5-x)};


		realkind Bxs = getBTz<0>(ix,0,0);;
		realkind Bys = 0;
		realkind Bzs = 0;
		realkind Exs = 0;
		realkind Eys = 0;
		realkind Ezs = 0;



		realkind Ext[2];

		realkind Eyt[3];

		realkind Ezt[3];


		realkind Byt[2];
		realkind Bzt[2];

//		Byt[0] = getBTz<1>(ixm,0,0);
//		Bzt[0] = getBTz<2>(ixm,0,0);

		Byt[0] = getBTz<1>(ix,0,0);
		Bzt[0] = getBTz<2>(ix,0,0);

		Byt[1] = getBTz<1>(ixp,0,0);
		Bzt[1] = getBTz<2>(ixp,0,0);

//		Ayt[0] = getATz<1>(ixm,0,0);
//		Azt[0] = getATz<2>(ixm,0,0);
		Eyt[0] = getETz<1>(ixm,0,0);
		Ezt[0] = getETz<2>(ixm,0,0);

//		Ayt[1] = getATz<1>(ix,0,0);
//		Azt[1] = getATz<2>(ix,0,0);
		Ext[0] = getETz<0>(ix,0,0);
		Eyt[1] = getETz<1>(ix,0,0);
		Ezt[1] = getETz<2>(ix,0,0);

//		Ayt[2] = getATz<1>(ixp,0,0);
//		Azt[2] = getATz<2>(ixp,0,0);
		Ext[1] = getETz<0>(ixp,0,0);
		Eyt[2] = getETz<1>(ixp,0,0);
		Ezt[2] = getETz<2>(ixp,0,0);



//		for(int i=0;i<2;i++)
//			Bzs += s1x[i]*(Ayt[i+1]-Ayt[i]);
//		for(int i=0;i<2;i++)
//			Bys += s1x[i]*(Ayt[i+1]-Ayt[i]);


		for(int i=0;i<2;i++)
			Bys += Byt[i]*s1x[i];
		for(int i=0;i<2;i++)
			Bzs += Bzt[i]*s1x[i];


		for(int i=0;i<2;i++)
			Exs += Ext[i]*s1x[i];

		for(int i=0;i<3;i++)
			Eys += Eyt[i]*s2x[i];

		for(int i=0;i<3;i++)
			Ezs += Ezt[i]*s2x[i];



		Ex = Exs;
		Ey = Eys;
		Ez = Ezs;
		Bx = Bxs;
		By = Bys;
		Bz = Bzs;

//		double xt = (x+ix)*pdata->dxdi + pdata->xmin;
//		double yt = (y+iy)*pdata->dydi + pdata->ymin;
//
//
//		printf("Ex: %e , Ey: %e = %e, Ez: %e, Bx: %e, By: %e, Bz: %e\n",Ex,Ey,xt,Ez,Bx,By,Bz);
//




}


__global__
void transpose_aos(realkind* Ex,realkind* Ey, realkind* Ez,
				   realkind* Bx,realkind* By,realkind* Bz,
				   FieldValues* data_in,int alloc_size)
{
	int tidx = threadIdx.x;
	int gidx = tidx + blockIdx.x*blockDim.x;

	FieldValues local_data;

	while(gidx < alloc_size)
	{
		local_data = data_in[gidx];

		Ex[gidx] = local_data.vals[0];
		Ey[gidx] = local_data.vals[1];
		Ez[gidx] = local_data.vals[2];

		Bx[gidx] = local_data.vals[3];
		By[gidx] = local_data.vals[4];
		Bz[gidx] = local_data.vals[5];

		gidx += blockDim.x*gridDim.x;
	}
}

void FieldDataGPUSimple::copy_from(FieldDataCPU* src)
{

	// If field data is stored as a struct of arrays
	if(src->FieldType == FieldData_cpu)
	{
		printf("Copying soa fields to gpu\n");
		CUDA_SAFE_CALL(cudaMemcpy(Bx,src->Bx,alloc_size*sizeof(realkind),cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(By,src->By,alloc_size*sizeof(realkind),cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(Bz,src->Bz,alloc_size*sizeof(realkind),cudaMemcpyHostToDevice));

		CUDA_SAFE_CALL(cudaMemcpy(Ex,src->Ex,alloc_size*sizeof(realkind),cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(Ey,src->Ey,alloc_size*sizeof(realkind),cudaMemcpyHostToDevice));
		CUDA_SAFE_CALL(cudaMemcpy(Ez,src->Ez,alloc_size*sizeof(realkind),cudaMemcpyHostToDevice));


	}
	else if(src->FieldType == FieldData_cpu_aos)
	{
		// Field is stored as an array of structs
		CUDA_SAFE_CALL(cudaMemcpy(data,src->data,alloc_size*sizeof(FieldValues),cudaMemcpyHostToDevice));

		int cudaBlockSize = 512;
		int cudaGridSize = 96;
		// Transpose it
		CUDA_SAFE_KERNEL((transpose_aos<<<cudaBlockSize,cudaGridSize>>>(
				Ex,Ey,Ez,Bx,By,Bz,data,alloc_size)))

	}

	CUDA_SAFE_CALL(cudaMemcpy(q2m,src->q2m,nspecies*sizeof(realkind),cudaMemcpyHostToDevice));
}












