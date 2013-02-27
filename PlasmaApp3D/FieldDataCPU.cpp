/*-------------------------------------------------------------------------*/
/**
  @file		FieldDataCPU.cpp
*/
/*--------------------------------------------------------------------------*/
#include "FieldDataCPU.h"
#include "ShapeFunctions.h"
#include "vec_funcs.h"
#include "PlasmaData.h"
#include "ParallelInfo.h"
#include "mpi.h"
#include <immintrin.h>
#include <omp.h>


__inline__
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

__inline__
int part1by1(int ix)
{
	ix = (ix | (ix << 16)) & 0x030000FF;
	ix = (ix | (ix <<  8)) & 0x0300F00F;
	ix = (ix | (ix <<  4)) & 0x030C30C3;
	ix = (ix | (ix <<  2)) & 0x09249249;

    return ix;
}

__inline__
int zorder2D(int ix,int iy)
{


	return part1by1(ix) | (part1by1(iy) << 1);

}

void FieldDataCPU::allocate(int nx_in, int ny_in, int nz_in, int nspecies_in)
{
	nx = nx_in;
	ny = ny_in;
	nz = nz_in;
	nspecies = nspecies_in;

	alloc_size = zorder(nx,ny,nz);

	alloc_size = nx*ny*nz;

	Ex = (realkind*)malloc(alloc_size*sizeof(realkind));
	Ey = (realkind*)malloc(alloc_size*sizeof(realkind));
	Ez = (realkind*)malloc(alloc_size*sizeof(realkind));

	Bx = (realkind*)malloc(alloc_size*sizeof(realkind));
	By = (realkind*)malloc(alloc_size*sizeof(realkind));
	Bz = (realkind*)malloc(alloc_size*sizeof(realkind));

	q2m = (realkind*)malloc(nspecies*sizeof(realkind));
}


void FieldDataCPU::allocate(PlasmaData* pdata_in)
{
	pdata = pdata_in;
	allocate(pdata->nx,pdata->ny,pdata->nz,pdata->nspecies);

	for(int i=0;i<pdata->nspecies;i++)
		q2m[i] = pdata->qspecies[i]/pdata->mspecies[i];

	ndimensions = pdata->ndimensions;


}

void FieldDataCPU::broadcast(FieldData** all_fields,ParallelInfo* myInfo)
{
	if(myInfo->tid == 0)
	{
		MPI_Bcast(all_fields[0]->Ex,alloc_size,MPI_REALKIND,0,MPI_COMM_WORLD);
		MPI_Bcast(all_fields[0]->Ey,alloc_size,MPI_REALKIND,0,MPI_COMM_WORLD);
		MPI_Bcast(all_fields[0]->Ez,alloc_size,MPI_REALKIND,0,MPI_COMM_WORLD);

		MPI_Bcast(all_fields[0]->Bx,alloc_size,MPI_REALKIND,0,MPI_COMM_WORLD);
		MPI_Bcast(all_fields[0]->By,alloc_size,MPI_REALKIND,0,MPI_COMM_WORLD);
		MPI_Bcast(all_fields[0]->Bz,alloc_size,MPI_REALKIND,0,MPI_COMM_WORLD);

	}

	copy_from(all_fields[0]);

}

template<int icomponent>
realkind& FieldDataCPU::getET(int ix, int iy, int iz)
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

template<int icomponent>
realkind& FieldDataCPU::getBT(int ix, int iy, int iz)
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

template<int icomponent>
realkind& FieldDataCPU::getETz(int ix, int iy, int iz)
{
	realkind* result;

	ix = (((ix&(nx-1)) + nx)&(nx-1));
	iy = (((iy&(ny-1)) + ny)&(ny-1));
	iz = (((iz&(nz-1)) + nz)&(nz-1));

	int iout = zorder(ix,iy,iz);

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

template<int icomponent>
realkind& FieldDataCPU::getBTz(int ix, int iy, int iz)
{
	realkind* result;

	ix = (((ix&(nx-1)) + nx)&(nx-1));
	iy = (((iy&(ny-1)) + ny)&(ny-1));
	iz = (((iz&(nz-1)) + nz)&(nz-1));

	int iout = zorder(ix,iy,iz);

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

realkind& FieldDataCPU::getE(int ix,int iy,int iz,int icomponent)
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

realkind& FieldDataCPU::getB(int ix,int iy,int iz,int icomponent)
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

template<int icomponent, enum FieldData_deriv ideriv>
realkind FieldDataCPU::intrpET(realkind x, realkind y, realkind z, int icellx, int icelly, int icellz)
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

static int icache = 0;

//double avx_hadd2(__m256d a)
//{
//
//		__m256d temp2 = _mm256_hadd_pd(a,a);
//
//
//		__m128d hi128 = _mm256_extractf128_pd(temp2,1);
//		__m128d dotproduct = _mm_add_pd(*(__m128d*)&temp2,hi128);
//
//		return *(double*)&dotproduct;
//
//}

template<int icomponent, enum FieldData_deriv ideriv>
realkind FieldDataCPU::intrpBT(realkind x, realkind y, realkind z, int icellx, int icelly, int icellz)
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

				realkind Etemp = getBT<0>(icellx+i,icelly,icellz);

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
			for(int i=0;i<2;i++)
			{
				realkind xp, yp, zp;

				realkind Etemp = getBT<1>(icellx+i,icelly,icellz);

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
		case 2:

			// z component
			for(int i=0;i<2;i++)
			{
				realkind xp, yp, zp;

				realkind Etemp = getBT<2>(icellx+i,icelly,icellz);

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

					realkind Etemp = getBT<0>(icellx+i,icelly+j,icellz);

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

					realkind Etemp = getBT<1>(icellx+i,icelly+j,icellz);

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

					realkind Etemp = getBT<2>(icellx+i,icelly+j,icellz);

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



realkind FieldDataCPU::intrpE(realkind x, realkind y, realkind z,
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


realkind FieldDataCPU::intrpB(realkind x, realkind y, realkind z,
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

void FieldDataCPU::copy_from(FieldData* src)
{
	memcpy(Bx,src->Bx,alloc_size*sizeof(realkind));
	memcpy(By,src->By,alloc_size*sizeof(realkind));
	memcpy(Bz,src->Bz,alloc_size*sizeof(realkind));

	memcpy(Ex,src->Ex,alloc_size*sizeof(realkind));
	memcpy(Ey,src->Ey,alloc_size*sizeof(realkind));
	memcpy(Ez,src->Ez,alloc_size*sizeof(realkind));

	memcpy(q2m,src->q2m,nspecies*sizeof(realkind));
}

realkind FieldDataCPU::evaluate_energy(void)
{
	double Ex_t, Ey_t, Ez_t, Emag;
	double Bx_t, By_t, Bz_t, Bmag;

	Emag = 0;
	Bmag = 0;
	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int i=0;i<nx;i++)
			{
				Ex_t = getE(i,j,k,0);
				Ey_t = getE(i,j,k,1);
				Ez_t = getE(i,j,k,2);

				Emag += (Ex_t*Ex_t + 0 + 0);

				Bx_t = getB(i,j,k,0);
				By_t = getB(i,j,k,1);
				Bz_t = getB(i,j,k,2);

				Bmag += (Bx_t*Bx_t + By_t*By_t + Bz_t*Bz_t);

			}
		}
	}

	realkind energy = 0.5 *(Emag/(qe2me*qe2me*epsilon_naught))*pdata->dxdi*pdata->dydi*pdata->dzdi;

	return energy;


}

void Average_Fields(FieldData* a,FieldData* b,FieldData* c)
{
	printf("averaging fields\n");

	if((a->alloc_size != b->alloc_size)||(a->alloc_size != c->alloc_size))
	{
		printf("Error: Can't average fields of different sizes\n");
		exit(1);
	}

	omp_set_num_threads(a->pdata->num_cores);

	realkind* As[6] = {a->Ex,a->Ey,a->Ez,a->Bx,a->By,a->Bz};
	realkind* Bs[6] = {b->Ex,b->Ey,b->Ez,b->Bx,b->By,b->Bz};
	realkind* Cs[6] = {c->Ex,c->Ey,c->Ez,c->Bx,c->By,c->Bz};

#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int ntotal = a->nx * a->ny;

		int stride = (ntotal + nthreads - 1)/nthreads;
		int istart = tid * stride;
		int iend = std::min((tid + 1) * stride-1,ntotal-1);

		for(int i=0;i<3;i++)
		for(int j=istart;j<=iend;j++)
		{
			int ix,iy;
			iy = j/a->nx;
			ix = j - iy*a->nx;
			c->getE(ix,iy,0,i) = 0.5*(a->getE(ix,iy,0,i) + b->getE(ix,iy,0,i));
		}

		for(int i=0;i<3;i++)
		for(int j=istart;j<=iend;j++)
		{
			int ix,iy;
			iy = j/a->nx;
			ix = j - iy*a->nx;
			c->getB(ix,iy,0,i) = 0.5*(a->getB(ix,iy,0,i) + b->getB(ix,iy,0,i));
		}


//
//		for(int i=0;i<6;i++)
//			for(int j=istart;j<=iend;j++)
//				Cs[i][j] = 0.5*(As[i][j] + Bs[i][j]);

//		for(int i=0;i<6;i++)
//			for(int j=0;j<a->alloc_size;j++)
//				Cs[i][j] = 0.5*(As[i][j] + Bs[i][j]);

//		for(int i=0;i<a->alloc_size;i++)
//			c->Ex[i] = 0.5f * (a->Ex[i] + b->Ex[i]);
//		for(int i=0;i<a->alloc_size;i++)
//			c->Ey[i] = 0.5f * (a->Ey[i] + b->Ey[i]);
//		for(int i=0;i<a->alloc_size;i++)
//			c->Ez[i] = 0.5f * (a->Ez[i] + b->Ez[i]);
//		for(int i=0;i<a->alloc_size;i++)
//			c->Bx[i] = 0.5f * (a->Bx[i] + b->Bx[i]);
//		for(int i=0;i<a->alloc_size;i++)
//			c->By[i] = 0.5f * (a->By[i] + b->By[i]);
//		for(int i=0;i<a->alloc_size;i++)
//			c->Bz[i] = 0.5f * (a->Bz[i] + b->Bz[i]);

	}

		printf("averaged fields\n");
	// Filter The Electric Field
	if(a->pdata->ndimensions == 1)
	{

		for(int k=0;k<c->nz;k++)
//#pragma omp parallel for
		for(int j=0;j<c->ny;j++)

		for(int i=0;i<c->nx;i++)
		{
			c->getE(i,j,k,0) = 0.125*(a->getE(i+1,j,k,0)
					+ 2.0*a->getE(i,j,k,0)
					+ a->getE(i-1,j,k,0)
					+ b->getE(i+1,j,k,0)
					+ 2.0*b->getE(i,j,k,0)
					+ b->getE(i-1,j,k,0));
		}
	}




}

void Difference_Fields(const FieldDataCPU& a,const FieldDataCPU& b,FieldDataCPU& c)
{

	if((a.alloc_size != b.alloc_size)||(a.alloc_size != c.alloc_size))
	{
		printf("Error: Can't average fields of different sizes\n");
		exit(1);
	}

	for(int i=0;i<a.alloc_size;i++)
	{
		c.Ex[i] = (a.Ex[i] - b.Ex[i]);
	}

	for(int i=0;i<a.alloc_size;i++)
	{
		c.Ey[i] = (a.Ey[i] - b.Ey[i]);
	}

	for(int i=0;i<a.alloc_size;i++)
	{
		c.Ez[i] = (a.Ez[i] - b.Ez[i]);
	}



	for(int i=0;i<a.alloc_size;i++)
	{
		c.Bx[i] = (a.Bx[i] - b.Bx[i]);
	}

	for(int i=0;i<a.alloc_size;i++)
	{
		c.By[i] = (a.By[i] - b.By[i]);
	}

	for(int i=0;i<a.alloc_size;i++)
	{
		c.Bz[i] = (a.Bz[i] - b.Bz[i]);
	}


}


void FieldDataCPU::init_plot()
{
	plot_handle = gnuplot_init();
}

void FieldDataCPU::reset_plot()
{
	gnuplot_resetplot(plot_handle);
}

void FieldDataCPU::close_plot()
{
	gnuplot_close(plot_handle);
}


void FieldDataCPU::plot(PlasmaData* pdata,int position, int plane = 0,
		const int Fplot = 0,const int icomponent = 0)
{
	/*
	 * plane = 0: xy plane
	 * plane = 1: xz plane
	 * plane = 2: yz plane
	 * Fplot = 0: Efield
	 * Fplot = 1: Bfield
	 *
	 */
	int nxp, nyp;
	int i,j,k;

	int* i_out_t;
	int* j_out_t;

	float dxp,dyp;
	float x0,y0;

	float* x_vals;
	float* y_vals;
	float* z_vals;

	float* vals_in;

	if(plane == 0)
	{
		nxp = pdata->nx;
		nyp = pdata->ny;

		dxp = pdata->dxdi;
		dyp = pdata->dydi;

		x0 = pdata->xmin;
		y0 = pdata->ymin;

		i_out_t = &i;
		j_out_t = &j;

		gnuplot_cmd(plot_handle,"set xlabel \"x\"");
		gnuplot_cmd(plot_handle,"set ylabel \"y\"");
	}
	else if(plane == 1)
	{
		nxp = pdata->nx;
		nyp = pdata->nz;

		dxp = pdata->dxdi;
		dyp = pdata->dzdi;

		x0 = pdata->xmin;
		y0 = pdata->zmin;

		i_out_t = &i;
		j_out_t = &k;

		gnuplot_cmd(plot_handle,"set xlabel \"x\"");
		gnuplot_cmd(plot_handle,"set ylabel \"z\"");
	}
	else if(plane == 2)
	{
		nxp = pdata->ny;
		nyp = pdata->nz;

		dxp = pdata->dydi;
		dyp = pdata->dzdi;

		x0 = pdata->ymin;
		y0 = pdata->zmin;

		i_out_t = &j;
		j_out_t = &k;

		gnuplot_cmd(plot_handle,"set xlabel \"y\"");
		gnuplot_cmd(plot_handle,"set ylabel \"z\"");
	}

	int& i_out = *i_out_t;
	int& j_out = *j_out_t;

	x_vals = (float*)malloc(nxp*nyp*sizeof(float));
	y_vals = (float*)malloc(nxp*nyp*sizeof(float));
	z_vals = (float*)malloc(nxp*nyp*sizeof(float));

	if(plane == 0)
	{
		k = position;

		for(j=0;j<pdata->ny;j++)
		{
			for(i=0;i<pdata->nx;i++)
			{
				x_vals[i_out] = dxp*i_out + x0;
				y_vals[j_out] = dyp*j_out + y0;
				if(Fplot == 0)
				{
					if(icomponent < 3)
						z_vals[i_out+nxp*j_out] = getE(i,j,k,0);
					else if(icomponent == 3)
						z_vals[i_out+nxp*j_out] = sqrt(
								powf(getE(i,j,k,0),2.0)
								+ powf(getE(i,j,k,1),2.0)
								+ powf(getE(i,j,k,2),2.0));
				}
				else
				{
					z_vals[i_out+nxp*j_out] = getB(i,j,k,icomponent);
				}

			}
		}
	}
	else if(plane == 1)
	{
		j = position;

		for(k=0;k<pdata->nz;k++)
		{
			for(i=0;i<pdata->nx;i++)
			{
				x_vals[i_out] = dxp*i_out + x0;
				y_vals[j_out] = dyp*j_out + y0;
				if(Fplot == 0)
				{
					if(icomponent < 3)
						z_vals[i_out+nxp*j_out] = getE(i,j,k,0);
					else if(icomponent == 3)
						z_vals[i_out+nxp*j_out] = sqrt(
								powf(getE(i,j,k,0),2.0)
								+ powf(getE(i,j,k,1),2.0)
								+ powf(getE(i,j,k,2),2.0));
				}
				else
				{
					z_vals[i_out+nxp*j_out] = getB(i,j,k,icomponent);

				}
			}
		}
	}
	else if(plane == 2)
	{
		i = position;

		for(k=0;k<pdata->nz;k++)
		{
			for(j=0;j<pdata->ny;j++)
			{
				x_vals[i_out] = dxp*i_out + x0;
				y_vals[j_out] = dyp*j_out + y0;
				if(Fplot == 0)
				{
					if(icomponent < 3)
						z_vals[i_out+nxp*j_out] = getE(i,j,k,0);
					else if(icomponent == 3)
						z_vals[i_out+nxp*j_out] = sqrt(
								powf(getE(i,j,k,0),2.0)
								+ powf(getE(i,j,k,1),2.0)
								+ powf(getE(i,j,k,2),2.0));
				}
				else
				{
					z_vals[i_out+nxp*j_out] = getB(i,j,k,icomponent);
				}
			}
		}
	}

	gnuplot_plot_xyz(plot_handle,x_vals,y_vals,z_vals,nxp,nyp,"Fields");


	free(x_vals);
	free(y_vals);
	free(z_vals);










}




















