/*-------------------------------------------------------------------------*/
/**
  @file		CurrentTally.cu
*/
/*--------------------------------------------------------------------------*/
#include "CurrentTally.h"
#include "ShapeFunctions.h"

#ifdef GPU_CODE
__device__
void atomicAddD(double* address,double value)
{
	unsigned long long oldval, newval, readback;
   oldval = __double_as_longlong(*address);
   newval = __double_as_longlong(__longlong_as_double(oldval) + value);
   while ((readback=atomicCAS((unsigned long long *)address, oldval, newval)) != oldval)
	 {
	  oldval = readback;
	  newval = __double_as_longlong(__longlong_as_double(oldval) + value);
	 }
}
#endif


FUNCTION_TYPE
CurrentTally::CurrentTally(realkind* currentx_in,
						   realkind* currenty_in,
						   realkind* currentz_in,
						   int3 dims_in,
						   realkind spacingx,realkind spacingy,realkind spacingz,
						   int ndimensions_in)
{
	currentx = currentx_in;
	currenty = currenty_in;
	currentz = currentz_in;

	nx = dims_in.x;
	ny = dims_in.y;
	nz = dims_in.z;

	dx = spacingx;
	dy = spacingy;
	dz = spacingz;

	ndimensions = ndimensions_in;

}

FUNCTION_TYPE
CurrentTally::CurrentTally(){}

FUNCTION_TYPE
void CurrentTally::tally1d1v(const realkind px,
		 const realkind vx,
		 const int ix_in,
		 const realkind scale)
{
	int ix;
	realkind vol_inv =  scale;



	for(int i=0;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = realkind(0.25)*vx*S1_shape(xp)*vol_inv;

#ifdef GPU_CODE
	#ifdef DOUBLE_PRECISION
		atomicAddD(currentx+ix,temp);
	#else
		atomicAdd(currentx+ix,temp);
	#endif
#else
//		if(isnan(temp))
//			printf("Warning NAN current value at %i with %e, %e\n",ix,vx,vol_inv);

		currentx[ix] += temp;
#endif

	}



}

FUNCTION_TYPE
void CurrentTally::tally1d2v(const realkind px,
		 const realkind vx,const realkind vy,
		 const int ix_in,
		 const realkind scale)
{

	int ix;
	realkind vol_inv =  scale;



	for(int i=0;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = realkind(0.25)*vx*S1_shape(xp)*vol_inv;

		currentx[ix] += temp;


	}

	for(int i=0;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = 0.25*vy*S1_shape(xp)*vol_inv;

		currenty[ix] += temp;


	}



}

__host__ __device__
void CurrentTally::tally1d3v(const realkind px,
		 const realkind vx,const realkind vy,const realkind vz,
		 const int ix_in,
		 const realkind scale)
{

	int ix;
	realkind vol_inv =  scale;



	for(int i=0;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = 0.25*vx*S1_shape(xp)*vol_inv;

		currentx[ix] += temp;


	}

	for(int i=0;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = 0.25*vy*S1_shape(xp)*vol_inv;

		currenty[ix] += temp;


	}

	for(int i=0;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = 0.25*vz*S1_shape(xp)*vol_inv;

		currentz[ix] += temp;


	}



}

__host__ __device__
void CurrentTally::tally2d2v(const realkind px,const realkind py,
		 const realkind vx,const realkind vy,
		 const int ix_in,const int iy_in,
		 const realkind scale)
{

	int ix,iy;
	realkind vol_inv =  scale;


	// x component

		for(int j=-1;j<2;j++)
		{
			for(int i=0;i<2;i++)
			{
				realkind xp, yp;

				realkind temp;

				xp = i - px;
				yp = j + 0.5 - py;

				ix = ix_in + i;
				iy = iy_in + j;

				ix = ((ix%nx + nx)%nx);
				iy = ((iy%ny + ny)%ny);

				temp = 0.5*vx*S1_shape(xp)*S2_shape(yp)*vol_inv;

				currentx[ix + nx*(iy)] += temp;


			}
		}


		// y component

		for(int j=0;j<2;j++)
		{
			for(int i=-1;i<2;i++)
			{
				realkind xp, yp;

				realkind temp;

				xp = i + 0.5 - px;
				yp = j - py;

				ix = ix_in + i;
				iy = iy_in + j;

				ix = ((ix%nx + nx)%nx);
				iy = ((iy%ny + ny)%ny);

				temp = 0.5f*vy*S2_shape(xp)*S1_shape(yp)*vol_inv;

				currenty[ix + nx*(iy)] += temp;


			}
		}



}

__host__ __device__
void CurrentTally::tally2d3v(const realkind px,const realkind py,
		 const realkind vx,const realkind vy,const realkind vz,
		 const int ix_in,const int iy_in,
		 const realkind scale)
{

	int ix,iy;
	realkind vol_inv =  scale;


	// x component

		for(int j=-1;j<2;j++)
		{
			for(int i=0;i<2;i++)
			{
				realkind xp, yp;

				realkind temp;

				xp = i - px;
				yp = j + 0.5 - py;

				ix = ix_in + i;
				iy = iy_in + j;

				ix = ((ix%nx + nx)%nx);
				iy = ((iy%ny + ny)%ny);

				temp = 0.5f*vx*S1_shape(xp)*S2_shape(yp)*vol_inv;
//				printf("temp = %f\n",temp);
				currentx[ix + nx*(iy)] += temp;


			}
		}


		// y component

		for(int j=0;j<2;j++)
		{
			for(int i=-1;i<2;i++)
			{
				realkind xp, yp;

				realkind temp;

				xp = i + 0.5 - px;
				yp = j - py;

				ix = ix_in + i;
				iy = iy_in + j;

				ix = ((ix%nx + nx)%nx);
				iy = ((iy%ny + ny)%ny);

				temp = 0.5f*vy*S2_shape(xp)*S1_shape(yp)*vol_inv;

				currenty[ix + nx*(iy)] += temp;


			}
		}


		// z component
		for(int j=-1;j<2;j++)
		{
			for(int i=-1;i<2;i++)
			{
				realkind xp, yp;

				realkind temp;

				xp = i + 0.5 - px;
				yp = j + 0.5 - py;

				ix = ix_in + i;
				iy = iy_in + j;


				ix = ((ix%nx + nx)%nx);
				iy = ((iy%ny + ny)%ny);


				temp = 0.5f*vz*S2_shape(xp)*S2_shape(yp)*vol_inv;

				currentz[ix + nx*(iy)] += temp;


			}
		}



}

__host__ __device__
void CurrentTally::tally3d3v(const realkind px,const realkind py,const realkind pz,
		 const realkind vx,const realkind vy,const realkind vz,
		 const int ix,const int iy,const int iz,
		 const realkind scale)
{

}

__host__ __device__
void CurrentTally::tally(const realkind px, const realkind py, const realkind pz,
						 const realkind vx, const realkind vy, const realkind vz,
						 const int ix_in, const int iy_in, const int iz_in,
						 const realkind scale)
{

	switch(ndimensions)
	{
	case 1:
		tally1d1v(px,vx,ix_in,scale);
		break;
	case 2:
		tally2d3v(px,py,vx,vy,vz,ix_in,iy_in,scale);
		break;
	case 3:
		tally3d3v(px,py,pz,vx,vy,vz,
				ix_in,iy_in,iz_in,scale);
		break;
	default:
		break;
	}

//	int ix,iy,iz;
//	realkind vol_inv =  scale;
//	if(ndimensions == 1)
//	{
//
//		for(int i=0;i<2;i++)
//		{
//			realkind xp;
//
//			realkind temp;
//
//			xp = i - px;
//
//			ix = ix_in + i;
//
//
//			ix = (((ix%(nx)) + nx)%(nx));
//
//			temp = 0.25*vx*S1_shape(xp)*vol_inv;
//
//			currentx[ix] += temp;
//
//
//		}
//
//
//
//	}
//	else
//	{
//
//		// x component
//		for(int k=-1;k<2;k++)
//		{
//			for(int j=-1;j<2;j++)
//			{
//				for(int i=0;i<2;i++)
//				{
//					realkind xp, yp, zp;
//
//					realkind temp;
//
//					xp = i - px;
//					yp = j + 0.5 - py;
//					zp = k + 0.5 - pz;
//
//					ix = ix_in + i;
//					iy = iy_in + j;
//					iz = iz_in + k;
//
//					ix = ((ix%nx + nx)%nx);
//					iy = ((iy%ny + ny)%ny);
//					iz = ((iz%nz + nz)%nz);
//
//					temp = vx*S1_shape(xp)*S2_shape(yp)*S2_shape(zp)*vol_inv;
//
//					currentx[ix + nx*(iy + ny*(iz))] += temp;
//
//
//				}
//			}
//		}
//
//		// y component
//		for(int k=-1;k<2;k++)
//		{
//			for(int j=0;j<2;j++)
//			{
//				for(int i=-1;i<2;i++)
//				{
//					realkind xp, yp, zp;
//
//					realkind temp;
//
//					xp = i + 0.5 - px;
//					yp = j - py;
//					zp = k + 0.5 - pz;
//
//					ix = ix_in + i;
//					iy = iy_in + j;
//					iz = iz_in + k;
//
//					ix = ((ix%nx + nx)%nx);
//					iy = ((iy%ny + ny)%ny);
//					iz = ((iz%nz + nz)%nz);
//
//					temp = vy*S2_shape(xp)*S1_shape(yp)*S2_shape(zp)*vol_inv;
//
//					currenty[ix + nx*(iy + ny*(iz))] += temp;
//
//
//				}
//			}
//		}
//
//		// z component
//		for(int k=0;k<2;k++)
//		{
//			for(int j=-1;j<2;j++)
//			{
//				for(int i=-1;i<2;i++)
//				{
//					realkind xp, yp, zp;
//
//					realkind temp;
//
//					xp = i + 0.5 - px;
//					yp = j + 0.5 - py;
//					zp = k - pz;
//
//					ix = ix_in + i;
//					iy = iy_in + j;
//					iz = iz_in + k;
//
//					ix = ((ix%nx + nx)%nx);
//					iy = ((iy%ny + ny)%ny);
//					iz = ((iz%nz + nz)%nz);
//
//					temp = vz*S2_shape(xp)*S2_shape(yp)*S1_shape(zp)*vol_inv;
//
//					currentz[ix + nx*(iy + ny*(iz))] += temp;
//
//
//				}
//			}
//		}
//	}


}
