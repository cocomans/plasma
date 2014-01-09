/*-------------------------------------------------------------------------*/
/**
  @file		CurrentTallyCPU.cpp
*/
/*--------------------------------------------------------------------------*/
#include "CurrentTallyCPU.h"
#include "ShapeFunctions.h"

FUNCTION_TYPE
CurrentTallyCPU::CurrentTallyCPU(realkind* currentx_in,
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
CurrentTallyCPU::CurrentTallyCPU(){}

FUNCTION_TYPE
void CurrentTallyCPU::tally1d1v(const realkind px,
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

		temp = vx*S1_shape(xp)*vol_inv;

//		if(isnan(temp))
//			printf("Warning NAN current value at %i with %e, %e\n",ix,vx,vol_inv);
		currentx[ix] += temp;

	}



}

FUNCTION_TYPE
void CurrentTallyCPU::tally1d2v(const realkind px,
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

		temp = vx*S1_shape(xp)*vol_inv;

		currentx[ix] += temp;

	}

	for(int i=0;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = vy*S1_shape(xp)*vol_inv;

		currenty[ix] += temp;


	}



}

FUNCTION_TYPE
void CurrentTallyCPU::tally1d3v(const realkind px,
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

		temp = vx*S1_shape(xp)*vol_inv;

		currentx[ix] += temp;


	}

	for(int i=-1;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i +0.5 - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = vy*S2_shape(xp)*vol_inv;

		currenty[ix] += temp;


	}

	for(int i=-1;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i +0.5 - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = vz*S2_shape(xp)*vol_inv;

		currentz[ix] += temp;


	}



}

FUNCTION_TYPE
void CurrentTallyCPU::tally2d2v(const realkind px,const realkind py,
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

				temp = vx*S1_shape(xp)*S2_shape(yp)*vol_inv;


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

				temp = vy*S2_shape(xp)*S1_shape(yp)*vol_inv;

				currenty[ix + nx*(iy)] += temp;


			}
		}



}

FUNCTION_TYPE
void CurrentTallyCPU::tally2d3v(const realkind px,const realkind py,
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
				yp = j + HALF_r - py;

				ix = ix_in + i;
				iy = iy_in + j;

				ix = (((ix&(nx-1)) + nx)&(nx-1));
				iy = (((iy&(ny-1)) + ny)&(ny-1));
				temp = vx*S1_shape(xp)*S2_shape(yp)*vol_inv;
//				if(isnan(temp))
//					printf("Warning NAN current value at %i with %e, %e\n",ix,vx,vol_inv);

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

				xp = i + HALF_r - px;
				yp = j - py;

				ix = ix_in + i;
				iy = iy_in + j;

				ix = (((ix&(nx-1)) + nx)&(nx-1));
				iy = (((iy&(ny-1)) + ny)&(ny-1));
				temp = vy*S2_shape(xp)*S1_shape(yp)*vol_inv;
//				if(isnan(temp))
//					printf("Warning NAN current value at %i with %e, %e\n",ix,vx,vol_inv);

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

				xp = i + HALF_r - px;
				yp = j + HALF_r - py;

				ix = ix_in + i;
				iy = iy_in + j;


				ix = (((ix&(nx-1)) + nx)&(nx-1));
				iy = (((iy&(ny-1)) + ny)&(ny-1));


				temp = vz*S2_shape(xp)*S2_shape(yp)*vol_inv;
//				if(isnan(temp))
//					printf("Warning NAN current value at %i with %e, %e\n",ix,vx,vol_inv);

				currentz[ix + nx*(iy)] += temp;


			}
		}



}

FUNCTION_TYPE
void CurrentTallyCPU::tally3d3v(const realkind px,const realkind py,const realkind pz,
		 const realkind vx,const realkind vy,const realkind vz,
		 const int ix,const int iy,const int iz,
		 const realkind scale)
{

}

FUNCTION_TYPE
void CurrentTallyCPU::tally(const realkind px, const realkind py, const realkind pz,
						 const realkind vx, const realkind vy, const realkind vz,
						 const int ix_in, const int iy_in, const int iz_in,
						 const realkind scale)
{

	switch(ndimensions)
	{
	case 1:
		tally1d3v(px,vx,vy,vz,ix_in,scale);
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


}
