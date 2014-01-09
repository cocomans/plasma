#include "ChargeTallyGPU.h"
#include "ShapeFunctions.h"


ChargeTallyGPU::ChargeTallyGPU(float* charge_in,
						   int3 dims_in,
						   realkind spacingx,realkind spacingy,realkind spacingz,
						   int ndimensions_in)
{
	charge = charge_in;

	nx = dims_in.x;
	ny = dims_in.y;
	nz = dims_in.z;

	ndimensions = ndimensions_in;

}


ChargeTallyGPU::ChargeTallyGPU(){};



void ChargeTallyGPU::tally1d(const realkind px,
		 const int ix_in,
		 const realkind scale)
{
	int ix;
	realkind vol_inv =  scale;



	for(int i=-2;i<3;i++)
	{
		realkind xp;

		realkind temp;

		xp = i - px + realkind(0.5);

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = S2_shape(xp)*vol_inv;



//#ifdef DOUBLE_PRECISION
//	atomicAddD(((double*)charge)+ix,temp);
//#else
	atomicAdd(charge+ix,temp);
//#endif
	}



}




void ChargeTallyGPU::tally2d(const realkind px,const realkind py,
		 const int ix_in,const int iy_in,
		 const realkind scale)
{

	int ix,iy;
	realkind vol_inv =  scale;


	// x component

		for(int j=-2;j<3;j++)
		{
			for(int i=-2;i<3;i++)
			{
				realkind xp, yp;

				realkind temp;

				xp = i - px +0.5f;
				yp = j - py +0.5f;

				ix = ix_in + i;
				iy = iy_in + j;

				ix = ((ix&(nx-1) + nx)&(nx-1));
				iy = ((iy&(ny-1) + ny)&(ny-1));

				temp = 0.5f*S2_shape(xp)*S2_shape(yp)*vol_inv;

//		if(isnan(temp))
//			printf("Warning NAN charge value at %i with %e\n",ix,vol_inv);

		charge[ix + nx*(iy)] += temp;



			}
		}



}


void ChargeTallyGPU::tally3d(const realkind px,const realkind py,const realkind pz,
		 const int ix_in,const int iy_in,const int iz_in,
		 const realkind scale)
{

	int ix,iy,iz;
	realkind vol_inv =  scale;


	// x component
	for(int k=-2;k<3;k++)
	{
		for(int j=-2;j<3;j++)
		{
			for(int i=-2;i<3;i++)
			{
				realkind xp, yp, zp;

				realkind temp;
				xp = i - px +0.5;
				yp = j - py +0.5;
				zp = k - pz +0.5;

				ix = ix_in + i;
				iy = iy_in + j;
				iz = iz_in + k;

				ix = (((ix&(nx-1)) + nx)&(nx-1));
				iy = (((iy&(ny-1)) + ny)&(ny-1));
				iz = (((iz&(nz-1)) + nz)&(nz-1));

				temp = S2_shape(xp)*S2_shape(yp)*S2_shape(zp)*vol_inv;
//		if(isnan(temp))
//			printf("Warning NAN charge value at %i with %e\n",ix,vol_inv);

		charge[ix + nx*(iy+ny*iz)] += temp;



			}
		}
	}

}

void ChargeTallyGPU::tally(const realkind px, const realkind py, const realkind pz,
						 const int ix_in, const int iy_in, const int iz_in,
						 const realkind scale)
{

	switch(ndimensions)
	{
	case 1:
		tally1d(px,ix_in,1.0);
		break;
	case 2:
		tally2d(px,py,ix_in,iy_in,1.0);
		break;
	case 3:
		tally3d(px,py,pz,
				ix_in,iy_in,iz_in,1.0);
		break;
	default:
		break;
	}
}
