#include "StressTallyGPU.h"
#include "ShapeFunctions.h"

//FUNCTION_TYPE
//StressTally::StressTally(realkind* _S2xx,
//						 realkind* _S2xy,
//						 realkind* _S2xz,
//						 realkind* _S2yy,
//						 realkind* _S2yz,
//						 realkind* _S2zz,
//						 int _nx,int _ny,int _nz,
//						 int _ndimensions) :
//						 S2xx(_S2xx),S2xy(_S2xy),S2xz(_S2xz),
//						 S2yy(_S2yy),S2yz(_S2yz),S2zz(_S2zz),
//						 nx(_nx),ny(_ny),nz(_nz),ndimensions(_ndimensions)
//{
//
//}
//
//FUNCTION_TYPE
//StressTally::StressTally(realkind* _S2xx,
//			 int _nx,int _ny,int _nz) :
//			 S2xx(_S2xx),
//			 nx(_nx),ny(_ny),nz(_nz),ndimensions(1)
//{
//
//}


__device__
void StressTallyGPU::tally1d1v(const realkind px,
		 const realkind vx,
		 const int ix_in,
		 const realkind scale)
{
	int ix;
	realkind vol_inv =  scale;



	for(int i=0;i<2;i++)
	{
		realkind xp;

		float temp;

		xp = i - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = vx*vx*S1_shape(xp)*vol_inv;


//		if(isnan(temp))
//			printf("Warning NAN current value at %i with %e, %e\n",ix,vx,vol_inv);

//#ifdef DOUBLE_PRECISION
//		atomicAddD(((double*)S2xx)+ix,temp);
//
//#else
		atomicAdd(S2xx+ix,(float)temp);
//#endif

	}



}

__device__
void StressTallyGPU::tally1d3v(const realkind px,
						 const realkind vx, const realkind vy, const realkind vz,
						 const int ix_in,
						 const realkind scale)
{

	int ix;
	realkind vol_inv =  scale;

	for(int i=-1;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i + 0.5 - px;


		ix = ix_in + i;


		ix = ((ix&(nx-1) + nx)&(nx-1));

		temp = S2_shape(xp)*vol_inv;
//#ifdef DOUBLE_PRECISION
//		atomicAddD(((double*)S2xx)+ix,temp*vx*vx);
//		atomicAddD(((double*)S2xy)+ix,temp*vx*vy);
//		atomicAddD(((double*)S2xz)+ix,temp*vx*vz);
//		atomicAddD(((double*)S2yy)+ix,temp*vy*vy);
//		atomicAddD(((double*)S2yz)+ix,temp*vy*vz);
//		atomicAddD(((double*)S2zz)+ix,temp*vz*vz);
//#else
		atomicAdd(S2xx+ix,(float)temp*vx*vx);
		atomicAdd(S2xy+ix,(float)temp*vx*vy);
		atomicAdd(S2xz+ix,(float)temp*vx*vz);
		atomicAdd(S2yy+ix,(float)temp*vy*vy);
		atomicAdd(S2yz+ix,(float)temp*vy*vz);
		atomicAdd(S2zz+ix,(float)temp*vz*vz);
//#endif


	}







}

__device__
void StressTallyGPU::tally2d3v(const realkind px, const realkind py,
						 const realkind vx, const realkind vy, const realkind vz,
						 const int ix_in, const int iy_in,
						 const realkind scale)
{

	int ix,iy;
	realkind vol_inv =  scale;


	for(int j=0;j<2;j++)
	{
		for(int i=0;i<2;i++)
		{
			realkind xp, yp;

			float temp;

			xp = i - px;
			yp = j - py;

			ix = ix_in + i;
			iy = iy_in + j;

			ix = ((ix&(nx-1) + nx)&(nx-1));
			iy = ((iy&(ny-1) + ny)&(ny-1));

			int iout = ix + nx*iy;

			temp = S1_shape(xp)*S1_shape(yp)*vol_inv;


			S2xx[iout] += temp*vx*vx;
			S2xy[iout] += temp*vx*vy;
			S2xz[iout] += temp*vx*vz;
			S2yy[iout] += temp*vy*vy;
			S2yz[iout] += temp*vy*vz;
			S2zz[iout] += temp*vz*vz;


		}
	}

}


__device__
void StressTallyGPU::tally(const realkind px, const realkind py, const realkind pz,
		 const realkind vx, const realkind vy, const realkind vz,
		 const int ix, const int iy, const int iz,
		 const realkind scale)
{
	if(ndimensions == 1)
	{
		if(nVel == 1)
		{
			tally1d1v(px,vx,ix,scale);
		}
		else if(nVel == 3)
		{
			tally1d3v(px,vx,vy,vz,ix,scale);
		}
	}
	else if(ndimensions == 2)
	{
		if(nVel == 3)
		{
			tally2d3v(px,py,vx,vy,vz,ix,iy,scale);
		}
	}
}
