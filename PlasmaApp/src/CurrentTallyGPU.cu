#include "CurrentTallyGPU.cuh"
#include "ShapeFunctions.h"


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

FUNCTION_TYPE
CurrentTallyGPU::CurrentTallyGPU(float* currentx_in,
		float* currenty_in,
		float* currentz_in,
			int _nx, int _ny, int _nz,
			int _ix0,int _iy0,int _iz0,
			int ndimensions_in)
	{
		currentxf=currentx_in;
		currentyf=currenty_in;
		currentzf=currentz_in;
		nx=_nx;ny=_ny;nz=_nz;
		ix0=_ix0;iy0=_iy0;iz0=_iz0;
		ndimensions=ndimensions_in;
	}

FUNCTION_TYPE
void CurrentTallyGPU::tally1d1v(const realkind px,
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

		if(isnan(temp))
			printf("Warning NAN current value at %i with %e, %e\n",ix,vx,vol_inv);

//#ifdef DOUBLE_PRECISION
//	atomicAddD(((double*)currentxf)+ix,temp);
//#else
	atomicAdd(currentxf+ix,temp);
//#endif

	}



}

FUNCTION_TYPE
void CurrentTallyGPU::tally1d2v(const realkind px,
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

//#ifdef DOUBLE_PRECISION
//	atomicAddD(currentx+ix,temp);
//#else
//	atomicAdd(currentx+ix,temp);
//#endif


	}

	for(int i=0;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = vy*S1_shape(xp)*vol_inv;

//#ifdef DOUBLE_PRECISION
//	atomicAddD(currenty+ix,temp);
//#else
//	atomicAdd(currenty+ix,temp);
//#endif


	}



}

FUNCTION_TYPE
void CurrentTallyGPU::tally1d3v(const realkind px,
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

//#ifdef DOUBLE_PRECISION
//	atomicAddD(((double*)currentxf)+ix,temp);
//#else
	atomicAdd(currentxf+ix,temp);
//#endif


	}

	for(int i=-1;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i +0.5 - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = vy*S2_shape(xp)*vol_inv;

//#ifdef DOUBLE_PRECISION
//	atomicAddD(((double*)currentyf)+ix,temp);
//#else
	atomicAdd(currentyf+ix,temp);
//#endif


	}

	for(int i=-1;i<2;i++)
	{
		realkind xp;

		realkind temp;

		xp = i + 0.5 - px;

		ix = ix_in + i;


		ix = (((ix&(nx-1)) + nx)&(nx-1));

		temp = vz*S2_shape(xp)*vol_inv;

//#ifdef DOUBLE_PRECISION
//	atomicAddD(((double*)currentzf)+ix,temp);
//#else
	atomicAdd(currentzf+ix,temp);
//#endif


	}



}

FUNCTION_TYPE
void CurrentTallyGPU::tally2d2v(const realkind px,const realkind py,
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

				ix = ix_in + i - ix0;
				iy = iy_in + j - iy0;

				temp = 0.5*vx*S1_shape(xp)*S2_shape(yp)*vol_inv;

//#ifdef DOUBLE_PRECISION
//	atomicAddD(currentx+ix+nx*iy,temp);
//#else
//	atomicAdd(currentx+ix+nx*iy,temp);
//#endif


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

				ix = ix_in + i - ix0;
				iy = iy_in + j - iy0;

				temp = 0.5f*vy*S2_shape(xp)*S1_shape(yp)*vol_inv;

//#ifdef DOUBLE_PRECISION
//	atomicAddD(currenty+ix+nx*iy,temp);
//#else
//	atomicAdd(currenty+ix+nx*iy,temp);
//#endif


			}
		}



}

FUNCTION_TYPE
void CurrentTallyGPU::tally2d3v(const realkind px,const realkind py,
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

				ix = ix_in + i - ix0;
				iy = iy_in + j - iy0;

				temp = 0.5f*vx*S1_shape(xp)*S2_shape(yp)*vol_inv;

				if((ix < 0)||(ix >= nx)||(iy < 0)||(iy >= ny))
				{
					printf("Error cell index %i %i is out of bounds!!!! \n",ix,iy);
				}
				else
				{

//#ifdef DOUBLE_PRECISION/
//	atomicAddD(currentx+ix+nx*iy,temp);
//#else
	atomicAdd(currentxf+ix+nx*iy,temp);
//#endif

				}
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

				ix = ix_in + i - ix0;
				iy = iy_in + j - iy0;

				temp = 0.5f*vy*S2_shape(xp)*S1_shape(yp)*vol_inv;

				if((ix < 0)||(ix >= nx)||(iy < 0)||(iy >= ny))
				{
					printf("Error cell index %i %i is out of bounds!!!! \n",ix,iy);
				}
				else
				{

//#ifdef DOUBLE_PRECISION
//	atomicAddD(currentx+ix+nx*iy,temp);
//#else
	atomicAdd(currentyf+ix+nx*iy,temp);
//#endif

				}
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

				ix = ix_in + i - ix0;
				iy = iy_in + j - iy0;



				temp = 0.5f*vz*S2_shape(xp)*S2_shape(yp)*vol_inv;

				if((ix < 0)||(ix >= nx)||(iy < 0)||(iy >= ny))
				{
					printf("Error cell index %i %i is out of bounds!!!! \n",ix,iy);
				}
				else
				{
//#ifdef DOUBLE_PRECISION
//	atomicAddD(currentz+ix+nx*iy,temp);
//#else
	atomicAdd(currentzf+ix+nx*iy,float(temp));
//#endif

				}
			}
		}



}

FUNCTION_TYPE
void CurrentTallyGPU::tally3d3v(const realkind px,const realkind py,const realkind pz,
		 const realkind vx,const realkind vy,const realkind vz,
		 const int ix,const int iy,const int iz,
		 const realkind scale)
{

}

FUNCTION_TYPE
void CurrentTallyGPU::tally(const realkind px, const realkind py, const realkind pz,
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


}
