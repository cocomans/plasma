#include "StressTally.h"
#include "ShapeFunctions.h"

FUNCTION_TYPE
StressTally::StressTally(realkind* stress_in,
						   int3 dims_in,
						   realkind spacingx,realkind spacingy,realkind spacingz,
						   int ndimensions_in)
{
	stress = stress_in;

	nx = dims_in.x;
	ny = dims_in.y;
	nz = dims_in.z;

	dx = spacingx;
	dy = spacingy;
	dz = spacingz;

	ndimensions = ndimensions_in;

}

FUNCTION_TYPE
StressTally::StressTally(){};

FUNCTION_TYPE
void StressTally::tally1d1v(const realkind px,
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

		temp = realkind(0.25)*vx*vx*S1_shape(xp)*vol_inv;

#ifdef GPU_CODE
	#ifdef DOUBLE_PRECISION
		atomicAddD(stress+ix,temp);
	#else
		atomicAdd(stress+ix,temp);
	#endif
#else
//		if(isnan(temp))
//			printf("Warning NAN current value at %i with %e, %e\n",ix,vx,vol_inv);

		stress[ix] += temp;
#endif

	}



}

FUNCTION_TYPE
void StressTally::tally(const realkind px, const realkind py, const realkind pz,
						 const realkind vx, const realkind vy, const realkind vz,
						 const int ix_in, const int iy_in, const int iz_in,
						 const realkind scale)
{

	int ix,iy,iz;
	realkind vol_inv =  scale;
	if(ndimensions == 1)
	{
		for(int k=0;k<nz;k++)
		{
			for(int j=0;j<ny;j++)
			{

				for(int i=0;i<2;i++)
				{
					realkind xp, yp, zp;

					realkind temp;

					xp = i - px;
					yp = j - py;
					zp = k - pz;

					ix = ix_in + i;
					iy = j;
					iz = k;

					ix = ((ix&(nx-1) + nx)&(nx-1));
					iy = ((iy&(ny-1) + ny)&(ny-1));
					iz = ((iz&(nz-1) + nz)&(nz-1));

					temp = 0.25f*vx*vx*S1_shape(xp)*vol_inv;

					stress[ix + nx*(iy + ny*(iz))] += temp;


				}
			}
		}


	}
	else
	{
		for(int k=-1;k<2;k++)
		{
			for(int j=-1;j<2;j++)
			{
				for(int i=0;i<2;i++)
				{
					realkind xp, yp, zp;

					realkind temp;

					xp = i - px;
					yp = j + 0.5 - py;
					zp = k + 0.5 - pz;

					ix = ix_in + i;
					iy = iy_in + j;
					iz = iz_in + k;

					ix = ((ix%nx + nx)%nx);
					iy = ((iy%ny + ny)%ny);
					iz = ((iz%nz + nz)%nz);

					temp = vx*vx*S1_shape(xp)*S1_shape(yp)*S1_shape(zp)*vol_inv;

					stress[ix + nx*(iy + ny*(iz))] += temp;


				}
			}
		}


	}


}
