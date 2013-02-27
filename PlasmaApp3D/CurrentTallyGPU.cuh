#ifndef Current_Tally_GPU_H
#define Current_Tally_GPU_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "PlasmaData.h"
#include "CurrentTally.h"



class CurrentTallyGPU : public CurrentTally
{
public:
	__host__ __device__
	CurrentTallyGPU(realkind* currentx_in,
			realkind* currenty_in,
			realkind* currentz_in,
			int3 dims_in,
			realkind spacingx,realkind spacingy,realkind spacingz,
			int ndimensions_in);
	__host__ __device__
	CurrentTallyGPU();
	__device__
	void tally1d1v(const realkind px,
			 const realkind vx,
			 const int ix_in,
			 const realkind scale);
	__device__
	void tally1d2v(const realkind px,
			 const realkind vx,const realkind vy,
			 const int ix_in,
			 const realkind scale);
	__device__
	void tally1d3v(const realkind px,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,
			 const realkind scale);
	__device__
	void tally2d2v(const realkind px,const realkind py,
			 const realkind vx,const realkind vy,
			 const int ix_in,const int iy_in,
			 const realkind scale);
	__device__
	void tally2d3v(const realkind px,const realkind py,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,const int iy_in,
			 const realkind scale);
	__device__
	void tally3d3v(const realkind px,const realkind py,const realkind pz,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,const int iy_in,const int iz_in,
			 const realkind scale);


	__device__
	void tally(const realkind px, const realkind py, const realkind pz,
			 const realkind vx, const realkind vy, const realkind vz,
			 const int ix, const int iy, const int iz,
			 const realkind scale);


	realkind* currentx;
	realkind* currenty;
	realkind* currentz;

	int nx,ny,nz;
	int nptcls;
	int ndimensions;
	realkind dx,dy,dz;


};









#endif /* Current_Tally_GPU_H */
