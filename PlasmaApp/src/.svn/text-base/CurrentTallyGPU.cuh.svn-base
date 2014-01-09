#ifndef Current_Tally_GPU_H
#define Current_Tally_GPU_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "PlasmaData.h"


#ifdef GPU_CODE
#define FUNCTION_TYPE __attribute__((device))
#else
#define FUNCTION_TYPE __attribute__((host,device))
#endif


class CurrentTallyGPU
{
public:

	FUNCTION_TYPE
	CurrentTallyGPU(float* currentx_in,
			float* currenty_in,
			float* currentz_in,
			int _nx, int _ny, int _nz,
			int _ix0,int _iy0,int _iz0,
			int ndimensions_in);

	FUNCTION_TYPE
	CurrentTallyGPU(){};

	FUNCTION_TYPE
	~CurrentTallyGPU(){}

	FUNCTION_TYPE
	void tally1d1v(const realkind px,
			 const realkind vx,
			 const int ix_in,
			 const realkind scale);
	FUNCTION_TYPE
	void tally1d2v(const realkind px,
			 const realkind vx,const realkind vy,
			 const int ix_in,
			 const realkind scale);
	FUNCTION_TYPE
	void tally1d3v(const realkind px,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,
			 const realkind scale);
	FUNCTION_TYPE
	void tally2d2v(const realkind px,const realkind py,
			 const realkind vx,const realkind vy,
			 const int ix_in,const int iy_in,
			 const realkind scale);
	FUNCTION_TYPE
	void tally2d3v(const realkind px,const realkind py,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,const int iy_in,
			 const realkind scale);
	FUNCTION_TYPE
	void tally3d3v(const realkind px,const realkind py,const realkind pz,
			 const realkind vx,const realkind vy,const realkind vz,
			 const int ix_in,const int iy_in,const int iz_in,
			 const realkind scale);


	FUNCTION_TYPE
	void tally(const realkind px, const realkind py, const realkind pz,
			 const realkind vx, const realkind vy, const realkind vz,
			 const int ix, const int iy, const int iz,
			 const realkind scale);


	float* currentxf;
	float* currentyf;
	float* currentzf;

//	int nxt,nyt,nzt;
	int ix0,iy0,iz0;
//	int ndimensionst;
	int nx,ny,nz;
	int ndimensions;
	int nVel;


};









#endif /* Current_Tally_GPU_H */
