#ifndef MESH_INFO_H
#define MESH_INFO_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "PlasmaData.h"

class MeshInfo
{
	virtual ~MeshInfo(){};

	virtual void init(PlasmaData* pdata){};

	// Methods to get cell dimensions
	virtual realkind get_divs(int ix){return 0.0;};
	virtual realkind get_divs(int ix, int iy){return 0.0;};
	virtual realkind get_divs(int ix, int iy, int iz){return 0.0;};

	virtual typevecN<realkind,1> get_realpos(int ix, float px){};
	virtual typevecN<realkind,2> get_realpos(int ix, int iy, float px, float py){};
	virtual typevecN<realkind,3> get_realpos(int ix, int iy, int iz, float px, float py, float pz){};

	virtual void get_logpos(float x, float y, float z,
							float& px, float& py, float& pz,
							int& ix, int& iy, int& iz){};


};









#endif /* MESH_INFO_H */
