/*-------------------------------------------------------------------------*/
/**
	@file		ClusterInfo.cuh
	@author	J. Payne
	@date		03/01/2013
	@brief	Class to contain cell cluster information for GPU domain decomposition

*/
/*--------------------------------------------------------------------------*/
#ifndef CLUSTER_INFO_H
#define CLUSTER_INFO_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "PlasmaData.h"

class ClusterInfo
{
public:
	/// First Particle in the cluster list, will contain -1 if no particles are in this cluster
	int ifirst;
	/// Last particle in the cluster list
	int ilast;

	int clusterid;

};

#endif
