/*
 * StressTensor.cpp
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#include "StressTensor.h"
#include "Mesh.h"

namespace HoLo
{
typedef DeviceNS::CPU::real real;

void StressTensor::Allocate(SimParams* _params)
{
	params = _params;
	mesh = params->mesh;
	nVec = LF;

	nx = mesh->nX;
	ny = mesh->nY;
	nz = mesh->nZ;

	nBufx = mesh->nBufx;
	nBufy = mesh->nBufy;
	nBufz = mesh->nBufz;

	nTx = nx + 2*nBufx;
	nTy = ny + 2*nBufy;
	nTz = nz + 2*nBufz;

	nElements = nTx*nTy*nTz;

	map = new Epetra_BlockMap(nElements,2*nVec,0,*(params->comm));

	vec = new Epetra_Vector(*map);



}

void StressTensor::InitFiller()
{
	real k = 2.0*M_PI;

	for(int k=0;k<nz;k++)
	for(int j=0;j<ny;j++)
	for(int i=0;i<nx;i++)
	{
		real x = i/((real)nx-1.0);
		for(int l=0;l<nVec;l++)
			for(int m=0;m<2;m++)
			Get(l,m,i,j,k) = sin(k*x);

	}
}

double& StressTensor::Get(const int& iType,const int& iSpecies,args_ijk)
{
	const int hash = i+nBufx+nTx*(j+nBufy+nTy*(k+nBufz));
	return (*vec)[map->FirstPointInElement(hash)+iType + nVec*iSpecies];
}


} /* namespace HoLo */
