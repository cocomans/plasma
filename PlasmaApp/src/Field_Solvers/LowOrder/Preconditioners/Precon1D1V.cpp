/*
 * Precon1D1V.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: payne
 */

#include "Precon1D1V.h"
#include "Mesh.h"

namespace HoLo
{
namespace LowOrder
{

Precon1D1V::~Precon1D1V()
{

}

Precon1D1V::Precon1D1V(SimParams* _params,
		DeviceNS::CPU::Mesh* _mesh,
			Epetra_Comm* _comm)
{
	comm = _comm;
	mesh = _mesh;

	nReduced = 2;

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

	red_map = new Epetra_Map(nElements,nReduced,0,*(comm));
	r = new Epetra_Vector(*red_map);
	z = new Epetra_Vector(*red_map);

}

void Precon1D1V::UpdateMatrix(const Epetra_Vector &xg)
{

}
int Precon1D1V::ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const
{
	return 0;
}


} /* namespace LowOrder */
} /* namespace HoLo */
