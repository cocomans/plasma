/*
 * StateVectorT.cpp
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#include "StateVectorT.h"
#include "BoundaryValues.h"
#include <math.h>

namespace HoLo
{


template<class PDE>
void StateVector_T<PDE>::Allocate(Mesh* _mesh)
{
	mesh = _mesh;
	nDep = PDE::Dep::LF;
	nAux = PDE::Aux::LF;

	nx = mesh->nx;
	ny = mesh->ny;
	nz = mesh->nz;

	nBufx = mesh->nBufx;
	nBufy = mesh->nBufy;
	nBufz = mesh->nBufz;

	nTx = nx + 2*nBufx;
	nTy = ny + 2*nBufy;
	nTz = nz + 2*nBufz;

	nElements = nTx*nTy*nTz;

	Dep_map = new Epetra_BlockMap(nElements,nDep,0,*(params->comm));
	Aux_map = new Epetra_BlockMap(nElements,nAux,0,*(params->comm));

	Dep_vector = new Epetra_Vector(*Dep_map);
	Aux_vector = new Epetra_Vector(*Aux_map);


	ne = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::ne,i,j,k);}));

	ni = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::ni,i,j,k);}));

	phi = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::phi,i,j,k);}));

	pe_u = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::F,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::pe_u,i,j,k);}));

	pi_u = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::F,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::pi_u,i,j,k);}));

	Eu = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::F,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetAux(PDE::Aux::Eu,i,j,k);}));

	Au = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::F,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetAux(PDE::Aux::Au,i,j,k);}));

	pe_v = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::F,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::pe_v,i,j,k);}));

	pi_v = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::F,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::pi_v,i,j,k);}));

	Ev = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::F,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetAux(PDE::Aux::Ev,i,j,k);}));

	Av = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::F,
			Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::Av,i,j,k);}));

	pe_w = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::F,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::pe_w,i,j,k);}));

	pi_w = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::F,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::pi_w,i,j,k);}));

	Ew = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::F,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetAux(PDE::Aux::Ew,i,j,k);}));

	Aw = new Util::PDE::ScalarArray<real,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::C,
			Util::PDE::StagMesh::F,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetDep(PDE::Dep::Aw,i,j,k);}));

	Bu = new Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::F,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetAux(PDE::Aux::Bu,i,j,k);}));

	Bv = new Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::F,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetAux(PDE::Aux::Bv,i,j,k);}));

	Bw = new Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::C,SVGetter>(
			SVGetter([&](const int& i, const int j, const int k){return GetAux(PDE::Aux::Bw,i,j,k);}));


}

template<class PDE>
void StateVector_T<PDE>::InitFiller()
{
	real k = 2.0*M_PI;

	for(int k=0;k<nz;k++)
	for(int j=0;j<ny;j++)
	for(int i=0;i<nx;i++)
	{
		real x = i/((real)nx-1.0);
		for(int l=0;l<nDep;l++)
			GetDep(l,i,j,k) = sin(k*x);

		for(int l=0;l<nAux;l++)
			GetAux(l,i,j,k) = sin(k*x);
	}
}

template<class PDE>
double& StateVector_T<PDE>::GetDep(const int iType,const int& i,const int& j,const int& k)
{
	const int hash = i+nBufx+nTx*(j+nBufy+nTy*(k+nBufz));
	return (*Dep_vector)[Dep_map->FirstPointInElement(hash)+iType];
}

template<class PDE>
double& StateVector_T<PDE>::GetAux(const int iType,const int& i,const int& j,const int& k)
{
	const int hash = i+nBufx+nTx*(j+nBufy+nTy*(k+nBufz));
	return (*Aux_vector)[Aux_map->FirstPointInElement(hash)+iType];
}

template<class PDE>
double& StateVector_T<PDE>::Get(const int iType,const int& i,const int& j,const int& k)
{
	int iType2 = PDE::pdeVar2genVar[iType];
	int iStor = iType2 >> (8*sizeof(int)-1);
	iType2 = abs(iType2);

	double* result;

	if(iStor)
		result = &GetAux(iType2,i,j,k);
	else
		result = &GetDep(iType2,i,j,k);

	return *result;
}

template<class PDE>
void StateVector_T<PDE>::ApplyBC(HoLo::LowOrder::BoundaryValues* BCs)
{
	// east west BC's
	for(int l=0;l<SimParams::Vars::LF;l++)
	{
		BCs->ApplyBCs([&](int i, int j, int k){return Get(l,i,j,k);},l);
	}
}

} /* namespace HoLo */
