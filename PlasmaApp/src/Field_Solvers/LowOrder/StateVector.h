/*
 * StateVector.h
 *
 *  Created on: Dec 11, 2013
 *      Author: payne
 */

#ifndef STATEVECTOR_H_
#define STATEVECTOR_H_

#include "ScalarArray.h"
#include "SimParams.h"

namespace HoLo
{
typedef double real;
namespace LowOrder
{
class BoundaryValues;
}
class Mesh
{
public:
	int nx,ny,nz;
	int nBufx,nBufy,nBufz;
	int nSpatial;

	real Lx,Ly,Lz;
	template<class xGen>
	void GenMap(xGen xf)
	{

	}

	template<class xGen, class yGen>
	void GenMap(xGen xf,yGen yf)
	{

	}

	template<class xGen, class yGen,class zGen>
	void GenMap(xGen xf,yGen yf, zGen zf)
	{

	}

	real xMap(real u,real v=0.0, real w =0.0){return Lx/((real)nx)*u;}
	real yMap(real u,real v, real w =0.0){return Ly/((real)ny)*v;}
	real zMap(real u,real v, real w ){return Lz/((real)nz)*w;}

	real dxdu(int u,int v=0.0, int w =0.0){return Lx/((real)nx);}
	real dydv(int u,int v, int w =0.0){return Ly/((real)ny);}
	real dzdw(int u,int v, int w ){return Lz/((real)nz);}


};
class StateVector;

class SVGetter
{
public:

	virtual real operator()(const int &i, const int &j, const int&k);
};

template<class O>
class SVGetter_T : public SVGetter
{
public:
	const O op;

	SVGetter_T(O _op) : op(_op){};

	real operator()(const int &i, const int &j, const int&k)
	{
		return op(i,j,k);
	}
};

class StateVector
{
public:

	SimParams* params;
	Mesh* mesh;
	int nDep;
	int nAux;
	int nx,ny,nz;
	int nBufx,nBufy,nBufz;
	int nTx, nTy, nTz;

	int nAlloc_dep;
	int nAlloc_aux;
	int nElements;

	int pdeVar2genVar[SimParams::Vars::LF];




	virtual void Allocate(HoLo::Mesh* _mesh,SimParams* _params);

	virtual void Allocate(HoLo::Mesh* _mesh, StateVector* _old);

	virtual void PartialAllocate(SimParams* params);



	virtual void InitFiller();

	virtual real& GetDep(const int iType,const int& i,const int& j=0,const int& k=0);


	virtual real& GetAux(const int iType,const int& i,const int& j=0,const int& k=0);

	// Returns from general indexing located in simparams
	virtual real& Get(const int iType,const int& i,const int& j=0,const int& k=0);

	virtual void ApplyBC(HoLo::LowOrder::BoundaryValues* BCs);

	virtual void CopyDepFrom(const Epetra_Vector* _src);

	virtual void CopyAuxFrom(const Epetra_Vector* _src);

	virtual void CopyDepTo(Epetra_Vector* _src);

	virtual void CopyAuxTo(Epetra_Vector* _src);

	virtual void CopyFrom(StateVector* _src);

	virtual void Set(real val);




	Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::C,SVGetter> *ne,*ni,*phi;

	Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::C,SVGetter> *pe_u,*pi_u,*Eu,*Au;

	Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::C,SVGetter> *pe_v,*pi_v,*Ev,*Av;

	Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::F,SVGetter> *pe_w,*pi_w,*Ew,*Aw;

	Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::F,SVGetter> *Bu;

	Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::C,
	Util::PDE::StagMesh::F,SVGetter> *Bv;

	Util::PDE::ScalarArray<real,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::F,
	Util::PDE::StagMesh::C,SVGetter> *Bw;


	Epetra_Vector* Dep_vector;
	Epetra_Vector* Aux_vector;

	Epetra_BlockMap* Dep_map;
	Epetra_BlockMap* Aux_map;


};
}

#endif /* STATEVECTOR_H_ */
