/*
 * Preconditioner.h
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#ifndef PRECONDITIONER_H_
#define PRECONDITIONER_H_
#include "SimParams.h"

class MLSolver;
namespace DeviceNS
{
namespace CPU
{
	class Mesh;
}
}
namespace HoLo
{
namespace LowOrder
{

class Preconditioner:
		public Epetra_Operator
{
public:

	Preconditioner(){};
	virtual ~Preconditioner(){}

//	Preconditioner(const Teuchos::RCP<SimParams> &params,
//						Epetra_Comm* comm);

	virtual void UpdateMatrix(const Epetra_Vector &xg){};
	virtual int ApplyInverse( const Epetra_MultiVector &r, Epetra_MultiVector &z) const{return 0;};;

	// Inherited from Epetra_Operator
	virtual double NormInf() const;
	virtual int SetUseTranspose(bool UseTranspose);
	virtual int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
	virtual const char* Label() const;
	virtual bool UseTranspose() const;
	virtual bool HasNormInf() const;
	virtual const Epetra_Comm& Comm () const ;
	virtual const Epetra_Map& OperatorDomainMap () const;
	virtual const Epetra_Map& OperatorRangeMap () const;


	virtual void ImportParams(){};


	Epetra_Comm* comm;
	SimParams* simParams;
	DeviceNS::CPU::Mesh* mesh;

	MLSolver* mlSolver;


	Epetra_CrsMatrix* M;

	double me_h;
	double mi_h;
	double qe_h;
	double qi_h;
	double xi;
	double dx;
	double dt;
	bool pre;

	int nReduced;
	int nx,ny,nz;
	int nBufx,nBufy,nBufz;
	int nTx, nTy, nTz;

	int nElements;

	Epetra_Vector* r;
	Epetra_Vector* z;
	Epetra_Map* red_map;
};

} /* namespace LowOrder */
} /* namespace HoLo */
#endif /* PRECONDITIONER_H_ */
