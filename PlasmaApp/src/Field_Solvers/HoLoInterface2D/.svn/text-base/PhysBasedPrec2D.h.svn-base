#ifndef PDE_PRECONDITIONER_2D_H_
#define PDE_PRECONDITIONER_2D_H_

// Trilinos Includes
//#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Operator.h"
#include "NOX_Epetra.H"
#include "AztecOO.h"

// User Includes
#include "../HoLoInterface/SimParams.h"
//#include "DiffMath.h"
#include "MapManager2D.h"
#include "../HoLoInterface/MLSolver.h"

class PhysBasedPrec2D:
public Epetra_Operator
{
 public:
	PhysBasedPrec2D(const Teuchos::RCP<MapManager2D>& map_manager,
		const Teuchos::RCP<SimParams> &params,
		Epetra_Comm* comm);

  ~PhysBasedPrec2D();

  void UpdateMatrix(const Epetra_Vector &xg);
  int ApplyInverse( const Epetra_MultiVector &r, Epetra_MultiVector &z) const;;

  // Inherited from Epetra_Operator
  double NormInf() const;
  int SetUseTranspose(bool UseTranspose);
  int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
  const char* Label() const;
  bool UseTranspose() const;
  bool HasNormInf() const;
  const Epetra_Comm& Comm () const ;
  const Epetra_Map& OperatorDomainMap () const;
  const Epetra_Map& OperatorRangeMap () const;
   
 private:

  void ImportParams();
  
  Epetra_Comm* comm;
  Teuchos::RCP<SimParams> simParams;
  Teuchos::RCP<MapManager2D> mapManager;

//  Teuchos::RCP<MLSolver> mlSolver;
//
//  Teuchos::RCP<Epetra_Vector> ResSd;
//  Teuchos::RCP<DiffMath> ResOv;
//  Teuchos::RCP<Epetra_Vector> ResRed;
//
//  Teuchos::RCP<Epetra_Vector> SolSd;
//  Teuchos::RCP<DiffMath> SolOv;
//  Teuchos::RCP<Epetra_Vector> SolRed;

//  Teuchos::RCP<DiffMath> curr_X;

  Teuchos::RCP<Epetra_CrsMatrix> M;

  double me_h;
  double mi_h;
  double qe_h;
  double qi_h;
  double xi;
  double dx;
  double dt;
  bool pre;
};

#endif // PDE_PRECONDITIONER_2D_H_
