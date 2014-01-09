#ifndef LO_PROBLEM_2D_H_
#define LO_PROBLEM_2D_H_

// Trilinos Includes
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "NOX_Epetra.H"

// User Includes
#include "../HoLoInterface/SimParams.h"
#include "DiscretePDE2D.h"
#include "PhysBasedPrec2D.h"
#include "EpVecWrapper.h"

class LowOrderProblem2D:
public NOX::Epetra::Interface::Required,
public NOX::Epetra::Interface::Preconditioner
{
 public:

	LowOrderProblem2D ( const Teuchos::RCP<SimParams>& params, Epetra_Comm* comm );
  ~LowOrderProblem2D (  );

  void Init();//only used to load linear IC for timestepping tests

  Teuchos::RCP<Epetra_Vector> GetX();
  Teuchos::RCP<EpVecWrapper> GetXW();


  void SetXold(const Epetra_Vector &xold);
  Teuchos::RCP<Epetra_Vector> GetXold ();
  Teuchos::RCP<EpVecWrapper> GetXoldW ();


  void SetConsTerms(const Epetra_Vector &cons);
  Teuchos::RCP<Epetra_Vector> GetConsTerms ();
  Teuchos::RCP<EpVecWrapper> GetConsTermsW ();

  Teuchos::RCP<PhysBasedPrec2D> GetPreconditioner();
  Teuchos::RCP<MapManager2D> GetMapManager();

  // inherited functions
  bool computeF( const Epetra_Vector &x,
		 Epetra_Vector &F,
		 const FillType fillFlag );

  bool computePreconditioner( const Epetra_Vector &x,
			      Epetra_Operator &M,
			      Teuchos::ParameterList *precParams );
  


 private:
  Teuchos::RCP<MapManager2D> mapManager;
  Epetra_Comm*             comm;
  Teuchos::RCP<SimParams>  simParams;

  Teuchos::RCP<DiscretePDE2D>   pde;
  Teuchos::RCP<PhysBasedPrec2D> precond;

  Teuchos::RCP<EpVecWrapper> X;
  Teuchos::RCP<EpVecWrapper> Xold;
  Teuchos::RCP<EpVecWrapper> Res;
  Teuchos::RCP<EpVecWrapper> ConsTerms;




};

#endif // LO_PROBLEM_2D_H_
