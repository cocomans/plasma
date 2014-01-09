#ifndef LO_PROBLEM_H_
#define LO_PROBLEM_H_

// Trilinos Includes
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "NOX_Epetra.H"

// User Includes
#include "SimParams.h"
#include "DiscretePDE.h"
#include "PhysBasedPrec.h"

class LowOrderProblem:
public NOX::Epetra::Interface::Required,
public NOX::Epetra::Interface::Preconditioner
{
 public:

  LowOrderProblem ( const Teuchos::RCP<SimParams>& params, Epetra_Comm* comm );
  ~LowOrderProblem (  );

  void Init();//only used to load linear IC for timestepping tests

  Teuchos::RCP<Epetra_Vector> GetX();

  void SetXold(const Epetra_Vector &xold);
  Teuchos::RCP<Epetra_Vector> GetXold ();

  void SetConsTerms(const Epetra_Vector &cons);
  Teuchos::RCP<Epetra_Vector> GetConsTerms ();

  Teuchos::RCP<PhysBasedPrec> GetPreconditioner();
  Teuchos::RCP<MapManager> GetMapManager();

  // inherited functions
  bool computeF( const Epetra_Vector &x,
		 Epetra_Vector &F,
		 const FillType fillFlag );

  bool computePreconditioner( const Epetra_Vector &x,
			      Epetra_Operator &M,
			      Teuchos::ParameterList *precParams );
  
 private:
  Teuchos::RCP<MapManager> mapManager;  
  Epetra_Comm*             comm;
  Teuchos::RCP<SimParams>  simParams;

  Teuchos::RCP<DiscretePDE>   pde;
  Teuchos::RCP<PhysBasedPrec> precond;    

  Teuchos::RCP<Epetra_Vector> X;
  Teuchos::RCP<Epetra_Vector> Xold;
  Teuchos::RCP<Epetra_Vector> ConsTerms;

};

#endif // LO_PROBLEM_H_
