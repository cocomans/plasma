#ifndef JFNK_H_
#define JFNK_H_


// User Includes
#include "SimParams.h"
#include "LowOrderProblem.h"

class JFNK
{
 public:
  JFNK( const Teuchos::RCP<LowOrderProblem>& problem,
	const Teuchos::RCP<SimParams>& Params,
	Epetra_Comm* Comm );
  ~JFNK();

  void ImportParams();
  void Init();
  void UpdateTimeStep();
  void Solve();
  void PrintStatus();
  void DumpSolution(int timestep);

  void UpdateSolution(const Teuchos::RCP<Epetra_Vector>& v);
  void UpdateSolution(const Teuchos::RCP<Epetra_Vector>& v, double Target_Resid);
  void GetSolution(const Teuchos::RCP<Epetra_Vector>& X_out);
//  Epetra_Vector GetSolution();
  double ComputeResidualNorm();

  void ComputeResidual(const Teuchos::RCP<Epetra_Vector>& v);
  
 private:

  void BuildCombo();

  Epetra_Comm* comm;
  Teuchos::RCP<SimParams> simParams;
  Teuchos::RCP<LowOrderProblem> problem;

  Teuchos::RCP<NOX::Epetra::Vector> Q;

  Teuchos::RCP<NOX::Epetra::Interface::Required>       iReq;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec;

  Teuchos::RCP<PhysBasedPrec> precond;

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linearSolver;
  Teuchos::RCP<NOX::Epetra::Group> group;
  Teuchos::RCP<NOX::StatusTest::Combo> combo;
  Teuchos::RCP<NOX::Solver::Generic> nonLinearSolver;

  std::string dir;
  double relTol;
  int newtMax;

};

#endif // JFNK_H_
