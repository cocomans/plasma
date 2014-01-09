#ifndef JFNK_2D_H_
#define JFNK_2D_H_


// User Includes
#include "../HoLoInterface/SimParams.h"
#include "LowOrderProblem2D.h"
#include "EpVecWrapper.h"

class JFNK2D
{
 public:
	JFNK2D( const Teuchos::RCP<LowOrderProblem2D>& problem,
	const Teuchos::RCP<SimParams>& Params,
	Epetra_Comm* Comm );
  ~JFNK2D();

  void ImportParams();
  void Init();
  void UpdateTimeStep();
  void Solve();
  void PrintStatus();
  void DumpSolution(int timestep);

  void UpdateSolution(const Teuchos::RCP<EpVecWrapper>& v);
  void UpdateSolution(const Teuchos::RCP<Epetra_Vector>& v, double Target_Resid);
  void GetSolution(const Teuchos::RCP<EpVecWrapper>& X_out);
//  Epetra_Vector GetSolution();
  double ComputeResidualNorm();

  void ComputeResidual(const Teuchos::RCP<EpVecWrapper>& v);
  
 private:

  void BuildCombo();

  Epetra_Comm* comm;
  Teuchos::RCP<SimParams> simParams;
  Teuchos::RCP<LowOrderProblem2D> problem;

  Teuchos::RCP<NOX::Epetra::Vector> Q;

  Teuchos::RCP<NOX::Epetra::Interface::Required>       iReq;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec;

  Teuchos::RCP<PhysBasedPrec2D> precond;

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linearSolver;
  Teuchos::RCP<NOX::Epetra::Group> group;
  Teuchos::RCP<NOX::StatusTest::Combo> combo;
  Teuchos::RCP<NOX::Solver::Generic> nonLinearSolver;

  std::string dir;
  double relTol;
  int newtMax;

};

#endif // JFNK_H_
