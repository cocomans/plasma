/*
 * NonLinearSolver.h
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#ifndef NONLINEARSOLVER_H_
#define NONLINEARSOLVER_H_
#include "SimParams.h"

namespace HoLo
{
	namespace LowOrder
	{

	class Problem;
	class Preconditioner;

	class NonLinearSolver
	{
		NonLinearSolver(Problem* problem,
			SimParams* Params,
			Epetra_Comm* Comm );
		  ~NonLinearSolver();

		  void ImportParams();
		  void Init(Problem* problem);
		  void Solve();
		  void PrintStatus();
		  void DumpSolution(int timestep);

		  void UpdateSolution(const Teuchos::RCP<Epetra_Vector>& v);
		  void GetSolution(const Teuchos::RCP<Epetra_Vector>& X_out);
		//  Epetra_Vector GetSolution();
		  double ComputeResidualNorm();

		  void ComputeResidual(const Teuchos::RCP<Epetra_Vector>& v);

		 private:

		  void BuildCombo();

		  Epetra_Comm* comm;
		  SimParams* simParams;

		  NOX::Epetra::Vector* Q;

		  // These are really copies of problem
		  Problem* problem;
		  NOX::Epetra::Interface::Required*       iReq;
		  NOX::Epetra::Interface::Preconditioner* iPrec;

		  Preconditioner* precond;

		  NOX::Epetra::LinearSystemAztecOO* linearSolver;
		  NOX::Epetra::Group* group;
		  NOX::StatusTest::Combo* combo;
		  NOX::Solver::Generic* nonLinearSolver;

		  std::string dir;
		  double relTol;
		  int newtMax;
	};

	} /* namespace LowOrder */
}
#endif /* NONLINEARSOLVER_H_ */
