/*
 * NonLinearSolver.cpp
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#include "NonLinearSolver.h"
#include "Problem.h"
#include "Preconditioner.h"
namespace HoLo{
namespace LowOrder
{
NonLinearSolver::NonLinearSolver(Problem* _problem,
								SimParams* Params,
								Epetra_Comm* Comm )
{
  comm      = Comm;
  simParams = Params;
  problem   = _problem;

  ImportParams();

}

void NonLinearSolver::ImportParams()
{
  dir = simParams->GetParamsFlags()->get<std::string>("Output Directory");

  relTol = simParams->GetParamsPhys()->get<double>("Newton Relative Tolerance");
  newtMax = simParams->GetParamsPhys()->get<int>("Newton Max Iterations");
}

void NonLinearSolver::Init(Problem* problem)
{
	printf("jfnk init\n");

  // Create a clone of the solution vector
  Q =  new NOX::Epetra::Vector( Teuchos::rcp(problem->GetX()),
					     NOX::Epetra::Vector::CreateCopy,
					     NOX::DeepCopy );

  printf("jfnk Q\n");
  // Construct interface pointers
  iReq  = problem;
  iPrec = problem;
  precond = problem->GetPreconditioner();

  // Construct a linear solver
  linearSolver =
     new NOX::Epetra::LinearSystemAztecOO( simParams->GetPrintParams(),
							simParams->GetLSParams(),
							Teuchos::rcp(iReq),
							Teuchos::rcp(iPrec),
							Teuchos::rcp(precond),
  							*Teuchos::rcp(Q) ) ;

  // Create a group
  group =
    new NOX::Epetra::Group(
					*(simParams->GetParamsNLS()),
					Teuchos::rcp(iReq),
					*Teuchos::rcp(Q),
					Teuchos::rcp(linearSolver));
  printf("jfnk group\n");

  // Set up stopping criteria for NOX
  BuildCombo();

  printf("jfnk buildcombo\n");


  // Create the JFNK method
  nonLinearSolver = &(*(NOX::Solver::buildSolver(Teuchos::rcp( group),
		  Teuchos::rcp(combo),
					      simParams->GetParamsNLS() )));

  printf("jfnk nonlin\n");


}

void NonLinearSolver::Solve()
{
   nonLinearSolver->solve();
}


void NonLinearSolver::UpdateSolution(const Teuchos::RCP<Epetra_Vector>& v)
{

   Epetra_Vector& tmp = Q->getEpetraVector();
   tmp = *v;
//   tmp.PutScalar(1.0);
  nonLinearSolver->reset(*Q);

  printf("finished solution update jfnk\n");

}





void NonLinearSolver::GetSolution(const Teuchos::RCP<Epetra_Vector>& X_out)
{
  // Extract solution from group
  const NOX::Epetra::Group &finalGroup =
    dynamic_cast<const NOX::Epetra::Group&> (nonLinearSolver->getSolutionGroup());

  (*X_out) =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

//  return finalSolution;

}

double NonLinearSolver::ComputeResidualNorm()
{

  // Extract solution from group
  const NOX::Epetra::Group* finalGroup =
    dynamic_cast<const NOX::Epetra::Group*> (&(nonLinearSolver->getSolutionGroup()));
  NOX::Epetra::Group* group =  const_cast<NOX::Epetra::Group*>(finalGroup);

  group->computeF();


  printf("finished redisd norm calc jfnk\n");


  return (nonLinearSolver->getSolutionGroup()).getNormF();

}

void NonLinearSolver::ComputeResidual(const Teuchos::RCP<Epetra_Vector>& v)
{
  //  nonLinearSolver->getSolutionGroup().computeF();
  // I'm sorry for who ever has to look at this line...
  bool t = (nonLinearSolver->getSolutionGroup()).isF();
  std::cout << "isF() = " << t << std::endl;

  const NOX::Epetra::Group* finalGroup =
    dynamic_cast<const NOX::Epetra::Group*> (&(nonLinearSolver->getSolutionGroup()));
  NOX::Epetra::Group* group =  const_cast<NOX::Epetra::Group*>(finalGroup);
  std::cout << "ComputeF from Residual" << std::endl;

  group->computeF();

  *v = (dynamic_cast<const NOX::Epetra::Vector&>(group->getF() )).getEpetraVector();

//  std::cout << *v;
}

void NonLinearSolver::PrintStatus()
{
  // Print final status
//  if(comm->MyPID() == 0)
//    std::cout<< (nonLinearSolver->getList()).sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output");

}

void NonLinearSolver::DumpSolution(int timestep)
{
  // Extract solution
  const NOX::Epetra::Group &finalGroup =
    dynamic_cast<const NOX::Epetra::Group&> (nonLinearSolver->getSolutionGroup());

  const Epetra_Vector &finalSolution =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

 // Dump to files

//  char file_name[25];
//  sprintf(file_name, "./%s/output_pid_%d_ts_%d", dir.c_str(), comm->MyPID(), timestep);
//  std::ofstream filestream ( file_name );
//  filestream<< setprecision(16);
//  filestream << finalSolution;
//  filestream.close();

//  std::cout << finalSolution;

  // std::cout<< (nonLinearSolver->getList()).sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").get<int>("Total Number of Linear Iterations");


}


void NonLinearSolver::BuildCombo()
{
//
//  // Build a combo of stopping criteria in OR mode
  combo = new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR);

  Teuchos::RCP<NOX::StatusTest::Combo> tolCombo =
		  Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  printf("Relative JFNK TOl = %e\n",relTol);

  Teuchos::RCP<NOX::StatusTest::NormUpdate> relResUpdate =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate( 1.0e-5) );

  Teuchos::RCP<NOX::StatusTest::NormF> relRes =
    Teuchos::rcp(new NOX::StatusTest::NormF(*group, relTol) );

  tolCombo->addStatusTest(relResUpdate);
  tolCombo->addStatusTest(relRes);



  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(newtMax));



  combo->addStatusTest(tolCombo);

  combo->addStatusTest(maxiters);

}



} /* namespace LowOrder */
} /* namespace HoLo */
