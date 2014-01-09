#include "MLSolver.h"

using namespace ML_Epetra;

MLSolver::MLSolver(const Teuchos::RCP<Epetra_CrsMatrix>& M,
                   const Teuchos::RCP<SimParams> &Params,
                   Epetra_Comm* Comm)

{
    comm       = Comm;
    simParams  = Params;

    mlPrec =  new MultiLevelPreconditioner( *M, *(simParams->GetParamsML()) ) ;
}

MLSolver::~MLSolver(){}

void MLSolver::Solve(const Teuchos::RCP<Epetra_Vector>& ResRed , const Teuchos::RCP<Epetra_Vector>& SolRed)
{
  mlPrec->ApplyInverse(*ResRed, *SolRed);
}

void MLSolver::ComputePrec()
{
  mlPrec->ReComputePreconditioner();
}
