#ifndef MLSOLVER_H_
#define MLSOLVER_H_

// Trilinos Includes
#include "NOX_Epetra.H"
#include "ml_MultiLevelPreconditioner.h"

// User Includes
#include "SimParams.h"
#include "MapManager.h"

extern int userSmoother(ML_Smoother *data, int x_length, double x[],
			int rhs_length, double rhs[]);

using namespace ML_Epetra;

class MLSolver
{
 public:
  MLSolver( const Teuchos::RCP<Epetra_CrsMatrix>& M,
	    const Teuchos::RCP<MapManager> &map_manager,
	    const Teuchos::RCP<SimParams> &Params,
	    Epetra_Comm* Comm);
  ~MLSolver();

  void Solve(const Teuchos::RCP<Epetra_Vector>& ResRed , const Teuchos::RCP<Epetra_Vector>& SolRed);
  void ComputePrec();

  /* int userSmoother(ML_Smoother *data, int x_length, double x[], */
  /* 		   int rhs_length, double rhs[]); */

 private:
  
  Epetra_Comm*             comm;
  Teuchos::RCP<SimParams>  simParams;
  Teuchos::RCP<MapManager> mapManager;

  MultiLevelPreconditioner* mlPrec;

};

#endif // MLSOLVER_H_
