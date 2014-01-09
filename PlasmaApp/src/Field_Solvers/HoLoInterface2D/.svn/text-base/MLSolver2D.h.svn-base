#ifndef MLSOLVER_2D_H_
#define MLSOLVER_2D_H_

// Trilinos Includes
#include "NOX_Epetra.H"
#include "ml_MultiLevelPreconditioner.h"

// User Includes
#include "../HoLoInterface/SimParams.h"
#include "MapManager2D.h"

extern int userSmoother(ML_Smoother *data, int x_length, double x[],
			int rhs_length, double rhs[]);

using namespace ML_Epetra;

class MLSolver2D
{
 public:
  MLSolver2D( const Teuchos::RCP<Epetra_CrsMatrix>& M,
	    const Teuchos::RCP<MapManager2D> &map_manager,
	    const Teuchos::RCP<SimParams> &Params,
	    Epetra_Comm* Comm);
  ~MLSolver2D();

  void Solve(const Teuchos::RCP<Epetra_Vector>& ResRed , const Teuchos::RCP<Epetra_Vector>& SolRed);
  void ComputePrec();

  /* int userSmoother(ML_Smoother *data, int x_length, double x[], */
  /* 		   int rhs_length, double rhs[]); */

 private:
  
  Epetra_Comm*             comm;
  Teuchos::RCP<SimParams>  simParams;
  Teuchos::RCP<MapManager2D> mapManager;

  MultiLevelPreconditioner* mlPrec;

};

#endif // MLSOLVER_2D_H_
