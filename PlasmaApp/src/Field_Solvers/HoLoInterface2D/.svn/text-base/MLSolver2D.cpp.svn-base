#include "MLSolver2D.h"


using namespace ML_Epetra;

MLSolver2D::MLSolver2D(const Teuchos::RCP<Epetra_CrsMatrix>& M,
		   const Teuchos::RCP<MapManager2D> &map_manager,
		   const Teuchos::RCP<SimParams> &Params,
		   Epetra_Comm* Comm)

{
  comm       = Comm;
  simParams  = Params;
  mapManager = map_manager;

   simParams->GetParamsML()->set("smoother: user-defined function", userSmoother);
   simParams->GetParamsML()->set("smoother: type", "user-defined");

  mlPrec =  new MultiLevelPreconditioner( *M, *(simParams->GetParamsML()) ) ;
  
}

MLSolver2D::~MLSolver2D()
{
  //delete mlPrec;
}

void MLSolver2D::Solve(const Teuchos::RCP<Epetra_Vector>& ResRed , const Teuchos::RCP<Epetra_Vector>& SolRed)
{
  //  std::cout << "MLSolver::Solve(...)" << std::endl;
   mlPrec->ApplyInverse(*ResRed, *SolRed);

}

void MLSolver2D::ComputePrec()
{
  mlPrec->ReComputePreconditioner();
}







//int userSmoother(ML_Smoother *data, int x_length, double x[],
//                 int rhs_length, double rhs[])
//{
////  std::cout << "userSmoother(...)" << std::endl;
//
//  int i;
//  double *ap, omega = .6666; /* temp vector and damping factor */
//  double *diag;
//  ML_Operator *Amat;
//  ML_Smoother *smoo;
//
//  smoo    = (ML_Smoother *) data;
//  Amat = (ML_Operator *) ML_Get_MySmootherData(smoo);
//  ap = (double *) ML_allocate(Amat->outvec_leng * sizeof(double));
//  ML_Operator_Apply(Amat, x_length, x, rhs_length, ap);
//  ML_Operator_Get_Diag(Amat, x_length, &diag);
//
//  for (i = 0; i < x_length; i++)
//    {
////      std::cout << "diag = " << diag[i] << std::endl;
//      x[i] = x[i] + omega*(rhs[i] - ap[i])/diag[i];
//    }
//
//  ML_free(ap);
//
//  return 0;
//}
