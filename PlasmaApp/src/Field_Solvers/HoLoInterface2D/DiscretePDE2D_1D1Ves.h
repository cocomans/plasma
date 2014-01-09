#ifndef PDE_1D1Ves2D_PROBLEM_H_
#define PDE_1D1Ves2D_PROBLEM_H_

// Standard Includes
#include <math.h> 

// Trilinos Includes
#include "NOX_Epetra.H"

// User Includes
#include "DiscretePDE2D.h"

class DiscretePDE2D_1D1Ves : public DiscretePDE2D
{
 public:

	DiscretePDE2D_1D1Ves(const Teuchos::RCP<SimParams> &_params,
	      Epetra_Comm* _comm);


  ~DiscretePDE2D_1D1Ves(){};

  void EvaluateResidual(EpVecWrapper& x,
		  const EpVecWrapper& Xold,
		  const EpVecWrapper& ConsTerm,
		  EpVecWrapper& res);

 protected:

  void ImportParams();
  
  Epetra_Comm*             comm;
  Teuchos::RCP<SimParams>  simParams;

  double me_h;
  double mi_h;
  double qe_h;
  double qi_h;
  double xi;
  double dx;
  double dt;

};
#endif // PDE_1D1Ves2D_PROBLEM_H_
