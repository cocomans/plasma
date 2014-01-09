#ifndef PDE_PROBLEM2D_1D3VEM_H_
#define PDE_PROBLEM2D_1D3VEM_H_

#include "DiscretePDE2D.h"

class DiscretePDE2D_1D3Vem : public DiscretePDE2D
{
 public:

    DiscretePDE2D_1D3Vem(const Teuchos::RCP<SimParams> &_params,
		      Epetra_Comm* _comm);

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
  double zeta;
  double dx;
  double dt;

};
#endif // PDE_PROBLEM2D_1D3VEM_H_
