#ifndef PDE_PROBLEM_1D3VEM_D_H_
#define PDE_PROBLEM_1D3VEM_D_H_

#include "DiscretePDE.h"

class DiscretePDE_1D3Vem_D : public DiscretePDE
{
 public:

    DiscretePDE_1D3Vem_D(const Teuchos::RCP<MapManager> &map_manager,
		      const Teuchos::RCP<SimParams> &_params,
		      Epetra_Comm* _comm);

  void EvaluateResidual(const Epetra_Vector&, const Epetra_Vector&, const Epetra_Vector&, Epetra_Vector&);

 protected:

  void ImportParams();

  
  Epetra_Comm*             comm;
  Teuchos::RCP<SimParams>  simParams;
  Teuchos::RCP<MapManager> mapManager;
  Teuchos::RCP<StagMesh>   stagMesh;

  Teuchos::RCP<DiffMath> OvX;
  Teuchos::RCP<DiffMath> OvXold;

  double me_h;
  double mi_h;
  double qe_h;
  double qi_h;
  double xi;
  double zeta;
  double dx;
  double dt;

};
#endif // PDE_PROBLEM_1D3VEM_H_
