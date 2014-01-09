#ifndef PDE_1D1Ves_PROBLEM_H_
#define PDE_1D1Ves_PROBLEM_H_

// Standard Includes
#include <math.h> 

// Trilinos Includes
#include "NOX_Epetra.H"

// User Includes
#include "SimParams.h"
#include "DiffMath.h"
#include "StagMesh.h"
#include "MapManager.h"
#include "DiscretePDE.h"

class DiscretePDE_1D1Ves : public DiscretePDE
{
 public:

	DiscretePDE_1D1Ves(const Teuchos::RCP<MapManager> &map_manager,
	      const Teuchos::RCP<SimParams> &_params,
	      Epetra_Comm* _comm);


  ~DiscretePDE_1D1Ves(){};

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
  double dx;
  double dt;

};
#endif // PDE_1D1Ves_PROBLEM_H_
