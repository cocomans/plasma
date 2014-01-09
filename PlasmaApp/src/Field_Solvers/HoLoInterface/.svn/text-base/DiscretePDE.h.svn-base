#ifndef PDE_PROBLEM_H_
#define PDE_PROBLEM_H_

// Standard Includes

#include <math.h> 

// Trilinos Includes
#include "NOX_Epetra.H"

// User Includes
#include "SimParams.h"
#include "DiffMath.h"
#include "StagMesh.h"
#include "MapManager.h"

class DiscretePDE
{
 public:

//  DiscretePDE(const Teuchos::RCP<MapManager> &map_manager,
//	      const Teuchos::RCP<SimParams> &_params,
//	      Epetra_Comm* _comm)
//	      {
//			  comm       = _comm;
//			  simParams  = _params;
//			  mapManager = map_manager;
//
//			  stagMesh   = Teuchos::rcp(new StagMesh(mapManager, simParams, comm) );
//			  OvX        = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );
//			  OvXold     = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );
//
//			  ImportParams();
//	      };

  virtual ~DiscretePDE(){};

  virtual void EvaluateResidual(const Epetra_Vector&, const Epetra_Vector&, const Epetra_Vector&, Epetra_Vector&){};

 protected:

  virtual void ImportParams(){printf("ERROR ImportParams Not Implemented in PDE inherited class\n"); exit(1);};
  
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
#endif // PDE_PROBLEM_H_
