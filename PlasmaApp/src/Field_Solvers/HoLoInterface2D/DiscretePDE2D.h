#ifndef PDE_PROBLEM2D_H_
#define PDE_PROBLEM2D_H_

// Standard Includes

#include <math.h> 

// Trilinos Includes
#include "NOX_Epetra.H"

// User Includes
#include "../HoLoInterface/SimParams.h"
#include "MapManager2D.h"
#include "EpVecWrapper.h"

class DiscretePDE2D
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

  virtual ~DiscretePDE2D(){};

  virtual void EvaluateResidual(EpVecWrapper& x,
		  const EpVecWrapper& Xold,
		  const EpVecWrapper& ConsTerm,
		  EpVecWrapper& res){};

 protected:

  virtual void ImportParams(){printf("ERROR ImportParams Not Implemented in PDE inherited class\n"); exit(1);};
  
  Epetra_Comm*             comm;
  Teuchos::RCP<SimParams>  simParams;

//  Teuchos::RCP<EpVecWrapper> OvX;
//  Teuchos::RCP<EpVecWrapper> OvX;
//  Teuchos::RCP<EpVecWrapper> OvXold;

  double me_h;
  double mi_h;
  double qe_h;
  double qi_h;
  double xi;
  double dx;
  double dt;

};
#endif // PDE_PROBLEM2D_H_
