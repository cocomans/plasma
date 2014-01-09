#ifndef SMESH_H_
#define SMESH_H_

// Trilinos Includes
#include "NOX_Epetra.H"

// User Includes
#include "SimParams.h"
#include "MapManager.h"

class StagMesh
{
 public:
  StagMesh(const Teuchos::RCP<MapManager>& map_manager,
	   const Teuchos::RCP<SimParams> &Params,
	   Epetra_Comm* Comm);

  ~StagMesh();

  double GetXLocC(int elem);
  double GetXLocF(int elem);

 private:

  void ImportParams();

  Epetra_Comm*  comm;
  Teuchos::RCP<SimParams> simParams;
  Teuchos::RCP<MapManager> mapManager;

  double x0;
  double xf;
  double dx;

};

#endif // SMESH_H_
