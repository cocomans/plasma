#ifndef MAPMANAGER2D_H_
#define MAPMANAGER2D_H_

//Trilinos Includes
#include "NOX_Epetra.H"
#include <Teuchos_FancyOStream.hpp>

//User Includes
#include "../HoLoInterface/SimParams.h"

class MapManager2D
{
 public:

	MapManager2D(const Teuchos::RCP<SimParams> &params, Epetra_Comm* comm);
  ~MapManager2D();

  Teuchos::RCP<Epetra_Map> GetElemSdMap();
  Teuchos::RCP<Epetra_Map> GetPtSdMap();

  // Element Variables
  Teuchos::RCP<Epetra_Map> ElemSdMap;
  int GlbNumElemX;//user supplied
  int NumElemBufX;//user supplied
  int LocSdNumElemX;//derived x
  
  // Point Variables
  Teuchos::RCP<Epetra_Map> PtSdMap;
  int NumUnk;//user supplied
  int GlbNumPtX;//derived
//  int LocSdNumPtX;//derived
//  int LocOvNumPtX;//derived
//  int NumPtBufX;//derived

//  int* ListRedUnk;//user supplied
//  int  RedNumUnk;//user supplied
//  int  GlbRedNumPtX;//derived
//  int  LocSdRedNumPtX;//derived
//  int  LocOvRedNumPtX;//derived
//  int  RedNumPtBufX;//derived

 private:

  Teuchos::RCP<Epetra_Import> PtSd_to_PtSd;


  void ImportParameters();
  void BuildElementMaps();
  void BuildPointMaps();
  void BuildReducedMaps();
  void BuildImporters();

  Teuchos::RCP<SimParams> simParams;

  Epetra_Comm* comm;
  int NumProc;
  int ProcId;
  int index;

};

#endif // MAPMANAGER2D_H_
