#ifndef MAPMANAGER_H_
#define MAPMANAGER_H_

//Trilinos Includes
#include "NOX_Epetra.H"
#include <Teuchos_FancyOStream.hpp>

//User Includes
#include "SimParams.h"

class MapManager
{
 public:

  MapManager(const Teuchos::RCP<SimParams> &params, Epetra_Comm* comm);
  ~MapManager();

  Teuchos::RCP<Epetra_Map> GetElemSdMap();
  Teuchos::RCP<Epetra_Map> GetElemOvMap();
  Teuchos::RCP<Epetra_Map> GetPtSdMap();
  Teuchos::RCP<Epetra_Map> GetPtOvMap();
  Teuchos::RCP<Epetra_Map> GetRedPtSdMap();
  Teuchos::RCP<Epetra_Map> GetRedPtOvMap();

  Teuchos::RCP<Epetra_Import> Get_PtSd_to_PtOv();
  Teuchos::RCP<Epetra_Import> Get_PtSt_to_RedPtSd();
  Teuchos::RCP<Epetra_Import> Get_RedPtSd_to_RedPtOv();
  Teuchos::RCP<Epetra_Import> Get_PtOv_to_RedPtSd();

  int LocElem_to_LocPtSd(int i, int t);
  int LocElem_to_LocPtOv(int i, int t);
  int LocElem_to_LocRedPtSd(int i, int t);
  int LocElem_to_LocRedPtOv(int i, int t);

  // Element Variables
  Teuchos::RCP<Epetra_Map> ElemSdMap;
  Teuchos::RCP<Epetra_Map> ElemOvMap;
  int GlbNumElemX;//user supplied
  int NumElemBufX;//user supplied
  int LocSdNumElemX;//derived x
  int LocOvNumElemX;//derived x
  
  // Point Variables
  Teuchos::RCP<Epetra_Map> PtSdMap;
  Teuchos::RCP<Epetra_Map> PtOvMap;
  int NumUnk;//user supplied
  int GlbNumPtX;//derived
  int LocSdNumPtX;//derived
  int LocOvNumPtX;//derived
  int NumPtBufX;//derived

  // Reduced Variables
  Teuchos::RCP<Epetra_Map> RedPtSdMap;
  Teuchos::RCP<Epetra_Map> RedPtOvMap;
  int* ListRedUnk;//user supplied
  int  RedNumUnk;//user supplied
  int  GlbRedNumPtX;//derived
  int  LocSdRedNumPtX;//derived
  int  LocOvRedNumPtX;//derived
  int  RedNumPtBufX;//derived

  // Importers
  Teuchos::RCP<Epetra_Import> PtSd_to_PtSd;
  Teuchos::RCP<Epetra_Import> PtSd_to_PtOv;
  Teuchos::RCP<Epetra_Import> PtSd_to_RedPtSd;
  Teuchos::RCP<Epetra_Import> RedPtSd_to_PtSd;
  Teuchos::RCP<Epetra_Import> RedPtSd_to_RedPtOv;
  Teuchos::RCP<Epetra_Import> PtOv_to_RedPtSd;

 private:

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

#endif // MAPMANAGER_H_
