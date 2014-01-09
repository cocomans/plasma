#include "MapManager2D.h"

MapManager2D::MapManager2D(const Teuchos::RCP<SimParams> &params, Epetra_Comm* Comm)
{
  simParams = params;
  comm = Comm;
  NumProc = comm->NumProc(); 
  ProcId = comm->MyPID();


  ImportParameters();
  BuildElementMaps();
  BuildPointMaps();
  BuildReducedMaps();
  BuildImporters();
  
}

MapManager2D::~MapManager2D()
{
//  delete [] ListRedUnk;
}

void MapManager2D::ImportParameters()
{
  index = simParams->GetParamsMap()->get<int>("Index");
  GlbNumElemX = simParams->GetParamsMap()->get<int>("Number Elements");
  NumElemBufX  = simParams->GetParamsMap()->get<int>("Element Overlap");
  NumUnk       = simParams->GetParamsMap()->get<int>("Number Unknowns");
//  RedNumUnk    = simParams->GetParamsMap()->get<int>("Number Reduced Unknowns");
//  ListRedUnk = new int[RedNumUnk];
  ///////////////////////////////////////////////////////
  //    WARNING: MUST CHANGE IF
  //    REDUCED VARIABLES CHANGE
  ///////////////////////////////////////////////////////
//  ListRedUnk[0] = 1;    //  Array that points to the reduced unknown
  // Example for multiple reduces: ne and phi
  // ListRedUnk[0] = n_e
  // ListRedUnk[1] = phi
}

void MapManager2D::BuildElementMaps()
{
  // Build Standard Element Map
  ElemSdMap = Teuchos::rcp( new Epetra_Map(GlbNumElemX, index, *comm) );

  // compute derived values
  LocSdNumElemX = ElemSdMap->NumMyElements();
  
  // std::cout << *ElemSdMap;
  // std::cout << *ElemOvMap;
}

void MapManager2D::BuildPointMaps()
{
  // compute derived values
  GlbNumPtX   = (GlbNumElemX) * NumUnk;
//  LocSdNumPtX = LocSdNumElemX * NumUnk;
//  NumPtBufX   = NumElemBufX * NumUnk;
  
  // Build the Standard Point Map
  PtSdMap = Teuchos::rcp( new Epetra_Map(-1,GlbNumPtX, index, *comm) );
  //std::cout << *PtSdMap;
  //std::cout << *PtOvMap;
   
}

void MapManager2D::BuildReducedMaps()
{


}

void MapManager2D::BuildImporters()
{

  PtSd_to_PtSd       = Teuchos::rcp(new Epetra_Import(*PtSdMap, *PtSdMap) );


}


Teuchos::RCP<Epetra_Map> MapManager2D::GetPtSdMap(){return PtSdMap;}



