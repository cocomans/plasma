#include "MapManager.h"

MapManager::MapManager(const Teuchos::RCP<SimParams> &params, Epetra_Comm* Comm)
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

MapManager::~MapManager()
{
  delete [] ListRedUnk;
}

void MapManager::ImportParameters()
{
  index = simParams->GetParamsMap()->get<int>("Index");
  GlbNumElemX = simParams->GetParamsMap()->get<int>("Number Elements");
  NumElemBufX  = simParams->GetParamsMap()->get<int>("Element Overlap");
  NumUnk       = simParams->GetParamsMap()->get<int>("Number Unknowns");
  RedNumUnk    = simParams->GetParamsMap()->get<int>("Number Reduced Unknowns");
  ListRedUnk = new int[RedNumUnk];
  ///////////////////////////////////////////////////////
  //    WARNING: MUST CHANGE IF
  //    REDUCED VARIABLES CHANGE
  ///////////////////////////////////////////////////////
  ListRedUnk[0] = 1;    //  Array that points to the reduced unknown
  // Example for multiple reduces: ne and phi
  // ListRedUnk[0] = n_e
  // ListRedUnk[1] = phi
}

void MapManager::BuildElementMaps()
{
  // Build Standard Element Map
  ElemSdMap = Teuchos::rcp( new Epetra_Map(GlbNumElemX, index, *comm) );

  // compute derived values
  LocSdNumElemX = ElemSdMap->NumMyElements();
  LocOvNumElemX = LocSdNumElemX  + 2 * NumElemBufX;

  // Build Overlap Element Map
  int* GlobalIdList = new int[LocOvNumElemX];
  int count = 0;

  // If only 1 Proc
  if( NumProc == 1)
    {
      for(int i=0; i<NumElemBufX; i++)
  	{
  	  GlobalIdList[i] = (GlbNumElemX -1) - (NumElemBufX-1) + i;
  	}
      count = 0;
      for(int i=NumElemBufX; i<LocOvNumElemX-NumElemBufX; i++)
  	{
  	  GlobalIdList[i] = count;
  	  count++;
  	}
      count = 0;
      for(int i= LocOvNumElemX - NumElemBufX; i<LocOvNumElemX; i++)
  	{
  	  GlobalIdList[i] = count;
  	  count++;
  	}
    }
  // If More than 1 proc
  else
    {
      // Left boundary
      if(ProcId == 0 )
	{
	  for(int i=0; i<NumElemBufX; i++)
	    {
	      GlobalIdList[i] = (GlbNumElemX -1) - (NumElemBufX-1) + i;
	    }
	  count = 0;
	  for(int i=NumElemBufX; i<LocOvNumElemX; i++)
	    {
	      GlobalIdList[i] = count;
	      count++;
	    }
	}
      // Internal Procs
      else if(ProcId != 0 && ProcId != NumProc - 1)
	{
	  for(int i=0; i<LocOvNumElemX; i++)
	    {
	      GlobalIdList[i] = ElemSdMap->MinMyGID() - NumElemBufX + i;
	    }
	}
      // Right boundary
      else if(ProcId == NumProc - 1 && NumProc != 0)
	{
	  for(int i=0; i<LocOvNumElemX - NumElemBufX; i++)
	    {
	      GlobalIdList[i] = ElemSdMap->MinMyGID() - NumElemBufX + i;
	    }
	  count = 0;
	  for(int i= LocOvNumElemX - NumElemBufX; i<LocOvNumElemX; i++)
	    {
	      GlobalIdList[i] = count;
	      count++;
	    }
	}
    }
  
  ElemOvMap = Teuchos::rcp( new Epetra_Map(-1, LocOvNumElemX, GlobalIdList, index, *comm) );

  // std::cout << *ElemSdMap;
  // std::cout << *ElemOvMap;
}

void MapManager::BuildPointMaps()
{
  // compute derived values
  GlbNumPtX   = GlbNumElemX * NumUnk;
  LocSdNumPtX = LocSdNumElemX * NumUnk;
  LocOvNumPtX = LocOvNumElemX * NumUnk;
  NumPtBufX   = NumElemBufX * NumUnk;
  
  // Build the Standard Point Map
  PtSdMap = Teuchos::rcp( new Epetra_Map(-1, LocSdNumPtX, index, *comm) );
  
  // Build the Overlap Point Map
  int* GlobalIdList = new int[LocOvNumPtX];
  int count = 0;

  // If only 1 Proc
  if( NumProc == 1)
    {
      for(int i=0; i<NumPtBufX; i++)
  	{
  	  GlobalIdList[i] = (GlbNumPtX -1) - (NumPtBufX-1) + i;
  	}
      count = 0;
      for(int i=NumPtBufX; i<LocOvNumPtX-NumPtBufX; i++)
  	{
  	  GlobalIdList[i] = count;
  	  count++;
  	}
      count = 0;
      for(int i= LocOvNumPtX - NumPtBufX; i<LocOvNumPtX; i++)
  	{
  	  GlobalIdList[i] = count;
  	  count++;
  	}
    }
  // If More than 1 proc
  else
    {
      // Left boundary
      if(ProcId == 0 )
	{
	  for(int i=0; i<NumPtBufX; i++)
	    {
	      GlobalIdList[i] = (GlbNumPtX -1) - (NumPtBufX-1) + i;
	    }
	  count = 0;
	  for(int i=NumPtBufX; i<LocOvNumPtX; i++)
	    {
	      GlobalIdList[i] = count;
	      count++;
	    }
	}
      // Internal Procs
      else if(ProcId != 0 && ProcId != NumProc - 1)
	{
	  for(int i=0; i<LocOvNumPtX; i++)
	    {
	      GlobalIdList[i] = PtSdMap->MinMyGID() - NumPtBufX + i;
	    }
	}
      // Right boundary
      else if(ProcId == NumProc - 1 && NumProc != 0)
	{
	  for(int i=0; i<LocOvNumPtX - NumPtBufX; i++)
	    {
	      GlobalIdList[i] = PtSdMap->MinMyGID() - NumPtBufX + i;
	    }
	  count = 0;
	  for(int i= LocOvNumPtX - NumPtBufX; i<LocOvNumPtX; i++)
	    {
	      GlobalIdList[i] = count;
	      count++;
	    }
	}
    }
  
  PtOvMap = Teuchos::rcp( new Epetra_Map(-1, LocOvNumPtX, GlobalIdList, index, *comm) );
  
  //std::cout << *PtSdMap;
  //std::cout << *PtOvMap;
   
}

void MapManager::BuildReducedMaps()
{
  // compute derived values
  //  int  GlbRedNumPtX = GlbNumElemX * RedNumUnk;
  int  LocSdRedNumPtX = LocSdNumElemX * RedNumUnk;
  int  LocOvRedNumPtX = LocOvNumElemX * RedNumUnk;
  //  int  RedNumPtBufX = NumElemBufX * RedNumUnk;

  // Build Standard Reduced Point Map
  int* GlobalIdList = new int[LocSdRedNumPtX];
  int count = 0;
  
  count = 0;
  for(int i=0; i<LocSdNumElemX; i++)
    {
      int firstPoint = ElemSdMap->GID(i) * NumUnk;
      for(int k = 0; k<RedNumUnk; k++)
	{
	  GlobalIdList[count] = firstPoint + ListRedUnk[k];
	  count++;
	}
    }

  RedPtSdMap = Teuchos::rcp( new Epetra_Map(-1, LocSdRedNumPtX, GlobalIdList, index, *comm) );
  
  // Build Overlab Reduced Point Map
  delete [] GlobalIdList;
  GlobalIdList = new int[LocOvRedNumPtX];

  count = 0;
  for(int i=0; i<LocOvNumElemX; i++)
    {
      int firstPoint = ElemOvMap->GID(i) * NumUnk;
      for(int k = 0; k<RedNumUnk; k++)
	{
	  GlobalIdList[count] = firstPoint + ListRedUnk[k];
	  count++;
	}
    }
  RedPtOvMap = Teuchos::rcp( new Epetra_Map(-1, LocOvRedNumPtX, GlobalIdList, index, *comm) );
  
  //  std::cout << *RedPtOvMap;

}

void MapManager::BuildImporters()
{

  PtSd_to_PtSd       = Teuchos::rcp(new Epetra_Import(*PtSdMap, *PtSdMap) );

  PtSd_to_PtOv       = Teuchos::rcp(new Epetra_Import(*PtOvMap, *PtSdMap) );

  PtSd_to_RedPtSd    = Teuchos::rcp(new Epetra_Import(*RedPtSdMap, *PtSdMap) );

  RedPtSd_to_PtSd    = Teuchos::rcp(new Epetra_Import(*PtSdMap, *RedPtSdMap) );

  PtOv_to_RedPtSd    = Teuchos::rcp(new Epetra_Import(*RedPtSdMap, *PtOvMap) );

  RedPtSd_to_RedPtOv = Teuchos::rcp(new Epetra_Import(*RedPtOvMap, *RedPtSdMap) );

}

Teuchos::RCP<Epetra_Map> MapManager::GetElemSdMap(){return ElemSdMap;}
Teuchos::RCP<Epetra_Map> MapManager::GetElemOvMap(){return ElemOvMap;}
Teuchos::RCP<Epetra_Map> MapManager::GetPtSdMap(){return PtSdMap;}
Teuchos::RCP<Epetra_Map> MapManager::GetPtOvMap(){return PtOvMap;}
Teuchos::RCP<Epetra_Map> MapManager::GetRedPtSdMap(){return RedPtSdMap;}
Teuchos::RCP<Epetra_Map> MapManager::GetRedPtOvMap(){return RedPtOvMap;}


Teuchos::RCP<Epetra_Import> MapManager::Get_PtSd_to_PtOv(){return PtSd_to_PtOv;}
Teuchos::RCP<Epetra_Import> MapManager::Get_PtSt_to_RedPtSd(){return PtSd_to_RedPtSd;}
Teuchos::RCP<Epetra_Import> MapManager::Get_RedPtSd_to_RedPtOv(){return RedPtSd_to_RedPtOv;}


// Index Functions
int MapManager::LocElem_to_LocPtSd(int i, int t)
{
  return i * NumUnk + t;
}

int MapManager::LocElem_to_LocPtOv(int i, int t)
{
  return (i + NumElemBufX) * NumUnk + t;
}

int MapManager::LocElem_to_LocRedPtSd(int i, int t)
{
  return i * RedNumUnk + t;
}

int MapManager::LocElem_to_LocRedPtOv(int i, int t)
{
  return (i + NumElemBufX) * RedNumUnk + t;
}
