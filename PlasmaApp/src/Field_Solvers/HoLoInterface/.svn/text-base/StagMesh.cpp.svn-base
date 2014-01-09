#include "StagMesh.h"

StagMesh::StagMesh(const Teuchos::RCP<MapManager>& map_manager,
		   const Teuchos::RCP<SimParams> &Params,
		   Epetra_Comm* Comm)
{
  
  comm = Comm;
  simParams = Params;
  mapManager = map_manager;

  ImportParams();


}

StagMesh::~StagMesh()
{
}

void StagMesh::ImportParams()
{

  x0 = simParams->GetParamsPhys()->get<double>("x0");
  xf = simParams->GetParamsPhys()->get<double>("xf");
  dx = simParams->GetParamsPhys()->get<double>("dx");

}

double StagMesh::GetXLocF(int elem)
{
  int glbElem = mapManager->ElemSdMap->GID(elem);
  return (double)glbElem * dx;
}  

double StagMesh::GetXLocC(int elem)
{
    int glbElem = mapManager->ElemSdMap->GID(elem);
    return ((double)glbElem + 0.5) * dx;
}


// int StagMesh::GetPoint(int elem, int type)
// {
//   return elem*m_Num_Unkwn + type;
// }

// int StagMesh::GetPointOverlap(int elem, int type)
// {
//   return (elem + m_overlap)*m_Num_Unkwn + type;
// }
