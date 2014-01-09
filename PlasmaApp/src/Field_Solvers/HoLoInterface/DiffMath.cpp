#include "DiffMath.h"

DiffMath::DiffMath( const Teuchos::RCP<MapManager>& map_manager, const Teuchos::RCP<SimParams> &Params, Epetra_Comm* Comm )
{
  comm       = Comm;
  simParams  = Params;
  mapManager = map_manager;
  X          = Teuchos::rcp( new Epetra_Vector(*(mapManager->GetPtOvMap() )));
  avg        = Teuchos::rcp( new Epetra_Vector(*(mapManager->GetRedPtSdMap() )));
}

DiffMath::~DiffMath(){}

void DiffMath::SetVector( const Epetra_Vector& source_vec )
{
  X->Import(source_vec, *mapManager->Get_PtSd_to_PtOv(), Insert);
}

double& DiffMath::Get(int t, int i)
{
  //  int W_loc = mapManager->LocElem_to_LocPtOv(i - 1, t);
  int C_loc = mapManager->LocElem_to_LocPtOv(i, t);
  //  int E_loc = mapManager->LocElem_to_LocPtOv(i + 1, t);
  return (*X)[C_loc];
}

double DiffMath::GetSM(int t, int i)
{
  int W_loc = mapManager->LocElem_to_LocPtOv(i - 1, t);
  int C_loc = mapManager->LocElem_to_LocPtOv(i    , t);
  int E_loc = mapManager->LocElem_to_LocPtOv(i + 1, t);
  return 0.25*( (*X)[E_loc] + 2.0*(*X)[C_loc] + (*X)[W_loc] );
}

double DiffMath::GetC2F(int t, int i)
{
  int W_loc = mapManager->LocElem_to_LocPtOv(i - 1, t);
  int C_loc = mapManager->LocElem_to_LocPtOv(i    , t);
  //  int E_loc = mapManager->LocElem_to_LocPtOv(i + 1, t);
  return 0.5*( (*X)[C_loc] + (*X)[W_loc] );
}

double DiffMath::GetF2C(int t, int i)
{
  //  int W_loc = mapManager->LocElem_to_LocPtOv(i - 1, t);
  int C_loc = mapManager->LocElem_to_LocPtOv(i    , t);
  int E_loc = mapManager->LocElem_to_LocPtOv(i + 1, t);
  return 0.5*( (*X)[C_loc] + (*X)[E_loc] );
}

double DiffMath::Dx(int t, int i, double dx)
{
  int W_loc = mapManager->LocElem_to_LocPtOv(i - 1, t);
  //  int C_loc = mapManager->LocElem_to_LocPtOv(i    , t);
  int E_loc = mapManager->LocElem_to_LocPtOv(i + 1, t);
  return ( (*X)[E_loc] - (*X)[W_loc] ) / (2.0*dx);
}

double DiffMath::Dxx(int t, int i, double dx)
{
  int W_loc = mapManager->LocElem_to_LocPtOv(i - 1, t);
  int C_loc = mapManager->LocElem_to_LocPtOv(i    , t);
  int E_loc = mapManager->LocElem_to_LocPtOv(i + 1, t);
  return ( (*X)[E_loc] - 2.0*(*X)[C_loc] + (*X)[W_loc] ) / (dx*dx);
}

double DiffMath::Dx(double W, double C, double E, double dx)
{
  return ( E - W ) / (2.0*dx);
}

double DiffMath::DxC2F(int t, int i, double dx)
{
  int W_loc = mapManager->LocElem_to_LocPtOv(i - 1, t);
  int C_loc = mapManager->LocElem_to_LocPtOv(i    , t);
  //  int E_loc = mapManager->LocElem_to_LocPtOv(i + 1, t);
  return ( (*X)[C_loc] - (*X)[W_loc] ) / dx;
}

double DiffMath::DxC2F(double W, double C, double E, double dx)
{
  return ( C - W ) / dx;
}

double DiffMath::DxF2C(int t, int i, double dx)
{
  //  int W_loc = mapManager->LocElem_to_LocPtOv(i - 1, t);
  int C_loc = mapManager->LocElem_to_LocPtOv(i    , t);
  int E_loc = mapManager->LocElem_to_LocPtOv(i + 1, t);
  return ( (*X)[E_loc] - (*X)[C_loc] ) / dx;
}

double DiffMath::DxF2C(double W, double C, double E, double dx)
{
  return ( E - C ) / dx;
}


// GLOBAl Operations
double DiffMath::Mean(int t)
{
  
  // Loop over processor local elements
  int N = mapManager->LocSdNumElemX;
  for(int i = 0; i<N; i++)
    {   
      // Construct the reduced residual
      int j_red = mapManager->LocElem_to_LocRedPtSd(i, 0);
      (*avg)[j_red + 0 ] = Get(t, i);
    }

  // std::cout<< *avg << std::endl;
  // std::cout<< *X << std::endl;
  
  double tmp_avg[1];
  avg->MeanValue(tmp_avg);
  double my_avg = tmp_avg[0];

  //  std::cout<< "my val = " << my_avg << std::endl;

  return my_avg;
}
