#include "PhysBasedPrec2D.h"
PhysBasedPrec2D::PhysBasedPrec2D(const Teuchos::RCP<MapManager2D> &map_manager,
			     const Teuchos::RCP<SimParams> &Params,
			     Epetra_Comm* Comm)
{
  comm = Comm;
  simParams = Params;
  mapManager = map_manager;

  ImportParams();

//  SolOv = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );
//  ResOv = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );
//
//  SolSd = Teuchos::rcp(new Epetra_Vector(*(mapManager->PtSdMap)));
//  ResSd  = Teuchos::rcp(new Epetra_Vector(*(mapManager->PtSdMap)));
//
//  ResRed = Teuchos::rcp(new Epetra_Vector(*(mapManager->RedPtSdMap)));
//  SolRed = Teuchos::rcp(new Epetra_Vector(*(mapManager->RedPtSdMap)));

//  curr_X = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );

  // Fill Matrix with dummy values
//  printf("Putting Problem\n");
//  SolSd->PutScalar(1.0);
//  printf("Building CrsMatrix\n");

//  M = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *(mapManager->RedPtSdMap), 3 ) );
//  UpdateMatrix(*SolSd);
//  printf("Building MLSolver\n");
//  mlSolver = Teuchos::rcp(new MLSolver(M, mapManager, simParams, comm) );
}

PhysBasedPrec2D::~PhysBasedPrec2D()
{
}

void PhysBasedPrec2D::ImportParams()
{


  me_h   = simParams->GetParamsPhys()->get<double>("me_h");
  mi_h   = simParams->GetParamsPhys()->get<double>("mi_h");
  qe_h   = simParams->GetParamsPhys()->get<double>("qe_h");
  qi_h   = simParams->GetParamsPhys()->get<double>("qi_h");

  xi    = simParams->GetParamsPhys()->get<double>("xi");
  dx    = simParams->GetParamsPhys()->get<double>("dx");
  dt    = simParams->GetParamsPhys()->get<double>("dt");

  pre   = simParams->GetParamsFlags()->get<bool>("pre");

}

////////////////////////////////////////////////////////////////////
// USER: Edit this routine to update matrix M using the current
// newton iterate x;
///////////////////////////////////////////////////////////////////
void PhysBasedPrec2D::UpdateMatrix(const Epetra_Vector &x)
{
	using namespace ES1D1V;

//
//  curr_X->SetVector(x);
//
//  // Build some temp variables
//  int* colIndices = new int[3];
//  double* vals = new double[3];
//  int numEntries = 0;
//  double dx2 = dx*dx;
//
//  // Loop over elements
//  int N = mapManager->LocSdNumElemX;
//  for(int i = 0; i<N; i++)
//    {
//      int j = mapManager->LocElem_to_LocPtSd(i, 0);
//      double ne = curr_X->GetC2F(n_e,i);
//      double Sxxe = curr_X->GetC2F(Sxx_e,i);
//
//      int row = mapManager->RedPtSdMap->GID( mapManager->LocElem_to_LocRedPtSd(i,0) );
//
//      colIndices[0] = mapManager->RedPtOvMap->GID(  mapManager->LocElem_to_LocRedPtOv(i - 1, 0)  );
//      colIndices[1] = mapManager->RedPtOvMap->GID( mapManager->LocElem_to_LocRedPtOv(i    , 0)   );
//      colIndices[2] = mapManager->RedPtOvMap->GID( mapManager->LocElem_to_LocRedPtOv(i + 1, 0)   );
//
//      vals[0] = - (dt*dt*xi*xi*me_h/4.0) * 1.0 / dx2;
//      vals[1] = ( xi*xi*me_h + dt*dt * qe_h*qe_h /4.0) + (dt*dt*xi*xi*me_h*Sxxe/4.0)*2.0/dx2;
//      vals[2] = - (dt*dt*xi*xi*me_h/4.0)*1.0 / dx2;
//
//      numEntries = 3;
//      M->InsertGlobalValues(row, numEntries, vals, colIndices);
//
//    }
//  //  std::cout<< *M;
//  M->FillComplete();
}

////////////////////////////////////////////////////////////////////
// USER: Edit this routine to reduce the residual to a reduced RHS
// for the ML solver.
///////////////////////////////////////////////////////////////////
int PhysBasedPrec2D::ApplyInverse( const Epetra_MultiVector &r, Epetra_MultiVector &z) const
{

  // Build an enum of variable names

	using namespace ES1D1V;

//  // Test flag for preconditioning
//  if(pre)
//    {
//
//      // Construct overlap res and fill z with res values
//      ResOv->SetVector(*r(0));
//
//      // Loop over processor local elements
//      int N = mapManager->LocSdNumElemX;
//      for(int i = 0; i<N; i++)
//	{
//	  int j = mapManager->LocElem_to_LocPtSd(i, 0);
//	  double ne = curr_X->GetC2F(n_e,i);
//	  double Sxxe = curr_X->GetC2F(Sxx_e,i);
//
//	  // Construct the reduced residual
//	  int j_red = mapManager->LocElem_to_LocRedPtSd(i,0);
//	  (*ResRed)[j_red + pe_red] = xi*xi*(ResOv->Get(p_e, i) - dt*me_h*Sxxe/2.0*ResOv->DxC2F(n_e,i,dx)) + dt*qe_h*ne/2.0*ResOv->Get(Ex,i);
//	}
//
//      // Call ML on reduced RHS
//      mlSolver->ComputePrec();
//      mlSolver->Solve(ResRed, SolRed);
//
//      // Propogate the reduced solution back into standard vector
//      z(0)->Import(*SolRed, *mapManager->RedPtSd_to_PtSd, Insert);
//
//      // Put z into a DiffMath wrapper
//      SolOv->SetVector(*z(0));
//
//      // Back substitute variables
//      for(int i = 0; i<N; i++)
//	{
//	  int j = mapManager->LocElem_to_LocPtSd(i, 0);
//	  // Variable 1 is dependent
//	  (*z(0))[j + n_e] = ResOv->Get(n_e, i) - dt/2.0*SolOv->DxF2C(p_e, i, dx);
//	  (*z(0))[j + Ex]   = ( ResOv->Get(Ex, i) - (dt*qe_h)/(2.0)*SolOv->Get(p_e, i) ) / (xi*xi);
//	}
//    }
  return 0;
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
// Epetra_Operator Inherited
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
const Epetra_Comm& PhysBasedPrec2D::Comm () const{  return *comm;}


const Epetra_Map& PhysBasedPrec2D::OperatorDomainMap () const
{
  return *mapManager->PtSdMap;
}


const Epetra_Map& PhysBasedPrec2D::OperatorRangeMap () const
{
  return *mapManager->PtSdMap;

}

int PhysBasedPrec2D::SetUseTranspose(bool UseTranspose)
{
  return 1;
}

int PhysBasedPrec2D::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return 1;
}


double PhysBasedPrec2D::NormInf() const
{
  return 1;
}

const char* PhysBasedPrec2D::Label() const
{
  return "User Preconditioner";
}

bool PhysBasedPrec2D::UseTranspose() const
{
  return false;
}

bool PhysBasedPrec2D::HasNormInf() const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// SCRATCH WORK BELOW THIS LINE
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


      // Teuchos::RCP<Epetra_LinearProblem> Problem = Teuchos::rcp(new Epetra_LinearProblem(M.get(), SolRed.get(), ResRed.get()) );
      // Teuchos::RCP<AztecOO> Solver = Teuchos::rcp(new AztecOO( *Problem) );
      // Solver->SetAztecOption( AZ_precond, AZ_none );
      // Solver->SetAztecOption(AZ_solver, AZ_gmres);
      // Solver->Iterate(100,1E-9);

      // std::cout << "b = \n";
      // std::cout << *ResRed;
      // std::cout << "x = \n";
      // std::cout << *SolRed;
      // put solver here
      // M->Apply(*SolRed, *ResRed);
      // std::cout << "M*x = \n";
      // std::cout << *ResRed;


  // // test matrix
  //  std::cout << *M;
  // Epetra_Vector x(*mapManager->RedPtSdMap);
  // Epetra_Vector y(*mapManager->RedPtSdMap);
  // x.PutScalar(1.0);
  // y.PutScalar(23.0);
  // std::cout << y;
  // M->Apply(x,y);
  // std::cout << y;

// void PhysBasedPrec::InitializeGraph()
// {

//   graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *(mapManager->RedPtSdMap), *(mapManager->RedPtOvMap), 3));
  
//   int N = mapManager->LocSdNumElemX;
//   int* colIndices = new int[3];

//   for(int i = 0; i < N; i++)
//     {
//       int row = mapManager->LocElem_to_LocRedPtSd(i,0);
//       colIndices[0] = mapManager->LocElem_to_LocRedPtOv(i - 1,0);
//       colIndices[1] = mapManager->LocElem_to_LocRedPtOv(i , 0);
//       colIndices[2] = mapManager->LocElem_to_LocRedPtOv(i + 1,0);
//       graph->InsertMyIndices(row, 3, colIndices);
//     }
//   graph->FillComplete();

// }
