#include "LowOrderProblem.h"

#include "DiscretePDE_1D1Ves.h"
#include "DiscretePDE_1D3Vem.h"

LowOrderProblem::LowOrderProblem( const Teuchos::RCP<SimParams> &Params,
				  Epetra_Comm* Comm )
{
  comm       = Comm;
  simParams  = Params;
  printf("Building MapManager\n");
  mapManager = Teuchos::rcp(new MapManager(simParams, comm));

  int nSpatial,nVel;
  printf("Getting Params\n");

  nSpatial = simParams->GetParamsPhys()->get<int>("nSpatial");
  nVel = simParams->GetParamsPhys()->get<int>("nVel");
  printf("Building PDEs\n");
  if((nSpatial == 1)&&(nVel == 1))
	  pde     = Teuchos::rcp( new DiscretePDE_1D1Ves( mapManager, simParams, comm ) );
  else if((nSpatial == 1)&&(nVel == 3))
	  pde     = Teuchos::rcp( new DiscretePDE_1D3Vem( mapManager, simParams, comm ) );

  printf("Building Precon\n");
  precond = Teuchos::rcp(new PhysBasedPrec( mapManager, simParams, comm ) );
  printf("Building Vectors\n");
  X = Teuchos::rcp( new Epetra_Vector( *(mapManager->GetPtSdMap() ) ) );
  Xold = Teuchos::rcp( new Epetra_Vector( *(mapManager->GetPtSdMap() ) ) );
  ConsTerms = Teuchos::rcp( new Epetra_Vector( *(mapManager->GetPtSdMap() ) ) );
  
}

LowOrderProblem::~LowOrderProblem(){}

void LowOrderProblem::Init()
{


	X -> PutScalar(1.0);
	Xold ->PutScalar(1.0);
	ConsTerms->PutScalar(0.0);
}


bool LowOrderProblem::computeF (const Epetra_Vector &x,
				Epetra_Vector &F,
				const FillType fillFlag)
{
  //Edit this function in DiscretePDE.cpp
  pde->EvaluateResidual(x, *(Xold),  *(ConsTerms), F);
  return true;
}

bool LowOrderProblem::computePreconditioner(const Epetra_Vector &x,
				       Epetra_Operator &M,
				       Teuchos::ParameterList *precParams)
{
  PhysBasedPrec  *Prec = dynamic_cast<PhysBasedPrec*>(&M);
  Prec->UpdateMatrix(x);
  return true;
}


Teuchos::RCP<Epetra_Vector> LowOrderProblem::GetX()
{
  return X;
}

Teuchos::RCP<Epetra_Vector> LowOrderProblem::GetXold()
{
  return Xold;
}

void LowOrderProblem::SetXold(const Epetra_Vector &xold)
{
  *Xold = xold;
}

Teuchos::RCP<Epetra_Vector> LowOrderProblem::GetConsTerms()
{
  return ConsTerms;
}

void LowOrderProblem::SetConsTerms(const Epetra_Vector &con)
{
  *ConsTerms = con;
}

Teuchos::RCP<PhysBasedPrec> LowOrderProblem::GetPreconditioner()
{
  return precond;
}

Teuchos::RCP<MapManager> LowOrderProblem::GetMapManager()
{
  return mapManager;
}




