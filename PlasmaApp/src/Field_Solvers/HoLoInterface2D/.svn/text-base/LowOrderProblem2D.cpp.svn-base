#include "LowOrderProblem2D.h"

#include "DiscretePDE2D_1D1Ves.h"
#include "DiscretePDE2D_1D3Vem.h"
#include "DiscretePDE2D_2D3Vem.h"

LowOrderProblem2D::LowOrderProblem2D( const Teuchos::RCP<SimParams> &Params,
				  Epetra_Comm* Comm )
{
  comm       = Comm;
  simParams  = Params;
  printf("Building MapManager\n");
  mapManager = Teuchos::rcp(new MapManager2D(simParams, comm));

  int nSpatial,nVel;
  printf("Getting Params\n");

  nSpatial = simParams->GetParamsPhys()->get<int>("nSpatial");
  nVel = simParams->GetParamsPhys()->get<int>("nVel");
  printf("Building PDEs\n");
  if((nSpatial == 1)&&(nVel == 1))
	  pde     = Teuchos::rcp( new DiscretePDE2D_1D1Ves( simParams, comm ) );
  else if((nSpatial == 1)&&(nVel == 3))
	  pde     = Teuchos::rcp( new DiscretePDE2D_1D3Vem( simParams, comm ) );
  else if((nSpatial == 2)&&(nVel == 3))
	  pde     = Teuchos::rcp( new DiscretePDE2D_2D3Vem( simParams, comm ) );


  printf("Building Precon\n");
  precond = Teuchos::rcp(new PhysBasedPrec2D( mapManager, simParams, comm ) );
  printf("Building Vectors\n");
  X = Teuchos::rcp( new EpVecWrapper((mapManager),simParams) );
  Xold = Teuchos::rcp( new EpVecWrapper((mapManager),simParams) );
  ConsTerms = Teuchos::rcp( new EpVecWrapper((mapManager),simParams) );
  Res = Teuchos::rcp( new EpVecWrapper((mapManager),simParams) );

  
}

LowOrderProblem2D::~LowOrderProblem2D(){}

void LowOrderProblem2D::Init()
{


	X -> X -> PutScalar(1.0);
	Xold -> X -> PutScalar(1.0);
	ConsTerms-> X -> PutScalar(0.0);
}


bool LowOrderProblem2D::computeF (const Epetra_Vector &x,
				Epetra_Vector &F,
				const FillType fillFlag)
{

	(*X) = x;
	(*Res) = F;

//	EpVecWrapper xt(x,simParams);
//	EpVecWrapper Ft(F,simParams);
  //Edit this function in DiscretePDE.cpp
  pde->EvaluateResidual(*X, *(Xold),  *(ConsTerms), *Res);

  F = *(Res->X);
  return true;
}

bool LowOrderProblem2D::computePreconditioner(const Epetra_Vector &x,
				       Epetra_Operator &M,
				       Teuchos::ParameterList *precParams)
{
//  PhysBasedPrec2D  *Prec = dynamic_cast<PhysBasedPrec2D*>(&M);
//  Prec->UpdateMatrix(x);
  return true;
}


Teuchos::RCP<Epetra_Vector> LowOrderProblem2D::GetX()
{
  return X->X;
}

Teuchos::RCP<Epetra_Vector> LowOrderProblem2D::GetXold()
{
  return Xold->X;
}

Teuchos::RCP<EpVecWrapper> LowOrderProblem2D::GetXW()
{
  return X;
}

Teuchos::RCP<EpVecWrapper> LowOrderProblem2D::GetXoldW()
{
  return Xold;
}

void LowOrderProblem2D::SetXold(const Epetra_Vector &xold)
{
  *Xold = xold;
}

Teuchos::RCP<Epetra_Vector> LowOrderProblem2D::GetConsTerms()
{
  return ConsTerms->X;
}

Teuchos::RCP<EpVecWrapper> LowOrderProblem2D::GetConsTermsW()
{
  return ConsTerms;
}

void LowOrderProblem2D::SetConsTerms(const Epetra_Vector &con)
{
  *ConsTerms = con;
}

Teuchos::RCP<PhysBasedPrec2D> LowOrderProblem2D::GetPreconditioner()
{
  return precond;
}

Teuchos::RCP<MapManager2D> LowOrderProblem2D::GetMapManager()
{
  return mapManager;
}




