#include "DiscretePDE2D_2D3Vem.h"

DiscretePDE2D_2D3Vem::DiscretePDE2D_2D3Vem(const Teuchos::RCP<SimParams> &_params,
	      Epetra_Comm* _comm)
{
  comm       = _comm;
  simParams  = _params;
//  mapManager = map_manager;
//
//  stagMesh   = Teuchos::rcp(new StagMesh(mapManager, simParams, comm) );
//  OvX        = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );
//  OvXold     = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );

  ImportParams();
};

//DiscretePDE::~DiscretePDE(){}

////////////////////////////////////////////////////////////////////
// USER: Edit any constants you may need in SimParams.cpp  //
///////////////////////////////////////////////////////////////////
void DiscretePDE2D_2D3Vem::ImportParams()
{

  me_h   = simParams->GetParamsPhys()->get<double>("me_h");
  mi_h   = simParams->GetParamsPhys()->get<double>("mi_h");
  qe_h   = simParams->GetParamsPhys()->get<double>("qe_h");
  qi_h   = simParams->GetParamsPhys()->get<double>("qi_h");

  zeta = simParams->GetParamsPhys()->get<double>("zeta");
  xi    = simParams->GetParamsPhys()->get<double>("xi");
  dx    = simParams->GetParamsPhys()->get<double>("dx");
  dt    = simParams->GetParamsPhys()->get<double>("dt");
}



////////////////////////////////////////////////////////////////////
// USER: This is where users must build the main residual calc.   //
// Propogate the enum below with names for each of your unknowns. //
// The number of these unknowns must match "Number Unknowns"      //
///////////////////////////////////////////////////////////////////
void DiscretePDE2D_2D3Vem::EvaluateResidual(EpVecWrapper& x,
				    const EpVecWrapper& xold,
				    const EpVecWrapper& gamma,
				    EpVecWrapper& res )
{
  // Build an enum of variable names

	using namespace EM2D3V;


//	printf("Evaluating 2D3V resid\n");

  //Global Operations (ie. averages)
  double jmeanx1 = qe_h*x.mean(px_e) + qi_h*x.mean(px_i);
  double jmeany1 = qe_h*x.mean(py_e) + qi_h*x.mean(py_i);
  double jmeanz1 = qe_h*x.mean(pz_e) + qi_h*x.mean(pz_i);

//  double Te_half = 0.25*me_h*(
//		  x.mean(Sxx_e) + x.mean(Syy_e) + x.mean(Szz_e)
//		  +
//		  xold.mean(Sxx_e) + xold.mean(Syy_e) + xold.mean(Szz_e));
//
//  double Ti_half = 0.25*mi_h*(
//		  x.mean(Sxx_i) + x.mean(Syy_i) + x.mean(Szz_i)
//		  +
//		  xold.mean(Sxx_i) + xold.mean(Syy_i) + xold.mean(Szz_i));

  // Half Time step quantities
//#pragma omp parallel for
//  for(int j = 0; j<x.ny; j++)
//  for(int i = 0; i<x.nx; i++)
//	  x(Phi,i,j) = 0.25*(x(Chi,i,j) + xold(Chi,i,j))*xi*xi + 0.5*xold(Phi,i,j);

#pragma omp parallel for
  for(int j = 0; j<x.ny; j++)
  for(int i = 0; i<x.nx; i++)
  {
//	  x(Ex,i,j) = -x.DxC2F(Phi,i,j) - (x(Ax,i,j) - xold(Ax,i,j))/dt;
//	  x(Ey,i,j) = -x.DyC2F(Phi,i,j) - (x(Ay,i,j) - xold(Ay,i,j))/dt;
//	  x(Ez,i,j) = -(x(Az,i,j) - xold(Az,i,j))/dt;
//
//	  x(Bx,i,j) = x.DyC2V(Az,i,j);
//	  x(By,i,j) = -x.DxC2V(Az,i,j);
//	  x(Bz,i,j) = x.DxC2F(Ay,i,j) - x.DyC2F(Ax,i,j);
//
//	  printf("Ex[%i, %i] = %e\n",i,j,x(Ex,i,j));

	  // Stress Tensors at t+1/2
	  double n0 = xold(n_e,i,j);
	  double n1 = x(n_e,i,j);
	  // Store it in gamma, not being used for anything else...
	  res(Sxx_e,i,j) = 0.5*(n1*x(Sxx_e,i,j) + n0*xold(Sxx_e,i,j));
	  res(Sxy_e,i,j) = 0.5*(n1*x(Sxy_e,i,j) + n0*xold(Sxy_e,i,j));
	  res(Sxz_e,i,j) = 0.5*(n1*x(Sxz_e,i,j) + n0*xold(Sxz_e,i,j));
	  res(Syy_e,i,j) = 0.5*(n1*x(Syy_e,i,j) + n0*xold(Syy_e,i,j));
	  res(Syz_e,i,j) = 0.5*(n1*x(Syz_e,i,j) + n0*xold(Syz_e,i,j));

	  n0 = xold(n_i,i,j);
	  n1 = x(n_i,i,j);
	  // Store it in gamma, not being used for anything else...
	  res(Sxx_i,i,j) = 0.5*(n1*x(Sxx_i,i,j) + n0*xold(Sxx_i,i,j));
	  res(Sxy_i,i,j) = 0.5*(n1*x(Sxy_i,i,j) + n0*xold(Sxy_i,i,j));
	  res(Sxz_i,i,j) = 0.5*(n1*x(Sxz_i,i,j) + n0*xold(Sxz_i,i,j));
	  res(Syy_i,i,j) = 0.5*(n1*x(Syy_i,i,j) + n0*xold(Syy_i,i,j));
	  res(Syz_i,i,j) = 0.5*(n1*x(Syz_i,i,j) + n0*xold(Syz_i,i,j));
  }

  // Phi, Ex, Ey, Ez, Bx, By, Bz are all half-time step quantities for x
  // The stress tensors are now time averaged and multiplied by density


#pragma omp parallel for
	for(int j = 0; j<x.ny; j++)
	for(int i = 0; i<x.nx; i++)
	{

//		// Electron equations
//		res(n_e,i,j) = (x(n_e,i,j) - xold(n_e,i,j))
//				+ dt*(x.DxF2C(px_e,i,j) + x.DyF2C(py_e,i,j))
//				- gamma(n_e,i,j);
//
//		// Forward Euler
//		res(px_e,i,j) = me_h*(x(px_e,i,j) - xold(px_e,i,j)) +
//				dt*(
//					me_h*res.DxC2F(Sxx_e,i,j) +
//					me_h*res.DyC2xF(Sxy_e,i,j) -
//					qe_h*0.5*(
//							0.5*(x.GetC2xF(n_e,i,j) + xold.GetC2xF(n_e,i,j))*(x(Ex,i,j)+xold(Ex,i,j)) +
//							x.GetyF2xF(py_e,i,j)*(x.GetV2xF(Bz,i,j) + xold.GetV2xF(Bz,i,j)) -
//							x.GetC2xF(pz_e,i,j)*(x.GetV2xF(By,i,j)+xold.GetV2xF(By,i,j)))
//
//				) - gamma(px_e,i,j);
//
//		res(py_e,i,j) = me_h*(x(py_e,i,j) - xold(py_e,i,j)) +
//				dt*(
//					me_h*res.DxC2yF(Sxy_e,i,j) +
//					me_h*res.DyC2F(Syy_e,i,j) -
//					qe_h*0.5*(
//							0.5*(x.GetC2yF(n_e,i,j) + xold.GetC2yF(n_e,i,j))*(x(Ey,i,j) + xold(Ey,i,j)) +
//							x.GetC2yF(pz_e,i,j)*(x.GetV2yF(Bx,i,j) + xold.GetV2yF(Bx,i,j)) -
//							x.GetxF2yF(px_e,i,j)*(x.GetV2yF(Bz,i,j) + xold.GetV2yF(Bz,i,j)))
//
//				) - gamma(py_e,i,j);
//
//		res(pz_e,i,j) = me_h*(x(pz_e,i,j) - xold(pz_e,i,j)) +
//				dt*(
//					me_h*res.Dx(Sxz_e,i,j) +
//					me_h*res.Dy(Syz_e,i,j) -
//					qe_h*(
//							0.5*(x.Get(n_e,i,j) + xold.Get(n_e,i,j))*(x(Ez,i,j) + xold(Ez,i,j)) +
//							x.GetxF2C(px_e,i,j)*(x.GetV2C(By,i,j) + xold.GetV2C(By,i,j)) -
//							x.GetyF2C(py_e,i,j)*(x.GetV2C(Bx,i,j) + xold.GetV2C(Bx,i,j)))
//
//				) - gamma(pz_e,i,j);
//
//
//		// Ion equations
//		res(n_i,i,j) = (x(n_i,i,j) - xold(n_i,i,j))
//				+ dt*(x.DxF2C(px_i,i,j) + x.DyF2C(py_i,i,j))
//				- gamma(n_i,i,j);
//
//		// Forward Euler
//		res(px_i,i,j) = mi_h*(x(px_i,i,j) - xold(px_i,i,j)) +
//				dt*(
//					mi_h*res.DxC2F(Sxx_i,i,j) +
//					mi_h*res.DyC2xF(Sxy_i,i,j) -
//					qi_h*0.5*(
//							0.5*(x.GetC2xF(n_i,i,j) + xold.GetC2xF(n_i,i,j))*(x(Ex,i,j)+xold(Ex,i,j)) +
//							x.GetyF2xF(py_i,i,j)*(x.GetV2xF(Bz,i,j) + xold.GetV2xF(Bz,i,j)) -
//							x.GetC2xF(pz_i,i,j)*(x.GetV2xF(By,i,j)+xold.GetV2xF(By,i,j)))
//
//				) - gamma(px_i,i,j);
//
//		res(py_i,i,j) = mi_h*(x(py_i,i,j) - xold(py_i,i,j)) +
//				dt*(
//					mi_h*res.DxC2yF(Sxy_i,i,j) +
//					mi_h*res.DyC2F(Syy_i,i,j) -
//					qi_h*0.5*(
//							0.5*(x.GetC2yF(n_i,i,j) + xold.GetC2yF(n_i,i,j))*(x(Ey,i,j)+xold(Ey,i,j)) +
//							x.GetC2yF(pz_i,i,j)*(x.GetV2yF(Bx,i,j)+xold.GetV2yF(Bx,i,j)) -
//							x.GetxF2yF(px_i,i,j)*(x.GetV2yF(Bz,i,j)+xold.GetV2yF(Bz,i,j)))
//
//				) - gamma(py_i,i,j);
//
//		res(pz_i,i,j) = me_h*(x(pz_i,i,j) - xold(pz_i,i,j)) +
//				dt*(
//					mi_h*res.Dx(Sxz_i,i,j) +
//					mi_h*res.Dy(Syz_i,i,j) -
//					qi_h*0.5*(
//							0.5*(x.Get(n_i,i,j) + xold.Get(n_i,i,j))*(x(Ez,i,j)+xold(Ez,i,j)) +
//							x.GetxF2C(px_i,i,j)*(x.GetV2C(By,i,j)+xold.GetV2C(By,i,j)) -
//							x.GetyF2C(py_i,i,j)*(x.GetV2C(Bx,i,j)+xold.GetV2C(Bx,i,j)))
//
//				) - gamma(pz_i,i,j);


//		// Field Equations
//		res(Chi,i,j) = x.Dxx(Chi,i,j) + x.Dyy(Chi,i,j)
//				- (qe_h*x.DxF2C(px_e,i,j) + qi_h*x.DxF2C(px_i,i,j) +
//				   qe_h*x.DyF2C(py_e,i,j) + qi_h*x.DyF2C(py_i,i,j));
//
//		res(Ax,i,j) = zeta*0.5*(x.Dxx(Ax,i,j) + x.Dyy(Ax,i,j)
//				+ xold.Dxx(Ax,i,j) + xold.Dyy(Ax,i,j))
//				+ (qe_h*x(px_e,i,j) + qi_h*x(px_i,i,j) -jmeanx1
//				- 0.5*(x.DxC2F(Chi,i,j) + xold.DxC2F(Chi,i,j)));
//
//		res(Ay,i,j) = zeta*0.5*(x.Dxx(Ay,i,j) + x.Dyy(Ay,i,j)
//				+ xold.Dxx(Ay,i,j) + xold.Dyy(Ay,i,j))
//				+ (qe_h*x(py_e,i,j) + qi_h*x(py_i,i,j) - jmeany1
//				- 0.5*(x.DyC2F(Chi,i,j) + xold.DyC2F(Chi,i,j)));
//
//		res(Az,i,j) = zeta*0.5*(x.Dxx(Az,i,j) + x.Dyy(Az,i,j)
//				+ xold.Dxx(Az,i,j) + xold.Dyy(Az,i,j))
//				+ (qe_h*x(pz_e,i,j) + qi_h*x(pz_i,i,j)
//				);

		res(Bx,i,j) = 0.5*dt*(x.DyC2V(Ez,i,j)+xold.DyC2V(Ez,i,j)) +
				  	  (x(Bx,i,j) - xold(Bx,i,j));

		res(By,i,j) = -0.5*dt*(x.DxC2V(Ez,i,j)+xold.DxC2V(Ez,i,j)) +
				  	  (x(By,i,j) - xold(By,i,j));

		res(Bz,i,j) = 0.5*dt*(x.DxC2F(Ey,i,j)+xold.DxC2F(Ey,i,j)) -
				      0.5*dt*(x.DyC2F(Ex,i,j)+xold.DyC2F(Ex,i,j)) +
						  	  (x(Bz,i,j) - xold(Bz,i,j));

		res(Ex,i,j) = (x(Ex,i,j) - xold(Ex,i,j)) +
					  dt/(xi*xi)*(
							 (qe_h*x(px_e,i,j) + qi_h*x(px_i,i,j))
							 -
							 0.5*zeta*(x.DyF2C(Bz,i,j)+xold.DyF2C(Bz,i,j))
					  );

		res(Ey,i,j) = (x(Ey,i,j) - xold(Ey,i,j)) +
					  dt/(xi*xi)*(
							 (qe_h*x(py_e,i,j) + qi_h*x(py_i,i,j))
							 +
							 0.5*zeta*(x.DxF2C(Bz,i,j)+xold.DxF2C(Bz,i,j))
					  );

		res(Ez,i,j) = (x(Ez,i,j) - xold(Ez,i,j)) +
					  dt/(xi*xi)*(
							 (qe_h*x(pz_e,i,j) + qi_h*x(pz_i,i,j))
							 +
							 0.5*zeta*(x.DyV2C(Bx,i,j)+xold.DyV2C(Bx,i,j))
							 -
							 0.5*zeta*(x.DxV2C(By,i,j)+xold.DxV2C(By,i,j))
					  );


	}


#pragma omp parallel for
  for(int j = 0; j<x.ny; j++)
  for(int i = 0; i<x.nx; i++)
  {
	  // Store it in gamma, not being used for anything else...
	  res(Sxx_e,i,j) = 0;
	  res(Sxy_e,i,j) = 0;
	  res(Sxz_e,i,j) = 0;
	  res(Syy_e,i,j) = 0;
	  res(Syz_e,i,j) = 0;

	  // Store it in gamma, not being used for anything else...
	  res(Sxx_i,i,j) = 0;
	  res(Sxy_i,i,j) = 0;
	  res(Sxz_i,i,j) = 0;
	  res(Syy_i,i,j) = 0;
	  res(Syz_i,i,j) = 0;
  }


}
