#include "DiscretePDE2D_1D1Ves.h"

DiscretePDE2D_1D1Ves::DiscretePDE2D_1D1Ves(const Teuchos::RCP<SimParams> &_params,
	      Epetra_Comm* _comm)
{
  comm       = _comm;
  simParams  = _params;

  ImportParams();
};

//DiscretePDE::~DiscretePDE(){}

////////////////////////////////////////////////////////////////////
// USER: Edit any constants you may need in SimParams.cpp  //
///////////////////////////////////////////////////////////////////
void DiscretePDE2D_1D1Ves::ImportParams()
{

  me_h   = simParams->GetParamsPhys()->get<double>("me_h");
  mi_h   = simParams->GetParamsPhys()->get<double>("mi_h");
  qe_h   = simParams->GetParamsPhys()->get<double>("qe_h");
  qi_h   = simParams->GetParamsPhys()->get<double>("qi_h");

  xi    = simParams->GetParamsPhys()->get<double>("xi");
  dx    = simParams->GetParamsPhys()->get<double>("dx");
  dt    = simParams->GetParamsPhys()->get<double>("dt");
}

////////////////////////////////////////////////////////////////////
// USER: This is where users must build the main residual calc.   //
// Propogate the enum below with names for each of your unknowns. //
// The number of these unknowns must match "Number Unknowns"      //
///////////////////////////////////////////////////////////////////

void DiscretePDE2D_1D1Ves::EvaluateResidual(EpVecWrapper& x,
		  const EpVecWrapper& Xold,
		  const EpVecWrapper& ConsTerm,
		  EpVecWrapper& res)
{
	using namespace ES1D1V;


  double jmean1 = qe_h*x.mean(p_e) + qi_h*x.mean(p_i);




//  std::cout << gamma ;

#pragma omp parallel for
  for(int i = 0; i<x.nx; i++)
    {   

//	  printf("ne = %i, pe = %i, E = %i, pi = %i, ni = %i\n",n_e,p_e,E,p_i,n_i);
      // Convert elem index to first point index

            
      // Construct nonlinear terms
      double tmp_W = 0.25*(x.Get(Sxx_e, i-1) + Xold.Get(Sxx_e, i-1)) * (x.Get(n_e,i-1)+Xold.Get(n_e,i-1));
      double tmp_C = 0.25*(x.Get(Sxx_e, i) + Xold.Get(Sxx_e, i)) * (x.Get(n_e,i)+Xold.Get(n_e,i));
      double tmp_E = 0.25*(x.Get(Sxx_e, i+1) + Xold.Get(Sxx_e, i+1)) * (x.Get(n_e,i+1)+Xold.Get(n_e,i+1));
      double Dx_Sxx_ne = x.DxC2F(tmp_W, tmp_C, tmp_E);


      tmp_W = 0.25*(x.Get(Sxx_i, i-1) + Xold.Get(Sxx_i, i-1)) * (x.Get(n_i,i-1)+Xold.Get(n_i,i-1));
      tmp_C = 0.25*(x.Get(Sxx_i, i) + Xold.Get(Sxx_i, i)) * (x.Get(n_i,i)+Xold.Get(n_i,i));
      tmp_E = 0.25*(x.Get(Sxx_i, i+1) + Xold.Get(Sxx_i, i+1)) * (x.Get(n_i,i+1)+Xold.Get(n_i,i+1));
      double Dx_Sxx_ni = x.DxC2F(tmp_W, tmp_C, tmp_E);


      double p_e1 = (x.Get(p_e, i) - Xold.Get(p_e, i))+0.5*(x.Get(p_e, i) + Xold.Get(p_e, i));
      double p_i1 = (x.Get(p_i, i) - Xold.Get(p_i, i))+0.5*(x.Get(p_i, i) + Xold.Get(p_i, i));


      // Fill each unknown
      res(n_e,i) = (x.Get(n_e, i) - Xold.Get(n_e, i)) + dt*( x.DxF2C(p_e, i))  - ConsTerm(n_e,i);

      res(p_e,i) = me_h*(x.Get(p_e, i) - Xold.Get(p_e, i)) + dt*(me_h*Dx_Sxx_ne
    		  - 0.25*(qe_h*x.GetC2F(n_e,i)+qe_h*Xold.GetC2F(n_e, i)) * (x.Get(Ex,i)+Xold.Get(Ex,i))
    		  ) - ConsTerm(p_e,i);
								   

      res(Ex,i) = xi*xi*(x.Get(Ex, i) - Xold.Get(Ex, i))
    				  + dt*(qe_h*x.Get(p_e, i) + qi_h*x.Get(p_i,i) - (jmean1));

      res(p_i,i) = mi_h*(x.Get(p_i, i) - Xold.Get(p_i, i)) + dt*(mi_h*Dx_Sxx_ni
    		  - 0.25*(qi_h*x.GetC2F(n_i,i)+qi_h*Xold.GetC2F(n_i, i)) * (x.Get(Ex,i)+Xold.Get(Ex,i))
    		 	 ) - ConsTerm(p_i,i);


      res(n_i,i) = (x.Get(n_i, i) - Xold.Get(n_i, i)) + dt*( x.DxF2C(p_i, i) )  - ConsTerm(n_i,i);


//      res(n_e,i) = 0.0;
//      res(n_i,i) = 0.0;
//      res(p_e,i) = 0.0;
//      res(p_i,i) = 0.0;
//      res[j+Sxx_e] = 0.0;
//      res[j+Sxx_i] = 0.0;
//      res[j+p_i] = 0.0;
//      res[j+n_i] = 0.0;
//      res[j+p_e] = 0.0;
//      res[j+n_e] = 0.0;

    }


}

