#include "DiscretePDE_1D1Ves.h"

DiscretePDE_1D1Ves::DiscretePDE_1D1Ves(const Teuchos::RCP<MapManager> &map_manager,
	      const Teuchos::RCP<SimParams> &_params,
	      Epetra_Comm* _comm)
{
  comm       = _comm;
  simParams  = _params;
  mapManager = map_manager;

  stagMesh   = Teuchos::rcp(new StagMesh(mapManager, simParams, comm) );
  OvX        = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );
  OvXold     = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );

  ImportParams();
};

//DiscretePDE::~DiscretePDE(){}

////////////////////////////////////////////////////////////////////
// USER: Edit any constants you may need in SimParams.cpp  //
///////////////////////////////////////////////////////////////////
void DiscretePDE_1D1Ves::ImportParams()
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

void DiscretePDE_1D1Ves::EvaluateResidual( const Epetra_Vector& x,
				    const Epetra_Vector& xold,
				    const Epetra_Vector& gamma,
				    Epetra_Vector& res )
{
	using namespace ES1D1V;

  // Build an enum of variable names


//  std::cout << "gamma = " << std::endl;
//  std::cout << gamma;

  // Load x into a discrete operator wrapper
  OvX->SetVector(x);

  // Load xold into a discrte operator wrapper
  OvXold->SetVector(xold);

  //Global Operations (ie. averages)
  double pe_avg = 0.0;
  double pi_avg = 0.0;

  double pe_avg_old = 0.0;
  double pi_avg_old = 0.0;

  // double pe_avg = OvX->Mean(p_e);
  // double pi_avg = OvX->Mean(p_i);

  // double pe_avg_old = OvXold->Mean(p_e);
  // double pi_avg_old = OvXold->Mean(p_i);

  double jmean1 = 0;


  // Loop over elements
  int N = mapManager->LocSdNumElemX;
  for(int i = 0; i<N; i++)
  {
      // Convert elem index to first point index
      int j = mapManager->LocElem_to_LocPtSd(i,0);

      jmean1 += qe_h*OvX->Get(p_e, i) + qi_h*OvX->Get(p_i,i);

  }

  jmean1 /= (double)N;

//  std::cout << gamma ;

#pragma omp parallel for
  for(int i = 0; i<N; i++)
    {   

//	  printf("ne = %i, pe = %i, E = %i, pi = %i, ni = %i\n",n_e,p_e,E,p_i,n_i);
      // Convert elem index to first point index
      int j = mapManager->LocElem_to_LocPtSd(i,0);
            
      // Construct nonlinear terms
      double tmp_W = 0.25*(OvX->Get(Sxx_e, i-1) + OvXold->Get(Sxx_e, i-1)) * (OvX->Get(n_e,i-1)+OvXold->Get(n_e,i-1));
      double tmp_C = 0.25*(OvX->Get(Sxx_e, i) + OvXold->Get(Sxx_e, i)) * (OvX->Get(n_e,i)+OvXold->Get(n_e,i));
      double tmp_E = 0.25*(OvX->Get(Sxx_e, i+1) + OvXold->Get(Sxx_e, i+1)) * (OvX->Get(n_e,i+1)+OvXold->Get(n_e,i+1));
      double Dx_Sxx_ne = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);


      tmp_W = 0.25*(OvX->Get(Sxx_i, i-1) + OvXold->Get(Sxx_i, i-1)) * (OvX->Get(n_i,i-1)+OvXold->Get(n_i,i-1));
      tmp_C = 0.25*(OvX->Get(Sxx_i, i) + OvXold->Get(Sxx_i, i)) * (OvX->Get(n_i,i)+OvXold->Get(n_i,i));
      tmp_E = 0.25*(OvX->Get(Sxx_i, i+1) + OvXold->Get(Sxx_i, i+1)) * (OvX->Get(n_i,i+1)+OvXold->Get(n_i,i+1));
      double Dx_Sxx_ni = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);


      double p_e1 = (OvX->Get(p_e, i) - OvXold->Get(p_e, i))+0.5*(OvX->Get(p_e, i) + OvXold->Get(p_e, i));
      double p_i1 = (OvX->Get(p_i, i) - OvXold->Get(p_i, i))+0.5*(OvX->Get(p_i, i) + OvXold->Get(p_i, i));


      // Fill each unknown
      res[j+n_e] = (OvX->Get(n_e, i) - OvXold->Get(n_e, i)) + dt*( OvX->DxF2C(p_e, i, dx))  - gamma[j+n_e];

      res[j+p_e] = me_h*(OvX->Get(p_e, i) - OvXold->Get(p_e, i)) + dt*(me_h*Dx_Sxx_ne
    		  - 0.25*(qe_h*OvX->GetC2F(n_e,i)+qe_h*OvXold->GetC2F(n_e, i)) * (OvX->Get(Ex,i)+OvXold->Get(Ex,i))
    		  ) - gamma[j+p_e];
								   

      res[j+Ex] = xi*xi*(OvX->Get(Ex, i) - OvXold->Get(Ex, i))
    				  + dt*(qe_h*OvX->Get(p_e, i) + qi_h*OvX->Get(p_i,i) - (jmean1));

      res[j+p_i] = mi_h*(OvX->Get(p_i, i) - OvXold->Get(p_i, i)) + dt*(mi_h*Dx_Sxx_ni
    		  - 0.25*(qi_h*OvX->GetC2F(n_i,i)+qi_h*OvXold->GetC2F(n_i, i)) * (OvX->Get(Ex,i)+OvXold->Get(Ex,i))
    		 	 ) - gamma[j+p_i];


      res[j+n_i] = (OvX->Get(n_i, i) - OvXold->Get(n_i, i)) + dt*( OvX->DxF2C(p_i, i, dx) )  - gamma[j+n_i];

      res[j+Sxx_e] = 0.0;
      res[j+Sxx_i] = 0.0;
//      res[j+p_i] = 0.0;
//      res[j+n_i] = 0.0;
//      res[j+p_e] = 0.0;
//      res[j+n_e] = 0.0;

    }


}

