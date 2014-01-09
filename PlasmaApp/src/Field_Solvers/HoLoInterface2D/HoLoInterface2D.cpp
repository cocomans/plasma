#include "HoLoInterface2D.h"
#include "../../RunData.h"
#include "../../OutputCtrl.h"
#include "../InitSolve.h"

HoLoInterface2D::HoLoInterface2D(PlasmaData* _pdata)
{

	pdata = _pdata;

#ifdef HAVE_MPI
//MPI_Init( &argc, &argv );
Epetra_MpiComm* Comm = new Epetra_MpiComm( MPI_COMM_SELF );
std::cout << *Comm << std::endl;
#else
Epetra_SerialComm* Comm = new Epetra_SerialComm();
std::cout << *Comm <<std::endl;
#endif
//////////////////////////////////////

// Initialize Simulation Parameters
printf("building SimParams\n");
simParams = Teuchos::rcp(new SimParams(Comm,pdata));

//	simParams->Print();

  comm      = Comm; // build comm here
  printf("Importing Params SimParams\n");
  ImportParams();
  printf("Building Problem\n");

  problem = Teuchos::rcp(new LowOrderProblem2D(simParams, comm));
  printf("Building JFNK\n");
  jfnk = Teuchos::rcp(new JFNK2D(problem, simParams, comm));



  X = problem->GetXW();
  Xold = problem->GetXoldW();
  ConsTerms = problem->GetConsTermsW();

  ResVec = Teuchos::rcp( new EpVecWrapper((problem->GetMapManager()),simParams) );


  HO_resid = 1.0;



}

HoLoInterface2D::~HoLoInterface2D()
{
}

void HoLoInterface2D::ImportParams()
{

  me_h   = simParams->GetParamsPhys()->get<double>("me_h");
  mi_h   = simParams->GetParamsPhys()->get<double>("mi_h");
  qe_h   = simParams->GetParamsPhys()->get<double>("qe_h");
  qi_h   = simParams->GetParamsPhys()->get<double>("qi_h");

  xi    = simParams->GetParamsPhys()->get<double>("xi");
  dx    = simParams->GetParamsPhys()->get<double>("dx");
  dt    = simParams->GetParamsPhys()->get<double>("dt");

  tol 	= 1e-6;

}


void HoLoInterface2D::init(
	PlasmaData* PData,
	FieldDataCPU* Fields_Old,
	NodeFieldData* Fields_Half,
	FieldDataCPU* Fields_Next,
	NodeHOMoments* Moments
)
{

   pdata        = PData;
   fields_old   = Fields_Old;
   fields_half  = Fields_Half;
   fields_next  = Fields_Next;
   moments_old  = Moments->moments_old;
   moments_next = Moments->moments_next;

  problem->Init();
  jfnk->Init();

  HoToLo(fields_old,moments_old, Xold);
  HoToLo(fields_old,moments_old, X);

  isConv = false;
}
//
//  We want to change the return argument: bool or int
//
void HoLoInterface2D::solve(PlasmaData* Pdata,
	     FieldDataCPU* Fields, //output
	     NodeHOMoments* Moments)
{


  // Solve!!!
  jfnk->Solve();


  // Return Fields to HO system
  LoToFields();
  // do some stuff to final state to put it back into fields next

  // Print (Optional)
//  jfnk->PrintStatus();
}


void HoLoInterface2D::InitStepSolve(PlasmaData* Pdata,
	     FieldDataCPU* Fields, //output
	     NodeHOMoments* Moments)
{

//	  if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1))
//	  {
//		  using namespace ES1D1V;
//
//		  LO_log.print("HoToLo\n");
//
//	   int Nelem = problem->GetMapManager()->LocSdNumElemX;
//	   for(int e = 0 ; e< Nelem ; e++)
//	     {
//	       int i = problem->GetMapManager()->ElemSdMap->GID(e);
//
//	       int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);
//	       (*Xold)[i_pt + n_e]   = moments_old->get_val(i, 0, 0, 0, HOMoments_charge);
//	//       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
//	       (*Xold)[i_pt + p_e]   = 2.0*moments_old->get_val(i, 0, 0, 0, HOMoments_currentx) - (*Xold)[i_pt + p_e];
//	       (*Xold)[i_pt + Sxx_e] = moments_old->get_val(i, 0, 0, 0, HOMoments_S2xx)/( (*Xold)[i_pt + n_e] );
//
//	       if(pdata->nspecies > 1){
//			   (*Xold)[i_pt + n_i]   = moments_old->get_val(i, 0, 0, 1, HOMoments_charge);
//			   (*Xold)[i_pt + p_i]   = 2.0*moments_old->get_val(i, 0, 0, 1, HOMoments_currentx) - (*Xold)[i_pt + p_i];
//			   (*Xold)[i_pt + Sxx_i] = moments_old->get_val(i, 0, 0, 1, HOMoments_S2xx)/((*Xold)[i_pt + n_i]);
//	       }
//	       else
//	       {
//			   (*Xold)[i_pt + n_i]   = 1.0;
//			   (*Xold)[i_pt + p_i]   = 0.0;
//			   (*Xold)[i_pt + Sxx_i] = 0.0;
//	       }
//
//	       (*Xold)[i_pt + Ex]     = fields_old->getE(i, 0, 0, 0);
//	     }
//	  }
//	jfnk->UpdateSolution(X);

	ConsTerms->PutScalar(0.0);
	HoToLo(fields_old,moments_old, Xold);

    // load new solution
    HoToLo(fields_next,moments_next,X);
    jfnk->UpdateSolution(X);

	// Solve!!!
	jfnk->Solve();

	// Return Fields to HO system
	LoToFields();
	// do some stuff to final state to put it back into fields next

	// Print (Optional)
	jfnk->PrintStatus();
}


realkind HoLoInterface2D::calc_residual(PlasmaData* pdata,
		  FieldDataCPU* Fields_next,
		  FieldDataCPU* Fields_old,
		  NodeHOMoments* moments)
{



	HoToLo(fields_old,moments_old, Xold);

	// load new solution
	HoToLo(fields_next,moments_next,X);
	jfnk->UpdateSolution(X);

	//    std::cout << *X;
	//    jfnk -> DumpSolution(0);

	// Check convergence
	//    LO_log << "Compute Residual Norm\n";
	//

	ConsTerms->PutScalar(0.0);
	jfnk->UpdateSolution(X);
	jfnk->ComputeResidual(ResVec);
	jfnk->UpdateSolution(X);
	jfnk->ComputeResidual(ConsTerms);
	jfnk->UpdateSolution(X);

	if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1))
	{
	  using namespace ES1D1V;
	  res = ResVec->subNorm2(1,Ex);


	}
	else if((pdata->ndimensions == 1)&&(pdata->nVelocity == 3))
	{
	  using namespace EM1D3V;
	  res = ResVec->subNorm2(5,Ex,Ey,Ez,Ay,Az);


	}
	else if((pdata->ndimensions == 2)&&(pdata->nVelocity == 3))
	{
	  using namespace EM2D3V;

	  res = ResVec->subNorm2(6,Ex,Ey,Ez,Bx,By,Bz);
	}




    //  Return
    return res;
};



void HoLoInterface2D::update_solution(void)
{
};

void HoLoInterface2D::LoToFields()
{
	//	Copy final solution to X
	jfnk->GetSolution(X);
	jfnk->UpdateSolution(X);
	jfnk->ComputeResidualNorm();

	if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1))
	{
	  using namespace ES1D1V;
		for(int j=0;j<pdata->ny;j++)
		for(int k=0;k<pdata->nz;k++)
		  for(int i=0;i<pdata->nx;i++)
		  {
			  fields_next->getE(i,j,k,0) = X->Get(Ex,i);
		  }


	}
	else if((pdata->ndimensions == 1)&&(pdata->nVelocity == 3))
	{
	  using namespace EM1D3V;
		for(int j=0;j<pdata->ny;j++)
			for(int k=0;k<pdata->nz;k++)
				for(int i=0;i<pdata->nx;i++)
				{
				  fields_next->getE(i,j,k,0) = X->Get(Ex,i);
				  fields_next->getE(i,j,k,1) = X->Get(Ey,i);
				  fields_next->getE(i,j,k,2) = X->Get(Ez,i);

				  fields_next->getA(i,j,k,1) = X->Get(Ay,i);
				  fields_next->getA(i,j,k,2) = X->Get(Az,i);

				  fields_next->getB(i,j,k,1) = -X->DxC2F(Az,i);
				  fields_next->getB(i,j,k,2) = X->DxC2F(Ay,i);

				}

	}
	else if((pdata->ndimensions == 2)&&(pdata->nVelocity == 3))
	{
		  using namespace EM2D3V;


		for(int k=0;k<pdata->nz;k++)
#pragma omp parallel for
			for(int j=0;j<pdata->ny;j++)
				for(int i=0;i<pdata->nx;i++)
				{
				  fields_next->getE(i,j,k,0) = X->Get(Ex,i,j) ;
				  fields_next->getE(i,j,k,1) = X->Get(Ey,i,j) ;
				  fields_next->getE(i,j,k,2) = X->Get(Ez,i,j) ;

				  fields_next->getA(i,j,k,0) = X->Get(Ax,i,j);
				  fields_next->getA(i,j,k,1) = X->Get(Ay,i,j);
				  fields_next->getA(i,j,k,2) = X->Get(Az,i,j);

				  fields_next->getB(i,j,k,0) = X->Get(Bx,i,j);
				  fields_next->getB(i,j,k,1) = X->Get(By,i,j);
				  fields_next->getB(i,j,k,2) = X->Get(Bz,i,j);

				  fields_next->getPhi(i,j,k) = 2.0*X->Get(Phi,i,j) - Xold->Get(Phi,i,j);

				  fields_next->getChi(i,j,k) = X->Get(Chi,i,j);

				}

	}



}
void HoLoInterface2D::LoToMoments()
{
	//	Copy final solution to X
	jfnk->GetSolution(X);

	if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1))
	{
	  using namespace ES1D1V;
		for(int j=0;j<pdata->ny;j++)
		for(int k=0;k<pdata->nz;k++)
		  for(int i=0;i<pdata->nx;i++)
		  {
			  fields_next->getE(i,j,k,0) = X->Get(Ex,i);
		  }


	}
	else if((pdata->ndimensions == 1)&&(pdata->nVelocity == 3))
	{
	  using namespace EM1D3V;
		for(int j=0;j<pdata->ny;j++)
			for(int k=0;k<pdata->nz;k++)
				for(int i=0;i<pdata->nx;i++)
				{
				  fields_next->getE(i,j,k,0) = X->Get(Ex,i);
				  fields_next->getE(i,j,k,1) = X->Get(Ey,i);
				  fields_next->getE(i,j,k,2) = X->Get(Ez,i);

				  fields_next->getA(i,j,k,1) = X->Get(Ay,i);
				  fields_next->getA(i,j,k,2) = X->Get(Az,i);

				  fields_next->getB(i,j,k,1) = -X->DxC2F(Az,i);
				  fields_next->getB(i,j,k,2) = X->DxC2F(Ay,i);

				}

	}
	else if((pdata->ndimensions == 2)&&(pdata->nVelocity == 3))
	{
		  using namespace EM2D3V;


		for(int k=0;k<pdata->nz;k++)
#pragma omp parallel for
			for(int j=0;j<pdata->ny;j++)
				for(int i=0;i<pdata->nx;i++)
				{
				  moments_next->get_val(i,j,k,0,HOMoments_charge) = X->Get(n_e,i,j);
				  moments_next->get_val(i,j,k,1,HOMoments_charge) = X->Get(n_i,i,j);

				  moments_next->get_val(i,j,k,0,HOMoments_currentx) = X->Get(px_e,i,j);
				  moments_next->get_val(i,j,k,0,HOMoments_currenty) = X->Get(py_e,i,j);
				  moments_next->get_val(i,j,k,0,HOMoments_currentz) = X->Get(pz_e,i,j);

				  moments_next->get_val(i,j,k,1,HOMoments_currentx) = X->Get(px_e,i,j);
				  moments_next->get_val(i,j,k,1,HOMoments_currenty) = X->Get(py_e,i,j);
				  moments_next->get_val(i,j,k,1,HOMoments_currentz) = X->Get(pz_e,i,j);





				}

	}


}


void HoLoInterface2D::HoToLo( FieldDataCPU* field_in, HOMomentsCPU* moments_in,
		Teuchos::RCP<EpVecWrapper>& v)
{
	  if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1))
	  {
		  using namespace ES1D1V;

		  LO_log.print("HoToLo\n");

	   for(int i = 0 ; i< pdata->nx ; i++)
	     {

		   v->Get(n_e,i)   = moments_in->get_val(i, 0, 0, 0, HOMoments_charge);
	//       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
		   v->Get(p_e,i)   = moments_in->get_val(i, 0, 0, 0, HOMoments_currentx);
		   v->Get(Sxx_e,i) = moments_in->get_val(i, 0, 0, 0, HOMoments_S2xx)/( v->Get(n_e,i));

	       if(pdata->nspecies > 1){
	    	   v->Get(n_i,i)   = moments_in->get_val(i, 0, 0, 1, HOMoments_charge);
	    	   v->Get(p_i,i)   = moments_in->get_val(i, 0, 0, 1, HOMoments_currentx);
	    	   v->Get(Sxx_i,i) = moments_in->get_val(i, 0, 0, 1, HOMoments_S2xx)/(v->Get(n_i,i));
	       }
	       else
	       {
	    	   v->Get(n_i,i)   = 1.0;
	    	   v->Get(p_i,i)   = 0.0;
	    	   v->Get(Sxx_i,i) = 0.0;
	       }

	       v->Get(Ex,i)     = field_in->getE(i, 0, 0, 0);
	     }
	  }
	  else if((pdata->ndimensions == 1)&&(pdata->nVelocity == 3))
	  {
		  using namespace EM1D3V;


		  LO_log.print("HoToLo\n");

	   for(int i = 0 ; i< pdata->nx ; i++)
	     {

	       v->Get(n_e,i)   = moments_in->get_val(i, 0, 0, 0, HOMoments_charge);
	//       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
	       v->Get(px_e,i)   = 0.5*(moments_in->get_val(i, 0, 0, 0, HOMoments_currentx)+moments_next->get_val(i, 0, 0, 0, HOMoments_currentx));
	       v->Get(py_e,i)   = 0.5*(moments_in->get_val(i, 0, 0, 0, HOMoments_currenty)+moments_next->get_val(i, 0, 0, 0, HOMoments_currenty));
	       v->Get(pz_e,i)   = 0.5*(moments_in->get_val(i, 0, 0, 0, HOMoments_currentz)+moments_next->get_val(i, 0, 0, 0, HOMoments_currentz));

	       v->Get(Sxx_e,i) = moments_in->get_val(i, 0, 0, 0, HOMoments_S2xx)/( v->Get(n_e,i) );
	       v->Get(Sxy_e,i) = moments_in->get_val(i, 0, 0, 0, HOMoments_S2xy)/( v->Get(n_e,i) );
	       v->Get(Sxz_e,i) = moments_in->get_val(i, 0, 0, 0, HOMoments_S2xz)/( v->Get(n_e,i) );

	       if(pdata->nspecies > 1){
	    	   v->Get(n_i,i)   = moments_in->get_val(i, 0, 0, 1, HOMoments_charge);
	    //       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
	    	   v->Get(px_i,i)   = 0.5*(moments_in->get_val(i, 0, 0, 1, HOMoments_currentx)+moments_next->get_val(i, 0, 0, 1, HOMoments_currentx));
	    	   v->Get(py_i,i)   = 0.5*(moments_in->get_val(i, 0, 0, 1, HOMoments_currenty)+moments_next->get_val(i, 0, 0, 1, HOMoments_currenty));
	    	   v->Get(pz_i,i)   = 0.5*(moments_in->get_val(i, 0, 0, 1, HOMoments_currentz)+moments_next->get_val(i, 0, 0, 1, HOMoments_currentz));

	    	   v->Get(Sxx_i,i) = moments_in->get_val(i, 0, 0, 1, HOMoments_S2xx)/( v->Get(n_i,i) );
	    	   v->Get(Sxy_i,i) = moments_in->get_val(i, 0, 0, 1, HOMoments_S2xy)/( v->Get(n_i,i) );
	    	   v->Get(Sxz_i,i) = moments_in->get_val(i, 0, 0, 1, HOMoments_S2xz)/( v->Get(n_i,i) );

	       }
	       else
	       {
	    	   v->Get(n_i,i)   = 1.0;
	    //       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
	    	   v->Get(px_i,i)  = 0.0;
	    	   v->Get(py_i,i)   = 0.0;
	    	   v->Get(pz_i,i)   = 0.0;

	    	   v->Get(Sxx_i,i) = 0.1;
	    	   v->Get(Sxy_i,i) = 0.1;
	    	   v->Get(Sxz_i,i) = 0.1;

	       }

	       v->Get(Ex,i)     = field_in->getE(i, 0, 0, 0);
	       v->Get(Ey,i)     = field_in->getE(i, 0, 0, 1);
	       v->Get(Ez,i)     = field_in->getE(i, 0, 0, 2);

	       v->Get(Bx,i)     = field_in->getB(i, 0, 0, 0);
	       v->Get(By,i)     = field_in->getB(i, 0, 0, 1);
	       v->Get(Bz,i)     = field_in->getB(i, 0, 0, 2);

	       v->Get(Ax,i)     = field_in->getA(i, 0, 0, 0);
	       v->Get(Ay,i)     = field_in->getA(i, 0, 0, 1);
	       v->Get(Az,i)     = field_in->getA(i, 0, 0, 2);

	     }
	  }
	  else if((pdata->ndimensions == 2)&&(pdata->nVelocity == 3))
	  {
		using namespace EM2D3V;


			LO_log.print("HoToLo\n");

#pragma omp parallel for
			for(int j=0;j<pdata->ny;j++)
			for(int i = 0 ; i< pdata->nx ; i++)
			{

				v->Get(n_e,i,j)   = moments_in->get_val(i, j, 0, 0, HOMoments_charge);
				//       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
				v->Get(px_e,i,j)   = (moments_in->get_val(i, j, 0, 0, HOMoments_currentx));
				v->Get(py_e,i,j)   = (moments_in->get_val(i, j, 0, 0, HOMoments_currenty));
				v->Get(pz_e,i,j)   = (moments_in->get_val(i, j, 0, 0, HOMoments_currentz));

				v->Get(Sxx_e,i,j) = moments_in->get_val(i, j, 0, 0, HOMoments_S2xx)/( v->Get(n_e,i,j) );
				v->Get(Sxy_e,i,j) = moments_in->get_val(i, j, 0, 0, HOMoments_S2xy)/( v->Get(n_e,i,j) );
				v->Get(Sxz_e,i,j) = moments_in->get_val(i, j, 0, 0, HOMoments_S2xz)/( v->Get(n_e,i,j) );
				v->Get(Syy_e,i,j) = moments_in->get_val(i, j, 0, 0, HOMoments_S2yy)/( v->Get(n_e,i,j) );
				v->Get(Syz_e,i,j) = moments_in->get_val(i, j, 0, 0, HOMoments_S2yz)/( v->Get(n_e,i,j) );
				v->Get(Szz_e,i,j) = moments_in->get_val(i, j, 0, 0, HOMoments_S2zz)/( v->Get(n_e,i,j) );

				if(pdata->nspecies > 1){
					v->Get(n_i,i,j)   = moments_in->get_val(i, j, 0, 1, HOMoments_charge);
					//       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
					v->Get(px_i,i,j)   = (moments_in->get_val(i, j, 0, 1, HOMoments_currentx));
					v->Get(py_i,i,j)   = (moments_in->get_val(i, j, 0, 1, HOMoments_currenty));
					v->Get(pz_i,i,j)   = (moments_in->get_val(i, j, 0, 1, HOMoments_currentz));

					v->Get(Sxx_i,i,j) = moments_in->get_val(i, j, 0, 1, HOMoments_S2xx)/( v->Get(n_i,i,j) );
					v->Get(Sxy_i,i,j) = moments_in->get_val(i, j, 0, 1, HOMoments_S2xy)/( v->Get(n_i,i,j) );
					v->Get(Sxz_i,i,j) = moments_in->get_val(i, j, 0, 1, HOMoments_S2xz)/( v->Get(n_i,i,j) );
					v->Get(Syy_i,i,j) = moments_in->get_val(i, j, 0, 1, HOMoments_S2yy)/( v->Get(n_i,i,j) );
					v->Get(Syz_i,i,j) = moments_in->get_val(i, j, 0, 1, HOMoments_S2yz)/( v->Get(n_i,i,j) );
					v->Get(Szz_i,i,j) = moments_in->get_val(i, j, 0, 1, HOMoments_S2zz)/( v->Get(n_i,i,j) );

				}
				else
				{
					v->Get(n_i,i,j)   = 1.0;

					v->Get(px_i,i,j)  = 0.0;
					v->Get(py_i,i,j)   = 0.0;
					v->Get(pz_i,i,j)   = 0.0;

					v->Get(Sxx_i,i,j) = 0.1;
					v->Get(Sxy_i,i,j) = 0.1;
					v->Get(Sxz_i,i,j) = 0.1;
					v->Get(Syy_i,i,j) = 0.1;
					v->Get(Syz_i,i,j) = 0.1;

				}

				v->Get(Ex,i,j)     = field_in->getE(i, j, 0, 0);
				v->Get(Ey,i,j)     = field_in->getE(i, j, 0, 1);
				v->Get(Ez,i,j)     = field_in->getE(i, j, 0, 2);

				v->Get(Bx,i,j)     = field_in->getB(i, j, 0, 0);
				v->Get(By,i,j)     = field_in->getB(i, j, 0, 1);
				v->Get(Bz,i,j)     = field_in->getB(i, j, 0, 2);

				v->Get(Ax,i,j)     = field_in->getA(i, j, 0, 0);
				v->Get(Ay,i,j)     = field_in->getA(i, j, 0, 1);
				v->Get(Az,i,j)     = field_in->getA(i, j, 0, 2);

				v->Get(Phi,i,j)     = field_in->getPhi(i, j, 0);
				v->Get(Chi,i,j)     = field_in->getChi(i, j, 0);

			}
	  }
}



// For standalone run
void HoLoInterface2D::SolveTimeStepping()
{
  


}
