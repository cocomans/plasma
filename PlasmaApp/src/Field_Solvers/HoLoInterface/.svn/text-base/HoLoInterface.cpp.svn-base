#include "HoLoInterface.h"
#include "../../RunData.h"
#include "../../OutputCtrl.h"
#include "../InitSolve.h"

HoLoInterface::HoLoInterface(PlasmaData* _pdata)
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

  problem = Teuchos::rcp(new LowOrderProblem(simParams, comm));
  printf("Building JFNK\n");
  jfnk = Teuchos::rcp(new JFNK(problem, simParams, comm));
  printf("Building DiffMath\n");
  OvX = Teuchos::rcp(new DiffMath(problem->GetMapManager(),simParams,comm));


  X = problem->GetX();
  Xold = problem->GetXold();
  ConsTerms = problem->GetConsTerms();

  HO_resid = 1.0;


}

HoLoInterface::~HoLoInterface()
{
}

void HoLoInterface::ImportParams()
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


void HoLoInterface::init(
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
void HoLoInterface::solve(PlasmaData* Pdata,
	     FieldDataCPU* Fields, //output
	     NodeHOMoments* Moments)
{


//  ConsTerms->PutScalar(0.0);
//  jfnk->UpdateSolution(X);
//  jfnk->ComputeResidual(ConsTerms);
//  jfnk->UpdateSolution(X);

//  res = jfnk->ComputeResidualNorm();


  printf("res = %e should be 0\n",res);
//  std::cout << *ConsTerms;

//	InitSolve(pdata,fields_next,moments_next);
//
//	double jmeanx1 = 0.0;
//	for(int i=0;i<pdata->nx;i++)
//	{
//	  jmeanx1 += qe_h*moments_next->get_val(i, 0, 0, 0, HOMoments_currentx) + qi_h*moments_next->get_val(i, 0, 0, 1, HOMoments_currentx);
//	}
//
//	jmeanx1 /= (double)(pdata->nx);
//
//	for(int i=0;i<pdata->nx;i++)
//	{
//		fields_next->getE(i,0,0,0) = fields_old->getE(i,0,0,0) - (qe_h*moments_next->get_val(i, 0, 0, 0, HOMoments_currentx) + qi_h*moments_next->get_val(i, 0, 0, 1, HOMoments_currentx) - jmeanx1)*dt/(xi*xi);
//		fields_next->getE(i,0,0,1) = -2.0*(fields_next->getA(i,0,0,1) - fields_old->getA(i,0,0,1))/dt - fields_old->getE(i,0,0,1);
//		fields_next->getE(i,0,0,2) = -2.0*(fields_next->getA(i,0,0,2) - fields_old->getA(i,0,0,2))/dt - fields_old->getE(i,0,0,2);
//
//	}
//
//	HoToLo(fields_next,moments_next,X);
//	jfnk->UpdateSolution(X);


//  throw(1);
  // Solve!!!
  jfnk->Solve();

  // Return Fields to HO system
  LoToFields();
  // do some stuff to final state to put it back into fields next

  // Print (Optional)
//  jfnk->PrintStatus();
}


void HoLoInterface::InitStepSolve(PlasmaData* Pdata,
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


realkind HoLoInterface::calc_residual(PlasmaData* pdata,
		  FieldDataCPU* Fields_next,
		  FieldDataCPU* Fields_old,
		  NodeHOMoments* moments)
{

	  if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1))
	  {
		  	using namespace ES1D1V;
			int Nelem = problem->GetMapManager()->LocSdNumElemX;
			double density_diff[2] = {0,0};
			double current_diff[2] = {0,0};

			OvX->SetVector(*X);
			for(int e = 0 ; e< Nelem ; e++)
			 {
			//			printf("Dumping element %i\n",e);
			   int i = problem->GetMapManager()->ElemSdMap->GID(e);
			   int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);

//			   LO_log.print("E-Field errors[%i]: %e, %e\n",i,Ey_err,Ez_err);

			   density_diff[0] += pow((moments_next->get_val(i,0,0,0,HOMoments_charge) - (*X)[i_pt + n_e])
			                           /(fabs(moments_next->get_val(i,0,0,0,HOMoments_charge)) + fabs((*X)[i_pt + n_e])),2.0);
			   density_diff[1] += pow((moments_next->get_val(i,0,0,1,HOMoments_charge) - (*X)[i_pt + n_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_charge)) + fabs((*X)[i_pt + n_i])),2.0);

			   current_diff[0] += pow((moments_next->get_val(i,0,0,0,HOMoments_currentx) - (*X)[i_pt + p_e])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,0,HOMoments_currentx)) + fabs((*X)[i_pt + p_e])),2.0);

			   current_diff[1] += pow((moments_next->get_val(i,0,0,1,HOMoments_currentx) - (*X)[i_pt + p_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_currentx)) + fabs((*X)[i_pt + p_i])),2.0);


			 }

			density_diff[0] = sqrt(density_diff[0])/((double)Nelem);
			density_diff[1] = sqrt(density_diff[1])/((double)Nelem);

			current_diff[0] = sqrt(current_diff[0])/((double)Nelem);
			current_diff[1] = sqrt(current_diff[1])/((double)Nelem);


			res = pow(density_diff[0],2) + pow(density_diff[1],2) + pow(current_diff[0],2) + pow(current_diff[1],2);

			res = sqrt(res)/4.0;



	  }
	  else
	  {

		  	using namespace EM1D3V;
			int Nelem = problem->GetMapManager()->LocSdNumElemX;
			double density_diff[2] = {0,0};
			double current_diff[6] = {0,0,0,0,0,0};

			OvX->SetVector(*X);
			for(int e = 0 ; e< Nelem ; e++)
			 {
			//			printf("Dumping element %i\n",e);
			   int i = problem->GetMapManager()->ElemSdMap->GID(e);
			   int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);

//			   LO_log.print("E-Field errors[%i]: %e, %e\n",i,Ey_err,Ez_err);

			   density_diff[0] += pow((moments_next->get_val(i,0,0,0,HOMoments_charge) - (*X)[i_pt + n_e])
			                           /(fabs(moments_next->get_val(i,0,0,0,HOMoments_charge)) + fabs((*X)[i_pt + n_e])),2.0);
			   density_diff[1] += pow((moments_next->get_val(i,0,0,1,HOMoments_charge) - (*X)[i_pt + n_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_charge)) + fabs((*X)[i_pt + n_i])),2.0);

			   current_diff[0] += pow((moments_next->get_val(i,0,0,0,HOMoments_currentx) - (*X)[i_pt + px_e])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,0,HOMoments_currentx)) + fabs((*X)[i_pt + px_e])),2.0);
			   current_diff[1] += pow((moments_next->get_val(i,0,0,0,HOMoments_currenty) - (*X)[i_pt + py_e])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,0,HOMoments_currenty)) + fabs((*X)[i_pt + py_e])),2.0);
			   current_diff[2] += pow((moments_next->get_val(i,0,0,0,HOMoments_currentz) - (*X)[i_pt + pz_e])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,0,HOMoments_currentz)) + fabs((*X)[i_pt + pz_e])),2.0);

			   current_diff[3] += pow((moments_next->get_val(i,0,0,1,HOMoments_currentx) - (*X)[i_pt + px_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_currentx)) + fabs((*X)[i_pt + px_i])),2.0);
			   current_diff[4] += pow((moments_next->get_val(i,0,0,1,HOMoments_currenty) - (*X)[i_pt + py_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_currenty)) + fabs((*X)[i_pt + py_i])),2.0);
			   current_diff[5] += pow((moments_next->get_val(i,0,0,1,HOMoments_currentz) - (*X)[i_pt + pz_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_currentz)) + fabs((*X)[i_pt + pz_i])),2.0);



			 }

			density_diff[0] = sqrt(density_diff[0]/((double)Nelem));
			density_diff[1] = sqrt(density_diff[1]/((double)Nelem));

			current_diff[0] = sqrt(density_diff[0]/((double)Nelem));
			current_diff[1] = sqrt(density_diff[1]/((double)Nelem));
			current_diff[2] = sqrt(density_diff[2]/((double)Nelem));

			current_diff[3] = sqrt(density_diff[3]/((double)Nelem));
			current_diff[4] = sqrt(density_diff[4]/((double)Nelem));
			current_diff[5] = sqrt(density_diff[5]/((double)Nelem));

			res = density_diff[0];

			res = fmax(res,density_diff[0]);


			for (int i=0;i<6;i++)
				res = fmax(res,current_diff[i]);


	  }

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
    res = fmin(jfnk->ComputeResidualNorm(),res);
    jfnk->UpdateSolution(X);
    jfnk->ComputeResidual(ConsTerms);
    jfnk->UpdateSolution(X);

	int Nelem = problem->GetMapManager()->LocSdNumElemX;
	double E_res = 0;

	for(int e = 0 ; e< Nelem ; e++)
	 {
	  	using namespace ES1D1V;

		   int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);

		   E_res += pow((*ConsTerms)[i_pt + Ex]/pdata->dt,2);
	 }

	E_res = sqrt(E_res)/((double)Nelem);
	res = fmin(res,E_res);


	HO_resid = jfnk->ComputeResidualNorm();

    res = fmin(HO_resid,res);
//    jfnk->UpdateSolution(X,fmax(fmin(res/100,1.0e-9),5.0e-11));






    //  Return
    return res;
};



void HoLoInterface::update_solution(void)
{
};

void HoLoInterface::LoToFields()
{
	//	Copy final solution to X
	jfnk->GetSolution(X);

//	LO_log.print("Finished getting solution\n");

	//	Copy E from X into FieldsNext

	  int ipiccard_outer = pdata->rstatus->ipiccard_outer;
	  int istep = pdata->rstatus->istep;

	  if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1))
	  {
			using namespace ES1D1V;
			LO_log.print("Loading ES fields\n");
			int Nelem = problem->GetMapManager()->LocSdNumElemX;

			for(int e = 0 ; e< Nelem ; e++)
			{
			//			printf("Dumping element %i\n",e);

			int i = problem->GetMapManager()->ElemSdMap->GID(e);

			int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);
			for(int j=0;j<pdata->ny;j++)
			for(int k=0;k<pdata->nz;k++)
			   fields_next->getE(i,j,k,0) 	= (*X)[i_pt + Ex];
			 }

//			int Nelem = problem->GetMapManager()->LocSdNumElemX;
			double density_diff[2] = {0,0};
			double current_diff[6] = {0,0,0,0,0,0};

			OvX->SetVector(*X);
			for(int e = 0 ; e< Nelem ; e++)
			 {
			//			printf("Dumping element %i\n",e);
			   int i = problem->GetMapManager()->ElemSdMap->GID(e);
			   int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);

//			   LO_log.print("E-Field errors[%i]: %e, %e\n",i,Ey_err,Ez_err);

			   density_diff[0] += pow((moments_next->get_val(i,0,0,0,HOMoments_charge) - (*X)[i_pt + n_e])
			                           /(fabs(moments_next->get_val(i,0,0,0,HOMoments_charge)) + fabs((*X)[i_pt + n_e])),2.0);
			   density_diff[1] += pow((moments_next->get_val(i,0,0,1,HOMoments_charge) - (*X)[i_pt + n_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_charge)) + fabs((*X)[i_pt + n_i])),2.0);

			   current_diff[0] += pow((moments_next->get_val(i,0,0,0,HOMoments_currentx) - (*X)[i_pt + p_e])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,0,HOMoments_currentx)) + fabs((*X)[i_pt + p_e])),2.0);

			   current_diff[1] += pow((moments_next->get_val(i,0,0,1,HOMoments_currentx) - (*X)[i_pt + p_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_currentx)) + fabs((*X)[i_pt + p_i])),2.0);


			 }

			double E_res = 0;

			for(int e = 0 ; e< Nelem ; e++)
			 {
			  	using namespace ES1D1V;

				   int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);

				   E_res += pow((*ConsTerms)[i_pt + Ex],2);
			 }

			E_res = sqrt(E_res)/((double)Nelem);
			density_diff[0] = sqrt(density_diff[0])/((double)Nelem);
			density_diff[1] = sqrt(density_diff[1])/((double)Nelem);

			current_diff[0] = sqrt(current_diff[0])/((double)Nelem);
			current_diff[1] = sqrt(current_diff[1])/((double)Nelem);
			current_diff[3] = res;
			current_diff[4] = HO_resid;
			current_diff[5] = E_res;

			LO_log.print("LO and HO moment differences [%i, %i]\n",istep,ipiccard_outer);
			LO_log.print("Density: %e %e\n",density_diff[0],density_diff[1]);
			LO_log.print("Current 0: %e %e %e\n",current_diff[0],current_diff[1],current_diff[2]);
			LO_log.print("Current 1: %e %e %e\n",current_diff[3],current_diff[4],current_diff[5]);

			pdata->output_ctrl->RecordMomentDiffs(density_diff,current_diff);


	  }
	  else
	  {
//		  LO_log.print("Loading EM fields\n");



		  /* @TODO This really should return B, which means we
		   * need to do a Curl(A), but we are lazy and want to get
		   * this somewhat working.
		   *
		   */
		  	using namespace EM1D3V;
			int Nelem = problem->GetMapManager()->LocSdNumElemX;
			double density_diff[2] = {0,0};
			double current_diff[6] = {0,0,0,0,0,0};

			OvX->SetVector(*X);

			  double jmeanx1 = 0;
			  double jmeany1 = 0;
			  double jmeanz1 = 0;

			  // Loop over elements
			  for(int i = 0; i<Nelem; i++)
			  {
			      // Convert elem index to first point index

			      jmeanx1 += qe_h*OvX->GetF2C(px_e, i) + qi_h*OvX->GetF2C(px_i,i);
			      jmeany1 += qe_h*OvX->Get(py_e, i) + qi_h*OvX->Get(py_i,i);
			      jmeanz1 += qe_h*OvX->Get(pz_e, i) + qi_h*OvX->Get(pz_i,i);


			  }

			  jmeanx1 /= (double)(Nelem);
			  jmeany1 /= (double)(Nelem);
			  jmeanz1 /= (double)(Nelem);
			for(int e = 0 ; e< Nelem ; e++)
			 {
			//			printf("Dumping element %i\n",e);

			   int i = problem->GetMapManager()->ElemSdMap->GID(e);

			   int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);
			   for(int j=0;j<pdata->ny;j++)
			   for(int k=0;k<pdata->nz;k++)
			   {
				   fields_next->getE(i,j,k,0) 	= (*X)[i_pt + Ex];//fields_old->getE(i,0,0,0) - dt*(qe_h*OvX->Get(px_e, i) + qi_h*OvX->Get(px_i,i) - jmeanx1)/(xi*xi);
				   fields_next->getE(i,j,k,1) 	= -2.0*(OvX->Get(Ay,i) - (*Xold)[i_pt + Ay])/dt - fields_old->getE(i,0,0,1);
				   fields_next->getE(i,j,k,2) 	= -2.0*(OvX->Get(Az,i) - (*Xold)[i_pt + Az])/dt - fields_old->getE(i,0,0,2);

				   fields_next->getB(i,j,k,1) 	= -OvX->DxC2F(Az,i,dx);
				   fields_next->getB(i,j,k,2) 	= OvX->DxC2F(Ay,i,dx);

				   fields_next->getA(i,j,k,0) 	= (*X)[i_pt + Ax];
				   fields_next->getA(i,j,k,1) 	= (*X)[i_pt + Ay];
				   fields_next->getA(i,j,k,2) 	= (*X)[i_pt + Az];


			   }

			   double Ey_err = (fields_next->getE(i,0,0,1) - (*X)[i_pt + Ey])/fabs(fields_next->getE(i,0,0,1));
			   double Ez_err = (fields_next->getE(i,0,0,2) - (*X)[i_pt + Ez])/fabs(fields_next->getE(i,0,0,2));

//			   LO_log.print("E-Field errors[%i]: %e, %e\n",i,Ey_err,Ez_err);

			   density_diff[0] += pow((moments_next->get_val(i,0,0,0,HOMoments_charge) - (*X)[i_pt + n_e])
			                           /(fabs(moments_next->get_val(i,0,0,0,HOMoments_charge)) + fabs((*X)[i_pt + n_e])),2.0);
			   density_diff[1] += pow((moments_next->get_val(i,0,0,1,HOMoments_charge) - (*X)[i_pt + n_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_charge)) + fabs((*X)[i_pt + n_i])),2.0);

			   current_diff[0] += pow((moments_next->get_val(i,0,0,0,HOMoments_currentx) - (*X)[i_pt + px_e])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,0,HOMoments_currentx)) + fabs((*X)[i_pt + px_e])),2.0);
			   current_diff[1] += pow((moments_next->get_val(i,0,0,0,HOMoments_currenty) - (*X)[i_pt + py_e])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,0,HOMoments_currenty)) + fabs((*X)[i_pt + py_e])),2.0);
			   current_diff[2] += pow((moments_next->get_val(i,0,0,0,HOMoments_currentz) - (*X)[i_pt + pz_e])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,0,HOMoments_currentz)) + fabs((*X)[i_pt + pz_e])),2.0);

			   current_diff[3] += pow((moments_next->get_val(i,0,0,1,HOMoments_currentx) - (*X)[i_pt + px_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_currentx)) + fabs((*X)[i_pt + px_i])),2.0);
			   current_diff[4] += pow((moments_next->get_val(i,0,0,1,HOMoments_currenty) - (*X)[i_pt + py_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_currenty)) + fabs((*X)[i_pt + py_i])),2.0);
			   current_diff[5] += pow((moments_next->get_val(i,0,0,1,HOMoments_currentz) - (*X)[i_pt + pz_i])
					   	   	   	   	   /(fabs(moments_next->get_val(i,0,0,1,HOMoments_currentz)) + fabs((*X)[i_pt + pz_i])),2.0);



			 }

			density_diff[0] = sqrt(density_diff[0]/((double)Nelem));
			density_diff[1] = sqrt(density_diff[1]/((double)Nelem));

			current_diff[0] = sqrt(density_diff[0]/((double)Nelem));
			current_diff[1] = sqrt(density_diff[1]/((double)Nelem));
			current_diff[2] = sqrt(density_diff[2]/((double)Nelem));

			current_diff[3] = sqrt(density_diff[3]/((double)Nelem));
			current_diff[4] = sqrt(density_diff[4]/((double)Nelem));
			current_diff[5] = sqrt(density_diff[5]/((double)Nelem));

//			res = density_diff[0];
//
//			res = fmax(HO_resid,density_diff[0]);
//
//
//			for (int i=0;i<6;i++)
//				res = fmax(res,current_diff[i]);


			LO_log.print("LO and HO moment differences [%i, %i]\n",istep,ipiccard_outer);
			LO_log.print("Density: %e %e\n",density_diff[0],density_diff[1]);
			LO_log.print("Current 0: %e %e %e\n",current_diff[0],current_diff[1],current_diff[2]);
			LO_log.print("Current 1: %e %e %e\n",current_diff[3],current_diff[4],current_diff[5]);


			pdata->output_ctrl->RecordMomentDiffs(density_diff,current_diff);

	  }
	//	FieldsNext->get(X
}


void HoLoInterface::HoToLo( FieldDataCPU* field_in, HOMomentsCPU* moments_in, const Teuchos::RCP<Epetra_Vector>& v)
{

  if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1))
  {
	  using namespace ES1D1V;

	  LO_log.print("HoToLo\n");
  
   int Nelem = problem->GetMapManager()->LocSdNumElemX;
   for(int e = 0 ; e< Nelem ; e++)
     {
       int i = problem->GetMapManager()->ElemSdMap->GID(e);

       int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);
       (*v)[i_pt + n_e]   = moments_in->get_val(i, 0, 0, 0, HOMoments_charge);
//       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
       (*v)[i_pt + p_e]   = moments_in->get_val(i, 0, 0, 0, HOMoments_currentx);
       (*v)[i_pt + Sxx_e] = moments_in->get_val(i, 0, 0, 0, HOMoments_S2xx)/( (*v)[i_pt + n_e] );

       if(pdata->nspecies > 1){
		   (*v)[i_pt + n_i]   = moments_in->get_val(i, 0, 0, 1, HOMoments_charge);
		   (*v)[i_pt + p_i]   = moments_in->get_val(i, 0, 0, 1, HOMoments_currentx);
		   (*v)[i_pt + Sxx_i] = moments_in->get_val(i, 0, 0, 1, HOMoments_S2xx)/((*v)[i_pt + n_i]);
       }
       else
       {
		   (*v)[i_pt + n_i]   = 1.0;
		   (*v)[i_pt + p_i]   = 0.0;
		   (*v)[i_pt + Sxx_i] = 0.0;
       }

       (*v)[i_pt + Ex]     = field_in->getE(i, 0, 0, 0);
     }
  }
  else
  {
	  using namespace EM1D3V;


	  LO_log.print("HoToLo\n");

   int Nelem = problem->GetMapManager()->LocSdNumElemX;
   for(int e = 0 ; e< Nelem ; e++)
     {
       int i = problem->GetMapManager()->ElemSdMap->GID(e);
      
       int i_pt = problem->GetMapManager()->LocElem_to_LocPtSd(e, 0);
       (*v)[i_pt + n_e]   = moments_in->get_val(i, 0, 0, 0, HOMoments_charge);
//       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
       (*v)[i_pt + px_e]   = 0.5*(moments_in->get_val(i, 0, 0, 0, HOMoments_currentx)+moments_next->get_val(i, 0, 0, 0, HOMoments_currentx));
       (*v)[i_pt + py_e]   = 0.5*(moments_in->get_val(i, 0, 0, 0, HOMoments_currenty)+moments_next->get_val(i, 0, 0, 0, HOMoments_currenty));
       (*v)[i_pt + pz_e]   = 0.5*(moments_in->get_val(i, 0, 0, 0, HOMoments_currentz)+moments_next->get_val(i, 0, 0, 0, HOMoments_currentz));

       (*v)[i_pt + Sxx_e] = moments_in->get_val(i, 0, 0, 0, HOMoments_S2xx)/( (*v)[i_pt + n_e] );
       (*v)[i_pt + Sxy_e] = moments_in->get_val(i, 0, 0, 0, HOMoments_S2xy)/( (*v)[i_pt + n_e] );
       (*v)[i_pt + Sxz_e] = moments_in->get_val(i, 0, 0, 0, HOMoments_S2xz)/( (*v)[i_pt + n_e] );

       if(pdata->nspecies > 1){
           (*v)[i_pt + n_i]   = moments_in->get_val(i, 0, 0, 1, HOMoments_charge);
    //       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
           (*v)[i_pt + px_i]   = 0.5*(moments_in->get_val(i, 0, 0, 1, HOMoments_currentx)+moments_next->get_val(i, 0, 0, 1, HOMoments_currentx));
           (*v)[i_pt + py_i]   = 0.5*(moments_in->get_val(i, 0, 0, 1, HOMoments_currenty)+moments_next->get_val(i, 0, 0, 1, HOMoments_currenty));
           (*v)[i_pt + pz_i]   = 0.5*(moments_in->get_val(i, 0, 0, 1, HOMoments_currentz)+moments_next->get_val(i, 0, 0, 1, HOMoments_currentz));

           (*v)[i_pt + Sxx_i] = moments_in->get_val(i, 0, 0, 1, HOMoments_S2xx)/( (*v)[i_pt + n_i] );
           (*v)[i_pt + Sxy_i] = moments_in->get_val(i, 0, 0, 1, HOMoments_S2xy)/( (*v)[i_pt + n_i] );
           (*v)[i_pt + Sxz_i] = moments_in->get_val(i, 0, 0, 1, HOMoments_S2xz)/( (*v)[i_pt + n_i] );

       }
       else
       {
           (*v)[i_pt + n_i]   = 1.0;
    //       printf("my moments[%i] = %e\n", i,1.0-moments_in->get_val(i, 0, 0, 0, HOMoments_charge));
           (*v)[i_pt + px_i]   = 0.0;
           (*v)[i_pt + py_i]   = 0.0;
           (*v)[i_pt + pz_i]   = 0.0;

           (*v)[i_pt + Sxx_i] = 0.1;
           (*v)[i_pt + Sxy_i] = 0.1;
           (*v)[i_pt + Sxz_i] = 0.1;

       }

       (*v)[i_pt + Ex]     = field_in->getE(i, 0, 0, 0);
       (*v)[i_pt + Ey]     = field_in->getE(i, 0, 0, 1);
       (*v)[i_pt + Ez]     = field_in->getE(i, 0, 0, 2);

       (*v)[i_pt + Bx]     = field_in->getB(i, 0, 0, 0);
       (*v)[i_pt + By]     = field_in->getB(i, 0, 0, 1);
       (*v)[i_pt + Bz]     = field_in->getB(i, 0, 0, 2);

       (*v)[i_pt + Ax]     = field_in->getA(i, 0, 0, 0);
       (*v)[i_pt + Ay]     = field_in->getA(i, 0, 0, 1);
       (*v)[i_pt + Az]     = field_in->getA(i, 0, 0, 2);

//       (*X)[i_pt + Ex]     = fields_old->getE(i, 0, 0, 0) - (qe_h*moments_next->get_val(i, 0, 0, 0, HOMoments_currentx) + qi_h*moments_next->get_val(i, 0, 0, 1, HOMoments_currentx) - jmeanx1)*dt;
//       (*X)[i_pt + Ey]     = -2.0*(fields_next->getA(i, 0, 0, 1) - fields_old->getA(i, 0, 0, 1))/pdata->dt - fields_old->getE(i, 0, 0, 1);
//       (*X)[i_pt + Ez]     = -2.0*(fields_next->getA(i, 0, 0, 2) - fields_old->getA(i, 0, 0, 2))/pdata->dt - fields_old->getE(i, 0, 0, 2);;


     }
  }
    
//   std::cout << *v << std::endl;

}



// For standalone run
void HoLoInterface::SolveTimeStepping()
{
  
//  jfnk->DumpSolution(0);
//
//  for(int i = 0; i < steps; i++)
//    {
//      if(comm->MyPID() == 0)
//      	std::cout<< std::endl << "-----Time Step "<<i+1<<"-----" << std::endl;
//
//      jfnk->Solve();
//      jfnk->DumpSolution(i+1);
//      jfnk->PrintStatus();
//
//      jfnk->UpdateTimeStep();
//
//    }

}
