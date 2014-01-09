//==================================================================
//
//	AUTHOR: WILLIAM T. TAITANO
//	This is the source file for the LO JFNK solver for the electro-static two fluid system
//
//==================================================================

//==================================================================
//	Trilinos related header-files
//==================================================================
#include "Teuchos_Comm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"
//#include "Teuchos_Version.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

//NOX related header files
#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_Vector.H"
#include "NOX_Common.H"

//Epetra comm and import related header files
#include "Epetra_Comm.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"

// Other Epetra related header files
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
//#include "Epetra_Version.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Time.h"

// Linear solver related header fiels
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"
//==================================================================
//	User-defined header-files
//==================================================================
#include 	"LOSolverJFNK.h"
#include 	"PlasmaUtility.h"

//==================================================================
//
//	Constructor
//
//==================================================================
LOSolverJFNK::LOSolverJFNK(	PlasmaData* pDataIn,
												HOMoments* curHOMoments, HOMoments* oldHOMoments,
												FieldData* 		curFieldData, 		FieldData* oldFieldData)
{
	//==================================================================
	//	Fill parent class (LOSolver) members
	//==================================================================
	this->pData 				= pDataIn;
	curLOMomentsES		= new LOMomentsES(pDataIn, curHOMoments);
	oldLOMomentsES 		= new LOMomentsES(pDataIn, oldHOMoments);
	consistencyTermES 	= new ConsistencyTermES(pDataIn);
	int 	nspecies 			= pData->nspecies;
	for (int s = 0; s < nspecies; s++)
	{
		consistencyTermES->consistencyCalcContinuity(s, curHOMoments, oldHOMoments);
		consistencyTermES->consistencyCalcMomentum(s, curHOMoments, oldHOMoments, curFieldData, oldFieldData);
	}
	curIntField 					= new FieldDataCPU();
	curIntField->allocate(pDataIn);
	for (int i = 0; i < pData->nx;i++)
	{
		curIntField->getE(i,0,0,0) 	= curFieldData->getE(i,0,0,0);
	}
}



//==================================================================
//
//	Destructor
//
//==================================================================

LOSolverJFNK::~LOSolverJFNK()
{
	delete(curLOMomentsES);
	delete(oldLOMomentsES);
	delete(consistencyTermES);
}



//==================================================================
//
//	Member Function: solve
//
//==================================================================
int 		LOSolverJFNK::solve(
						HOMoments* curHOMoments, 			// updates values
						HOMoments* oldHOMoments,
						FieldData* 		curFieldData,
						FieldData* 		oldFieldData)
{
	//==================================================================
	//	Declare variable
	//==================================================================
	int 	k_inner;	//<------------------------------------------------------------------------------------------------------------------------------	May not need to use
	//==================================================================
	//	Choose which solver based on the number of species specified by user (currently only 2 species)
	//==================================================================
	k_inner 		= jfnkSolverTwoSpecies(curHOMoments, oldHOMoments, curFieldData, oldFieldData);

	//==================================================================
	//	A return statement
	//==================================================================
	return 		k_inner;
}
//==================================================================
//
//	Member Function: jfnkSolverTwoSpecies
//						This is where all the Trilinos magic and interfacing with it happens.
//											RELEASE THE TRILINOS BEAST!!!
//
//==================================================================
int 		LOSolverJFNK::jfnkSolverTwoSpecies(
													HOMoments* curHOMoments,
													HOMoments* oldHOMoments,
													FieldData* 		curFieldData,
													FieldData* 		oldFieldData)
{
	//==================================================================
	//
	//															Step 0:
	//													Setup MPI Comm <-- A BIG MUST!!!
	//
	//==================================================================
	// The current implementation is serial, initialize serial Comm
	Epetra_SerialComm Comm;
	printf("MPIComm is off, currently running serial \n");
	//	Get the process ID and the total number of processors
	int 			MyPID 		= Comm.MyPID();
	int 			NumProc 	= Comm.NumProc();

	printf("passed step 0 \n");

	//==================================================================
	//
	//															Step 1:
	//											Initialize simulation parameters
	//
	//==================================================================
	//	Define some simulation parameters
	bool 		verbose 			= false; 						// verbosity flag for Trilinos solver
	bool 		restart_flag 		= false; 						// restart of solution
	int 			precon_flag 		= true;							// flag to use physics based preconditioner or not


	printf("passed step 1 \n");

	//==================================================================
	//
	//															Step 2:
	//								Create and setup the problem (nonlinear)
	//
	//==================================================================


	//	Setup Teuchos command line processor
	Teuchos::CommandLineProcessor 		clp(false);

	//	Setup the problem
	ES1D_2S_Problem			 					Problem = ES1D_2S_Problem(*pData,Comm,restart_flag,
																										curHOMoments, oldHOMoments,
																										curFieldData, oldFieldData);
	// Get the vector from the problem
	Teuchos::RCP<Epetra_Vector> 				soln 		= Problem.getSolution();
	NOX::Epetra::Vector 							noxSoln(soln, NOX::Epetra::Vector::CreateView);

	printf("passed step 2 \n");
	//==================================================================
	//
	//															Step 3:
	//										Setup nonlinear solver parameters
	//
	//==================================================================

	// Create the top level parameter list
	Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
			Teuchos::rcp(new Teuchos::ParameterList);
	Teuchos::ParameterList& nlParams = *nlParamsPtr.get();

	// Set the nonlinear solver method
	nlParams.set("Nonlinear Solver", "Line Search Based");
	//nlParams.set("Nonlinear Solver", "Trust Region Based");

	// Set the printing parameters in the "Printing" sublist
	Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
	if (verbose && MyPID == 0)	// If verbose option is true and host core
	{
		printParams.set("MyPID", MyPID);
		printParams.set("Output Precision", 3);
		printParams.set("Output Processor", 0);
//		printParams.set("Output Information",
//				NOX::Utils::OuterIteration +
//				NOX::Utils::OuterIterationStatusTest +
//				NOX::Utils::InnerIteration +
//				NOX::Utils::Parameters +
//				NOX::Utils::Details +
//				NOX::Utils::Warning);
		  printParams.set("Output Information",
		                  NOX::Utils::OuterIterationStatusTest +
		                  NOX::Utils::OuterIteration +
		                  NOX::Utils::Warning +
		                  NOX::Utils::LinearSolverDetails);
	}
	else
	{
		printParams.set("Output Information", NOX::Utils::Error);
	}
	// Create a print handle object
	NOX::Utils utils(printParams);

	// Sublist for line search
	Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
	searchParams.set("Method", "Full Step");

	// Sublist for direction
	Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
	dirParams.set("Method", "Newton");
	Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
	newtonParams.set("Forcing Term Method", "Constant");

	// Sublist for linear solver for the Newton method
	Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
	lsParams.set("Aztec Solver", "GMRES");
	lsParams.set("Max Iterations", 100);
	lsParams.set("Tolerance", 1e-4);
	lsParams.set("Output Frequency", 100);

	// Setup preconditioning flag parameter
//	if (precon_flag == 0)			// No preconditioner
//	{
		lsParams.set("Preconditioner", "None");
//	}
//	else if (precon_flag == 1)	// User-defined
//	{
//		lsParams.set("Preconditioner", "User Defined");
//		lsParams.set("Preconditioner Reuse Policy", "Recompute");
//		lsParams.set("Max Age Of Prec", 5);
//	}
//	else	// Finite-difference
//	{
//		lsParams.set("Preconditioner", "AztecOO");
//	}

	printf("passed step 3 \n");
	//==================================================================
	//
	//															Step 4:
	//								Create and setup nonlinear problem interface
	//							(this will serve as the interface between the application and Trilinos!!!)
	//
	//==================================================================

	Teuchos::RCP<ES1D_2S_Problem_Interface> interface 			= Teuchos::rcp(new ES1D_2S_Problem_Interface(Problem));
	// Matrix-Free (Epetra_Operator)
	Teuchos::RCP<NOX::Epetra::MatrixFree> 						MF 	= Teuchos::rcp(new NOX::Epetra::MatrixFree(
																									printParams, interface, noxSoln));

	// NOTE: Currently, no graph is created in the problem and FD Jac cannot be formed (will get seg-fault)
	// Finite Difference (Epetra_RowMatrix)
//	Teuchos::RCP<NOX::Epetra::FiniteDifference> 				FD 	= Teuchos::rcp(new NOX::Epetra::FiniteDifference(
//																									printParams, interface, noxSoln));

	Teuchos::RCP<NOX::Epetra::Interface::Required> 			iReq 	= interface;
	Teuchos::RCP<NOX::Epetra::Interface::Jacobian>			iJac 	= MF;
	Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>	iPrec = interface;

	//==================================================================
	//
	//															Step 4.1:
	//												Create pointers to the PBP
	//
	//
	//==================================================================

	printf("passed step 4.1 \n");
	//==================================================================
	//
	//															Step 5:
	//									Create and setup linear problem AztecOO for Newton
	//								(this will actually setup the linear solver to inver J*du = -F
	//
	//==================================================================

	// Create the Linear System
	Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys;
	linSys =
		Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
		iReq, noxSoln));

	printf("passed step 5 \n");

	//==================================================================
	//
	//															Step 6:
	//								Create the Group (trying to still figure out WTF this is...)
	//
	//==================================================================

	// Create the Group
	Teuchos::RCP<NOX::Epetra::Group> grpPtr =
				Teuchos::rcp(new NOX::Epetra::Group(printParams,
				iReq,
				noxSoln,
				linSys));
	// Get pointer for the group
	NOX::Epetra::Group& grp = *(grpPtr.get());

	printf("passed step 6 \n");
	//==================================================================
	//
	//															Step 7:
	//											Create convergence criteria
	//
	//==================================================================

	// Create the convergence tests
	Teuchos::RCP<NOX::StatusTest::NormF> absresid =
			Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-5, NOX::StatusTest::NormF::Unscaled));
	Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
			Teuchos::rcp(new NOX::StatusTest::MaxIters(25));
	Teuchos::RCP<NOX::StatusTest::FiniteValue> finiteval =
			Teuchos::rcp(new NOX::StatusTest::FiniteValue());
	Teuchos::RCP<NOX::StatusTest::Combo> combo =
			Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
	combo->addStatusTest(absresid);
	combo->addStatusTest(maxiters);
	combo->addStatusTest(finiteval);

	printf("passed step 7 \n");
	//==================================================================
	//
	//															Step 8:
	//											Create (Nonlinear) solver method
	//
	//==================================================================

	// Create the method
	Teuchos::RCP<NOX::Solver::Generic> solver =
			NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);

	printf("passed step 8 \n");
	//==================================================================
	//
	//															Step 9:
	//											Begin the nonlinear solve
	// 	(no time-integration is required as LO system is not explicitly integrated form time-step to time-step)
	//
	//==================================================================

	//================================
	//	Initialize the status of NOX solver
	//================================
	NOX::StatusTest::StatusType status = NOX::StatusTest::Unconverged;
	//================================
	//	Reset status of NOX solver before solve
	//================================
	status 		= NOX::StatusTest::Unconverged;
	//================================
	// 	Solve nonlinear system
	//================================
	printf("passed all the way before solver->sove() \n");
	status 		= solver->solve();

	printf("passed step 9 \n");
	//==================================================================
	//
	//															Step 10:
	//												Any post processing comes here
	//							Here, extract the trilinos solution and put into the LOMoment class values
	//
	//==================================================================

	printf("passed step 10 \n");
	//==================================================================
	//
	//															Step 11:
	//											Finally close out MPI (if MPI is used)
	//
	//==================================================================

	printf("passed step 11 \n");
	return 0;		// <-------------------------------------------------------------Really, I want to return the k_inner for Newton but I'll figure out later
}
