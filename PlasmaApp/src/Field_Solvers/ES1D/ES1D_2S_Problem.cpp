//==================================================================
//
//	AUTHOR: WILLIAM T. TAITANO
//	This is the source file for the electrostatic two species problem for Trilinos
//
//==================================================================

//==============================
//	Trilinos related header-files
//==============================
#include	"NOX_Common.H"
#include 	"Epetra_Comm.h"
#include 	"Epetra_Map.h"
#include 	"Epetra_Vector.h"
#include 	"Epetra_IntVector.h"
#include 	"Epetra_Import.h"
#include 	"Epetra_CrsGraph.h"
#include 	"Epetra_CrsMatrix.h"

//==============================
//	User-defined header-files
//==============================
#include 	"ES1D_2S_Problem.h"


//==================================================================
//
//	Constructor: Creates the Epetra objects (maps and vectos) and assign inputData_in to _inputData
//	as a pointer
//
//==================================================================
ES1D_2S_Problem::ES1D_2S_Problem(	PlasmaData plasmaData_in, Epetra_Comm &Comm_in, bool restart_in,
															HOMoments *curHOMoments, HOMoments *oldHOMoments,
															FieldData *curFieldData, FieldData *oldFieldData) :
				_Comm(&Comm_in),
				_restart(restart_in)
{
	//=================================================================
	//	First thing is first, set the pointer for _plasmaData to plasmaData_in
	//=================================================================
	_plasmaData 													= &plasmaData_in;

	//=================================================================
	//	Then the HOMoments
	//=================================================================
	_curHOMoments 											= curHOMoments;
	_oldHOMoments 											= oldHOMoments;
	//=================================================================
	//	Set NumGlobal quantities
	//=================================================================
	_NumSpecies 												= plasmaData_in.nspecies;		// # of charged species
	_NumGlobalNodes											= plasmaData_in.nx;				// # of cells
	_NumGlobalFaces 											= _NumGlobalNodes+1;			// # of faces
	_NumGlobalUnknowns 									= _NumSpecies*(_NumGlobalNodes + _NumGlobalFaces) + _NumGlobalFaces;	// tot # of unknowns

	//=================================================================
	//	Set MPI quantities
	//=================================================================
	_MyPID 															= _Comm->MyPID();
	_NumProc 														= _Comm->NumProc();

	//=================================================================
	//	Now construct the map for processors to the cell (node and face)
	//=================================================================
	_StandardNodeMap 										= new Epetra_Map(_NumGlobalNodes, 0, *_Comm);
	_StandardFaceMap 											= new Epetra_Map(_NumGlobalFaces, 0, *_Comm);

	//	Get the number of nodes owned by this processor
	_NumMyNodes 												= _StandardNodeMap->NumMyElements();
	_NumMyFaces 												= _StandardFaceMap->NumMyElements();

	//=================================================================
	//												Construct an overlap node map for the problem
	//												**** CURRENTLY ONLY FOR SERIAL!!!! ****
	//=================================================================
	_OverlapNodeMap 											= new Epetra_Map(*_StandardNodeMap);
	_OverlapFaceMap 											= new Epetra_Map(*_StandardFaceMap);

	//=================================================================
	//	Now create the unknowns maps for the node and face maps
	//=================================================================
	_NumMyNodeUnknowns 									= _NumSpecies*_NumMyNodes;
	_NumMyFaceUnknowns 									= _NumSpecies*_NumMyFaces + _NumMyFaces; 	// + _NumMyFaces for the addition of field!!!
	_NumMyUnknowns 											= _NumMyNodeUnknowns + _NumMyFaceUnknowns;

	int* 		StandardMyGlobalUnknowns 				= new int[_NumMyUnknowns];
	//=================================================================
	//	now loop through to assign the cell and face id's to the standardmyglobalunknowns on each processor
	//
	//	NOTE:	currently, the ordering of unknowns is as follows:
	//
	//	unknown 	= [n_e, n_i, nu_e, nu_i, E]^{T}
	//
	//	where: size(n_e) = size(n_i) = nx and size(nu_e) = size(nu_i) = size(E) = nx+1
	//
	//=================================================================
	int 		nx 		= _plasmaData->nx;
	for (int k = 0; k < _NumSpecies; k++)
	{
		// For node
		for (int i = 0; i <_NumMyNodes; i++)
		{
			// For n_e and n_i
			StandardMyGlobalUnknowns[k*nx + i] 	= k*nx + _StandardNodeMap->GID(i);
		}
		// For face
		for (int i =0; i <_NumMyFaces; i++)
		{
			// For nu_e and nu_i
			StandardMyGlobalUnknowns[(_NumSpecies-1)*nx + k*(nx+1) + i] =
										(_NumSpecies-1)*nx + k*(nx+1) + _StandardFaceMap->GID(i);
			// For E
			StandardMyGlobalUnknowns[(_NumSpecies-1)*nx + (_NumSpecies-1)*(nx+1) + i] =
										(_NumSpecies-1)*nx + (_NumSpecies-1)*(nx+1) + _StandardFaceMap->GID(i);
		}
	}

	_StandardMap 											= new Epetra_Map(-1, _NumMyUnknowns, StandardMyGlobalUnknowns, 0, *_Comm);

	delete [] StandardMyGlobalUnknowns;

	assert(_StandardMap->NumGlobalElements() == _NumGlobalUnknowns);

	//=================================================================
	//	Currently, just allow serial processing
	//=================================================================
	_OverlapMap 												= new Epetra_Map(*_StandardMap);

	//=================================================================
	//	Construct linear objects
	//=================================================================
	_Importer 													= new Epetra_Import(*_OverlapMap, *_StandardMap);
	_nodelImporter 											= new Epetra_Import(*_OverlapNodeMap, *_StandardNodeMap);
	_faceImporter 											= new Epetra_Import(*_OverlapFaceMap, *_StandardFaceMap);
	_initialSolution 											= Teuchos::rcp(new Epetra_Vector(*_StandardMap));
	_oldSolution 												= new Epetra_Vector(*_StandardMap);

	//=================================================================
	//	Create the nodal coordinates (serial only for now)
	//=================================================================
	_xpos_node 												= Teuchos::rcp(new Epetra_Vector(*_StandardNodeMap));
	double 	dx 											= _plasmaData->dxdi;
	for (int i = 0; i < _NumGlobalNodes; i++)
	{
		(*_xpos_node)[i] 	= (i + 0.5)*dx;
	}
	//=================================================================
	//	Create the face coordinates (serial only for now)
	//=================================================================
	_xpos_face 												= Teuchos::rcp(new Epetra_Vector(*_StandardFaceMap));
	for (int i = 0;i < _NumGlobalFaces; i++)
	{
		(*_xpos_face)[i] 		= (i)*dx;
	}

	//=================================================================
	//	Initialize the solution
	//=================================================================
	initializeSoln(curHOMoments, oldHOMoments, curFieldData, oldFieldData);
}














//==================================================================
//
//	Destructor
//
//==================================================================
ES1D_2S_Problem::~ES1D_2S_Problem()
{
	delete 		_oldSolution;
	delete 		_nodelImporter;
	delete 		_faceImporter;
	delete 		_Importer;
	delete 		_OverlapMap;
	delete 		_StandardMap;
	delete 		_StandardNodeMap;
	delete 		_StandardFaceMap;
	delete 		_OverlapNodeMap;
	delete 		_OverlapFaceMap;

}


























//==================================================================
//
//	Reset function: set the old solution to the new solution
//	Probably won't be using this function since the LO system is not strictly time integrated
//
//==================================================================
void ES1D_2S_Problem::reset(const Epetra_Vector &x)
{
	*_oldSolution 		= x;
}

















//==================================================================
//
//	Set initialSolution to desired initial condition
//
//==================================================================
void ES1D_2S_Problem::initializeSoln(HOMoments *curHOMoments, HOMoments *oldHOMoments,
															FieldData *curFieldData, FieldData *oldFieldData)
{
	Epetra_Vector 			&soln 			= *_initialSolution;
	Epetra_Vector	 			&xpos_node	= *_xpos_node;
	Epetra_Vector 			&xpos_face	= *_xpos_face;
	//==================================================================
	//	Extract geometrical info from PlasmaData
	//==================================================================
	int 							nx 				= _plasmaData->nx;
	//==================================================================
	//	Set initial condition to the soln vector
	//==================================================================
	if (_restart)
	{
		string 					name 			= "restartVec";
		Epetra_Vector 		*tmp_vec 		= NULL;
		soln 											= *tmp_vec;
	}
	else
	{
		for (int k = 0; k < _NumSpecies; k++)
		{
			// Loop over cells
			for (int i = 0; i < xpos_node.MyLength(); i++)
			{
				//	n
				soln[k*nx + i] 	= curHOMoments->get_val(i, 0, 0, k, HOMoments_charge);
			}	// Loop over cells
			// Loop over faces
			for (int i = 0; i< xpos_face.MyLength(); i++)
			{
				if (i == xpos_face.MyLength() - 1) // If right boundary face
				{
					//	nu
					soln[(_NumSpecies-1)*nx + k*(nx+1) + i] 	= curHOMoments->get_val(0, 0, 0, k, HOMoments_currentx);
					//	E
					soln[(_NumSpecies-1)*nx + (_NumSpecies-1)*(nx+1) + i] = curFieldData->getE(0, 0, 0, 0);
				}
				else	// other faces
				{
					//	nu
					soln[(_NumSpecies-1)*nx + k*(nx+1) + i] 	= curHOMoments->get_val(i, 0, 0, k, HOMoments_currentx);
					//	E
					soln[(_NumSpecies-1)*nx + (_NumSpecies-1)*(nx+1) + i] = curFieldData->getE(i, 0, 0, 0);
				}
			}	// Loop through cell faces
		}	// Loop through species
	}
	//==================================================================
	//	Now set the oldSolution
	//==================================================================
	for (int k = 0; k < _NumSpecies; k++)
	{
		// Loop over cells
		for (int i = 0; i < xpos_node.MyLength(); i++)
		{
			//	n
			(*_oldSolution)[k*nx + i] 	= oldHOMoments->get_val(i, 0, 0, k, HOMoments_charge);
		}	// Loop over cells
		// Loop over faces
		for (int i = 0; i< xpos_face.MyLength(); i++)
		{
			if (i == xpos_face.MyLength() - 1) // If right boundary face
			{
				//	nu
				(*_oldSolution)[(_NumSpecies-1)*nx + k*(nx+1) + i] 	= oldHOMoments->get_val(0, 0, 0, k, HOMoments_currentx);
				//	E
				(*_oldSolution)[(_NumSpecies-1)*nx + (_NumSpecies-1)*(nx+1) + i] = oldFieldData->getE(0, 0, 0, 0);
			}
			else	// other faces
			{
				//	nu
				(*_oldSolution)[(_NumSpecies-1)*nx + k*(nx+1) + i] 	= oldHOMoments->get_val(i, 0, 0, k, HOMoments_currentx);
				//	E
				(*_oldSolution)[(_NumSpecies-1)*nx + (_NumSpecies-1)*(nx+1) + i] = oldFieldData->getE(i, 0, 0, 0);
			}
		}	// Loop through cell faces
	}	// Loop through species
}



























//=================================================================
//
//	The evaluation of residual fnction
//
//	NOTE: I think there's something very very inefficient about this routine. I mean, how can it make any
//	sense to evaluate the "WHOLE" residual function vector??? It makes more sense to just perturb the
//	off-diagonal blocks and sum up to the current running J*v vector.... doesn't??? Ask how this function
//	perturbation is performed exactly...
//
//=================================================================
bool 		ES1D_2S_Problem::evaluate(	FillType fType,
															const Epetra_Vector *soln,
															Epetra_Vector *tmp_rhs,
															Epetra_RowMatrix *tmp_matrix)
{
	if (fType == MATRIX_ONLY)
	{
		cout << "This problem only works for finite-difference or matrix-free based Jacobians." <<
				"No analytic jacobian fill available !!" << endl;
		throw "Problem Error";
	}
	else
	{
		// Initialize the rhs residual vector
		_rhs 		= new Epetra_Vector(*_OverlapMap);
	}

	printf("came into evaluate ? \n");
	//=================================================================
	//	Create the overlapped solution and position vectors
	//=================================================================
	Epetra_Vector 			u(*_OverlapMap);
	Epetra_Vector 			uold(*_OverlapMap);
	Epetra_Vector 			xvec_node(*_OverlapNodeMap);
	Epetra_Vector 			xvec_face(*_OverlapFaceMap);

	//=================================================================
	//	Export solution to overlap vector
	//=================================================================
	uold.Import(*_oldSolution, *_Importer, Insert);
	xvec_node.Import(*_xpos_node, *_nodelImporter, Insert);
	xvec_face.Import(*_xpos_face, *_faceImporter, Insert);
	u.Import(*soln, *_Importer, Insert);
	printf("passed all the imports ? \n");
	//=================================================================
	//	Declare required variables
	//=================================================================
	int 			i;
	int 			OverlapNumMyNodes 	= _OverlapNodeMap->NumMyElements();
	int 			OverlapNumMyFaces	= _OverlapFaceMap->NumMyElements();

	int 			OverlapMinMyNodeGID;
	int 			OverlapMinMyFaceGID;

	int 			row1, row2;
	double 	dx 				= _plasmaData->dxdi;
	double 	dxRecip 			= 1.0/dx;
	double 	dxsRecip 		= dxRecip*dxRecip;

	double 	dt 					= _plasmaData->dt;
	double 	dtRecip 			= 1.0/dt;
//	double 	eps0 			= epsilon0;
	int 			nx 				= _plasmaData->nx;

	//=================================================================
	//	Declare some cell center and face solution quantities
	//=================================================================
	double 	n_i, n_old_i;
	double 	n_im1, n_old_im1;
	double 	n_face, n_face_old, n_face_half;
	double 	nu_ip1h, nu_im1h, E_ip1h;
	double 	nu_old_ip1h, E_old_ip1h;
	double 	E_ip1h_half;
	double 	S2_i, S2_im1, n_HO_i, n_HO_im1;
	double 	S2_old_i, S2_old_im1, n_HO_old_i, n_HO_old_im1;
	double 	S_tld_i, S_tld_im1;
	double 	gamma_nu;

	printf("passed all the declaration?\n");
	//=================================================================
	//	Specify the overlap (ghost cell) ids
	//	NOTE: My code is serial but I'm just including here in case of future application when we go parallel
	//=================================================================
	if (_MyPID == 0)			// If host core
	{
		OverlapMinMyNodeGID 	= _StandardNodeMap->MinMyGID();
		OverlapMinMyFaceGID 		= _StandardFaceMap->MinMyGID();
	}
	else							// if other
	{
		OverlapMinMyNodeGID 	= _StandardNodeMap->MinMyGID() - 1;
		OverlapMinMyFaceGID 		= _StandardFaceMap->MinMyGID() - 1;
	}
	printf("passed all the overlapmin/maxmynodefacegid ? \n");
	//=================================================================
	//	Zero out the object (RHS) that will be filled
	//=================================================================
	_rhs->PutScalar(0.0);


	printf("passed putscalar ? \n");




	//=================================================================
	//
	//								NOW CALCULATE THE CONSISTENCY TERMS!
	//
	//=================================================================

	ConsistencyTermES 		consistencyTermES			= ConsistencyTermES(_plasmaData);
	int nspecies 		= _plasmaData->nspecies;
	for (int s = 0; s < nspecies; s ++) {

		consistencyTermES.consistencyCalcContinuity(s,_curHOMoments,_oldHOMoments);
		printf("passed s = %d consistencycalcontinuity \n",i);
		consistencyTermES.consistencyCalcMomentum(s,_curHOMoments,_oldHOMoments,_curFieldData,_oldFieldData);
		printf("passed s = %d consistencycalmomentum \n",i);
	}
	printf("passed consistencyterm calc?\n");



	//=================================================================
	//
	//											BELOW IS FOR DENSITY RESIDUAL: 		DONE!!!!
	//
	//=================================================================
	//=================================================================
	//	Loop over the cells and calculate density equation residuals
	//=================================================================
	for (int k = 0; k < _NumSpecies; k++)
	{
		//	Declare mass
		double ms 					= mass_e*_plasmaData->mspecies[k];
		//	Declare offset value
		int 				offset 		= k*nx;
		int				offsetnu 	= (_NumSpecies-1)*nx + k*(nx+1);
		for (int i = 0; i < OverlapNumMyNodes; i++)
		{
			printf("Now at cell loop i = %d \n",i);
			//=================================================================
			//	Extract the cell center solutions and cell face flux
			//=================================================================
			n_i 						= u[offset + i];
			n_old_i 					= uold[offset + i];

			nu_ip1h 				= u[offsetnu + i+1];
			nu_im1h 				= u[offsetnu + i];
			//=================================================================
			//	Build residual for density
			//=================================================================
			(*_rhs)[offset + i] 	= ms*(dtRecip*(n_i - n_old_i) + dxRecip*(nu_ip1h - nu_im1h));
		}
	}









	//=================================================================
	//
	//						BELOW IS FOR MOMENTUM AND AMPERE RESIDUAL:
	//									CALCULATE CELL FACE DENSITY AT HALF TIME
	//									STORE HO SOLUTION THE PROBLEM CLASS AS A MEMBER VALUE
	//									TO EXTRACT THE STRESS TENSOR
	//
	//=================================================================

	//=================================================================
	//														MOMENTUM
	//	Loop over the faces and calculate momentum and ampere equation residual:			DONE!!!
	//=================================================================
	// Define some offsets
	int 			offsetE 			= (_NumSpecies-1)*nx + (_NumSpecies-1)*(nx+1);
	for (int i = 0; i < OverlapNumMyFaces; i++)
	{
		printf("Now at face loop i = %d \n",i);
		// Momentum Residual Calculation
		for (int k = 0; k < _NumSpecies; k++)
		{
			//=================================================================
			//	Calculate charge to mass ratio
			//	and extract mass of species
			//=================================================================
			double 		qs2ms		= qe2me*_plasmaData->qspecies[k]/_plasmaData->mspecies[k];
			double 		ms 			= mass_e*_plasmaData->mspecies[k];
			//=================================================================
			//	Calculate some offsets
			//=================================================================
			int 				offset 		= k*nx;
			int				offsetnu 	= (_NumSpecies-1)*nx + k*(nx+1);
			int 				offsetE 		= (_NumSpecies-1)*nx + (_NumSpecies-1)*(nx+1);
			//=================================================================
			//	Extract the cell face solution and build Residual for momentum and ampere
			//=================================================================
			nu_ip1h 						= u[offsetnu + i];
			nu_old_ip1h	 				= uold[offsetnu + i];
			E_ip1h	 						= u[offsetE + i];
			E_old_ip1h 					= uold[offsetE + i];
			if (i == 0)				// Left face
			{
				n_i 							= u[0];
				n_im1 						= u[nx];
				n_old_i 						= uold[0];
				n_old_im1 					= uold[nx];

				n_face 						= 0.5*(n_i + n_im1);
				n_face_old 				= 0.5*(n_old_i + n_old_im1);

				n_face_half 				= 0.5*(n_face + n_face_old);

				E_ip1h_half 				= 0.5*(E_ip1h + E_old_ip1h);

				S2_i 							= _curHOMoments->get_val(0,0,0,k,HOMoments_S2xx);
				S2_im1 						= _curHOMoments->get_val(nx,0,0,k,HOMoments_S2xx);
				n_HO_i 						= _curHOMoments->get_val(0,0,0,k,HOMoments_charge);
				n_HO_im1 					= _curHOMoments->get_val(nx,0,0,k,HOMoments_charge);

				gamma_nu 				= consistencyTermES.get_val(i,0,0,0,consistencyterm_momentum);

			}
			else if (i == nx)		// Right face
			{
				n_i 							= u[0];
				n_im1 						= u[nx];
				n_old_i 						= uold[0];
				n_old_im1 					= uold[nx];

				n_face 						= 0.5*(n_i + n_im1);
				n_face_old 				= 0.5*(n_old_i + n_old_im1);

				n_face_half 				= 0.5*(n_face + n_face_old);

				E_ip1h_half 				= 0.5*(E_ip1h + E_old_ip1h);

				S2_i 							= _curHOMoments->get_val(0,0,0,k,HOMoments_S2xx);
				S2_im1 						= _curHOMoments->get_val(nx,0,0,k,HOMoments_S2xx);
				n_HO_i 						= _curHOMoments->get_val(0,0,0,k,HOMoments_charge);
				n_HO_im1 					= _curHOMoments->get_val(nx,0,0,k,HOMoments_charge);
				gamma_nu 				= consistencyTermES.get_val(0,0,0,0,consistencyterm_momentum);
			}
			else
			{
				n_i 							= u[offset+i];
				n_im1 						= u[offset+i-1];
				n_old_i 						= uold[offset+i];
				n_old_im1 					= uold[offset+i-1];

				n_face 						= 0.5*(n_i + n_im1);
				n_face_old 				= 0.5*(n_old_i + n_old_im1);

				n_face_half 				= 0.5*(n_face + n_face_old);

				E_ip1h_half 				= 0.5*(E_ip1h + E_old_ip1h);

				S2_i 							= _curHOMoments->get_val(i,0,0,k,HOMoments_S2xx);
				S2_im1 						= _curHOMoments->get_val(i-1,0,0,k,HOMoments_S2xx);
				n_HO_i 						= _curHOMoments->get_val(i,0,0,k,HOMoments_charge);
				n_HO_im1 					= _curHOMoments->get_val(i-1,0,0,k,HOMoments_charge);
				gamma_nu 				= consistencyTermES.get_val(i,0,0,0,consistencyterm_momentum);
			}
			S_tld_i 					= 0.5*(S2_i + S2_im1)/n_HO_i;
			S_tld_im1 				= 0.5*(S2_old_i + S2_old_im1)/n_HO_im1;
			//=================================================================
			//	Actually compute the residual!!!		<------------------------------------------------- Make sure to add contribution from consistency term as well
			//=================================================================
			(*_rhs)[i] 		= ms*(	2*dtRecip*(nu_ip1h - nu_old_ip1h) +
											dxRecip*(n_i*S_tld_i - n_im1*S_tld_im1) -
											qs2ms*n_face_half*E_ip1h_half) -
											gamma_nu*n_face;
		}








		//=================================================================
		//																AMPERE
		// Ampere Residual calculate: 	DONE!!!
		//=================================================================
		//	Calculate the average current first
		//=================================================================
		double 	j_avg = 0.;
		for (int k = 0; k < _NumSpecies; k++)
		{
			double qs 			= _plasmaData->qspecies[k]*qe;
			int		offsetnu 	= (_NumSpecies-1)*nx + k*(nx+1);
			for (int i = 0; i < OverlapNumMyFaces - 1; i++)					// Note the OverlapNumMyFaces-1 is chosen to not include the right boundary value
			{
				nu_ip1h 	= u[offsetnu+i];
				j_avg += qs*nu_ip1h;
			}
		}
		j_avg 	/= double(nx);

		// 	Calculate the offset
		int 				offsetE 		= (_NumSpecies-1)*nx + (_NumSpecies-1)*(nx+1);
		//	Extract some solution
		E_ip1h 							= u[offsetE + i];
		E_old_ip1h 					= uold[offsetE + i];
		// actually compute the residual below
		(*_rhs)[offsetE + i]	= epsilon_naught*dtRecip*(E_ip1h - E_old_ip1h);	// the eps0*(E - E_old)/dt part
		for (int k = 0; k < _NumSpecies; k++) 	// the qs*nu_half part
		{
			//=================================================================
			//	Calculate some offsets
			//=================================================================
			int 				offset 		= k*nx;
			int				offsetnu 	= (_NumSpecies-1)*nx + k*(nx+1);
			//=================================================================
			//	Calculate the residual
			//=================================================================
			double 	qes 				= 		_plasmaData->qspecies[k]*qe;
			nu_ip1h 						= 		u[offsetnu + i];
			(*_rhs)[offsetE + i] 			+=	qes*nu_ip1h;
		}
		(*_rhs)[offsetE + i] 	-= j_avg; // the - j_avg part
	}


	//=================================================================
	//	Sync-up processors to be safe (although only serial, it's a good practice!)
	//=================================================================
	_Comm->Barrier();

	//=================================================================
	//	Do an assemble for overlap nodes
	//=================================================================
	tmp_rhs->Export(*_rhs, *_Importer, Insert);

	//=================================================================
	//	Delete the _rhs since we've already exported the solution to the tmp_rhs
	//=================================================================
	delete _rhs;		//<-- Shazam! Disappear bitch!

	return true;
}










//=================================================================
//
//	Function that sets up the preconditioner operator
//
//=================================================================
bool	ES1D_2S_Problem::setupPrecOperator(const Epetra_Vector *x, Epetra_Operator *Prec, Teuchos::ParameterList* precParams)
{
//	static bool setup_done = false;
//
//	if (setup_done)
//	{
//		return true;
//	}
////	bool setup_done = false;
//	MLPreconditioner	*MLPrec = dynamic_cast<MLPreconditioner *>(Prec);
//
////	printf("MLPrec->setup() entering \n");
//	MLPrec->setup(x);
////	printf("MLPrec->setup() success! \n");
//	setup_done 	= true;
	return true;
}

//=================================================================
//
//	Function that creates the preconditioning matrix/operator
//
//=================================================================
void ES1D_2S_Problem::createMatrix( const Epetra_Vector *x,std::vector<Epetra_CrsMatrix *> &A)
{

}










//=================================================================
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//				BELOW ARE GET FUNCTIONS TO ACCESS PRIVATE MEMBER VALUES
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//=================================================================












//=================================================================
//
//	Get initialSolution function
//
//=================================================================
Teuchos::RCP<Epetra_Vector> 	ES1D_2S_Problem::getSolution()
{
	return _initialSolution;
}

//=================================================================
//
//	Get Jacobian function
//
//=================================================================
Teuchos::RCP<Epetra_CrsMatrix> 	ES1D_2S_Problem::getJacobian()
{
	if (Teuchos::is_null(_A) )
	{
		return _A;
	}
	else
	{
		cout << "No valid Jacobian matrix for this problem. This is likely the"
				<< "result of overlapping NODES rather than ELEMENTS. \n" << endl;
		throw "MEB_Problem Error";
	}
}

//=================================================================
//
//	GetMesh function
//
//=================================================================
Teuchos::RCP<Epetra_Vector> 		ES1D_2S_Problem::getMesh()
{
	return _xpos_node;
}

//=================================================================
//
//	GetOldSolution function
//
//=================================================================
Epetra_Vector& 							ES1D_2S_Problem::getOldSoln()
{
	return *_oldSolution;
}

//=================================================================
//
//	Function that returns the Epetra_Comm
//
//=================================================================
Epetra_Comm& 							ES1D_2S_Problem::getComm()
{
	return *_Comm;
}

//=================================================================
//
//	Function that return the _StandardMap
//
//=================================================================
Epetra_Map& 							ES1D_2S_Problem::getStandardMap()
{
	return *_StandardMap;
}

//=================================================================
//
//	Function that return the _OverlapMap
//
//=================================================================
Epetra_Map& 							ES1D_2S_Problem::getOverlapMap()
{
	return *_OverlapMap;
}

//=================================================================
//
//	Function that returns the _StandardNodeMap
//
//=================================================================
Epetra_Map& 							ES1D_2S_Problem::getStandardNodeMap()
{
	return *_StandardNodeMap;
}

//=================================================================
//
//	Function that returns the _StandardFaceMap
//
//=================================================================
Epetra_Map& 							ES1D_2S_Problem::getStandardFaceMap()
{
	return *_StandardFaceMap;
}

//=================================================================
//
//	Function that returns the _OverlapNodeMap
//
//=================================================================
Epetra_Map& 							ES1D_2S_Problem::getOverlapNodeMap()
{
	return *_OverlapNodeMap;
}

//=================================================================
//
//	Function that returns the _OverlapFaceMap
//
//=================================================================
Epetra_Map& 							ES1D_2S_Problem::getOverlapFaceMap()
{
	return *_OverlapFaceMap;
}
//=================================================================
//
//	Function that returns the _Importer
//
//=================================================================
Epetra_Import& 							ES1D_2S_Problem::getImporter()
{
	return *_Importer;
}

//=================================================================
//
//	Function that returns the _nodelImporter
//
//=================================================================
Epetra_Import& 							ES1D_2S_Problem::getnodelImporter()
{
	return *_nodelImporter;
}

//=================================================================
//
//	Function that returns the _faceImporter
//
//=================================================================
Epetra_Import& 							ES1D_2S_Problem::getfaceImporter()
{
	return *_faceImporter;
}































