//==================================================================
//
//	AUTHOR: WILLIAM T. TAITANO
//	This is the source file for the 1D electro-static 2 species problem interface.
//
//==================================================================

//==================================================================
//	Standard header files
//==================================================================
#include 	<iostream>
//==================================================================
//	Trilinos related header files
//==================================================================

//==================================================================
//	User-defined header files
//==================================================================
#include 	"ES1D_2S_Problem_Interface.h"
#include 	"ES1D_2S_Problem.h"

//==================================================================
//
//	Constructor
//
//==================================================================
ES1D_2S_Problem_Interface::ES1D_2S_Problem_Interface( ES1D_2S_Problem 	&Problem) : _problem(Problem)
{

}
//==================================================================
//
//	Destructor
//
//==================================================================
ES1D_2S_Problem_Interface::~ES1D_2S_Problem_Interface()
{

}
//==================================================================
//
//	ComputeF function
//
//==================================================================
bool 	ES1D_2S_Problem_Interface::computeF(	const Epetra_Vector &x, Epetra_Vector &Fvec,
														FillType fillType)
{
	return _problem.evaluate(F_ONLY, &x, &Fvec, NULL);
}
//==================================================================
//
//	ComputeJacobian function
//
//==================================================================
bool 	ES1D_2S_Problem_Interface::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
{
	return _problem.evaluate(MATRIX_ONLY, &x, 0, 0);
}
//==================================================================
//
//	ComputePreconMatrix function
//
//==================================================================
bool 	ES1D_2S_Problem_Interface::computePrecMatrix(	const Epetra_Vector &x, Epetra_RowMatrix &M)
{
	cout << "Error: Problem_Interface::preconditionvector() - Use Explicit Jacobian only for this test problem!" << endl;
	throw 1;
}
//==================================================================
//
//	ComputePreconditioner function
//
//==================================================================
bool		ES1D_2S_Problem_Interface::computePreconditioner( 	const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList* precParams)
{
//	std::cout	<< "Problem_Interface::computePreconditioner" << &M << std::endl;

	cout << "No preconditioner still defined!!! ERROR!!!" << endl;
	throw 1;
//	return _problem.setupPrecOperator(&x, &M, precParams);
}
