/*-------------------------------------------------------------------------*/
/**
	@file		FieldData.h
	@author	J. Payne
	@date		12/21/2012
	@brief	Declares the FieldData class, a pure virtual class used to encapsulate
	the device specific and optimization specific implementations for field data
	interpolation.


*/
/*--------------------------------------------------------------------------*/
#ifndef FIELD_DATA_H
#define FIELD_DATA_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <gnuplot_i.h>
#include "PlasmaData.h"



typedef struct {
	realkind vals[6]; // Ex, Ey, Ez, Bx, By, Bz
} FieldValues;

/*-------------------------------------------------------------------------*/
/**
	@enum FieldData_deriv
	Enum for choosing derivative of a field data interpolation.

	\var FieldData_deriv::FieldData_deriv_f
	Return the result of the function itself, no derivitive.

	\var FieldData_deriv::FieldData_deriv_dfdx
	Return the 1st derivitive in the x direction.

	\var FieldData_deriv::FieldData_deriv_dfdy
	Return the 1st derivitive in the y direction.

	\var FieldData_deriv::FieldData_deriv_dfdz
	Return the 1st derivitive in the z direction.
*/
/*--------------------------------------------------------------------------*/
enum FieldData_deriv
{
	FieldData_deriv_f = 0,
	FieldData_deriv_dfdx = 1,
	FieldData_deriv_dfdy = 2,
	FieldData_deriv_dfdz = 3
};


enum FieldData_type
{
	FieldData_cpu = 0,
	FieldData_cpu_aos = 1,
	FieldData_gpu = 2,
	FieldData_mic = 3,
	FieldData_mic_aos = 4
};

#endif /* FIELD_DATA_H */
