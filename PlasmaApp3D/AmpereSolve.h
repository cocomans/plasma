#ifndef AMPERE_SOLVE_H
#define AMPERE_SOLVE_H
/*-------------------------------------------------------------------------*/
/**
  @file		AmpereSolve.h
  @author	J. Payne
  @date		2012
  @brief	Depricated Ampere solve function and residual calculation

	@depricated This file contains a function to solve amperes equation in 3D and
	a residual output. This is a depricated function and should no longer
	be used. The AmpereSolver class should be used instead.

*/
/*--------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>


class PlasmaData;
class HOMoments;
class FieldDataCPU;


void AmpereSolve(PlasmaData* pdata,
		FieldDataCPU* fields_next,
		FieldDataCPU* fields_old,
		HOMoments* moments);

double Calc_Residual(PlasmaData* pdata,
		FieldDataCPU* fields_next,
		FieldDataCPU* fields_old,
		HOMoments* moments);


#endif /* AMPERE_SOLVE_H */
