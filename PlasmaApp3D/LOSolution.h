#ifndef MOMENT_SOLUTION_H
#define MOMENT_SOLUTION_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>



class MomentSolution
{
public:


private:

  double* R;           // Cell density (Oth moment)
  double* RU;          // Face momentum (1st moment)
  double* S2;          // Cell stress (2nd moment)
  double* S3;          // Cell third moment

  double* U;           // Cell velocity
  double* T;           // Cell temperature
  double* P;           // Cell pressure

  double* E_TH;        // Cell thermal energy moment
  double* E_TOT;
  double* E;            // Face electric field

  double* R_f;         // Face density
  double* RU_c;        // Cell momentum


};








#endif /* MOMENT_SOLUTION_H */
