/*=========================================================================
                                                                                
Copyright (c) 2011, Los Alamos National Security, LLC

All rights reserved.

Copyright 2011. Los Alamos National Security, LLC. 
This software was produced under U.S. Government contract DE-AC52-06NA25396 
for Los Alamos National Laboratory (LANL), which is operated by 
Los Alamos National Security, LLC for the U.S. Department of Energy. 
The U.S. Government has rights to use, reproduce, and distribute this software. 
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,
EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
If software is modified to produce derivative works, such modified software 
should be clearly marked, so as not to confuse it with the version available 
from LANL.
 
Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
-   Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer. 
-   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution. 
-   Neither the name of Los Alamos National Security, LLC, Los Alamos National
    Laboratory, LANL, the U.S. Government, nor the names of its contributors
    may be used to endorse or promote products derived from this software 
    without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
                                                                                
=========================================================================*/

///////////////////////////////////////////////////////////////////////
//
// Tridiagonal matrix algorithm, simplified form of Gaussian elimination
// with backsubstitution.  Forward sweep solves for x1 in first equation
// and substitutes that into second equation, thus eliminating one unknown.
// This continues to the last equation when back substitution gets the result
//
// Thomas algorithm accomplishes this by altering the coefficients
//
///////////////////////////////////////////////////////////////////////

#include <math.h>
#include <iostream>

using namespace std;

double genRandom()
{
  return (double) rand() / (double) RAND_MAX;
}

void solveMatrix2(
                int n,          // Number of equations
                double *a,      // Sub diagonal values
                double *b,      // Main diagonal values
                double *c,      // Super diagonal values
                double *d,      // Right hand side
                double *x)      // Answer
{
  c[0] /= b[0];
  d[0] /= b[0];

  // Forward sweep to modify coefficients
  for (int i = 1; i < n; i++) {
    double m = 1.0 / (b[i] - c[i-1] * a[i]);
    c[i] = c[i] * m;
    d[i] = (d[i] - d[i-1] * a[i]) * m;
  }
 
  // Back substitution
  x[n - 1] = d[n - 1];
  for (int i = n - 2; i >= 0; i--)
    x[i] = d[i] - c[i] * x[i+1];
}

void solveMatrix(
                int n,          // Number of equations
                double *a,      // Sub diagonal values
                double *b,      // Main diagonal values
                double *c,      // Super diagonal values
                double *v,      // Right hand side
                double *x)      // Answer
{
  // Forward sweep to modify coefficients for elimination of unknowns
  for (int i = 1; i < n; i++) {
    double m = a[i] / b[i - 1];
    b[i] -= (m * c[i - 1]);
    v[i] -= (m * v[i - 1]);
  }
 
  // Back substitution
  x[n - 1] = v[n - 1] / b[n - 1];
  for (int i = n - 2; i >= 0; i--)
    x[i] = (v[i] - (c[i] * x[i + 1])) / b[i];
}


///////////////////////////////////////////////////////////////////////
//
// Perform a filtering operation on the variable of interest
// Weighted average of a cell with the cells on either side with wrap around
//
///////////////////////////////////////////////////////////////////////

void smoothCellFilter(
                int nx,
                double* val)
{
  double* sVal = new double[nx];
  double onefourth = 1.0 / 4.0;

  // First cell
  sVal[0] = onefourth * (val[1] + (2.0 * val[0]) + val[nx - 1]);

  // Last cell
  sVal[nx - 1] = onefourth * (val[0] + (2.0 * val[nx - 1]) + val[nx - 2]);

  // Other cells
  for (int cIndx = 1; cIndx < (nx - 1); cIndx++)
    sVal[cIndx] = onefourth * 
                  (val[cIndx + 1] + (2.0 * val[cIndx]) + val[cIndx - 1]);

  for (int cIndx = 0; cIndx < nx; cIndx++)
    val[cIndx] = sVal[cIndx];

  delete [] sVal;
}

///////////////////////////////////////////////////////////////////////
//
// Perform a filtering operation on the variable of interest
// Weighted average of a face where first face and last face are the same
//
///////////////////////////////////////////////////////////////////////

void smoothFaceFilter(
                int nfx,
                double* val)
{
  double* sVal = new double[nfx];
  double onefourth = 1.0 / 4.0;

  // First face
  sVal[0] = onefourth * (val[1] + (2.0 * val[0]) + val[nfx - 2]);

  // Last face
  sVal[nfx - 1] = onefourth * (val[1] + (2.0 * val[nfx - 1]) + val[nfx - 2]);

  // Other faces
  for (int cIndx = 1; cIndx < (nfx - 1); cIndx++)
    sVal[cIndx] = onefourth * 
                  (val[cIndx + 1] + (2.0 * val[cIndx]) + val[cIndx - 1]);

  for (int cIndx = 0; cIndx < nfx; cIndx++)
    val[cIndx] = sVal[cIndx];

  delete [] sVal;
}

///////////////////////////////////////////////////////////////////////
//
// Create a face quantity for the periodic system
//
///////////////////////////////////////////////////////////////////////

void makeFaceFlux(
                double* cellVal,
                double* faceVal,
                int numCell)
{
  // First and last face average of first and last cells
  faceVal[0] = (cellVal[0] + cellVal[numCell - 1]) / 2.0;
  faceVal[numCell] = faceVal[0];

  // Faces in the middle
  for (int i = 1; i < numCell; i++) {
    faceVal[i] = (cellVal[i-1] + cellVal[i]) / 2.0;
  }
}

///////////////////////////////////////////////////////////////////////
//
// Linear solver Gaussian elimination with back substitution
//
///////////////////////////////////////////////////////////////////////

void GaussElim(
                double** A, 
                int n, 
                double * rhs, 
                int m, 
                double* x,
                double &residual, 
                int &flag)
{
  // declare our variables...
  double **ab;
  int i,j,k,l,np1,nm1,kp1,bigindex;
  double sumr, rowmax, scale, big, tempb, tempi, qout, sum;
  ab = new double* [n];

  // initialize our 2-D matrix to solve with...
  for (i=0; i<n; i++) ab[i] = new double [n+1];

  np1 = n+1;
  nm1 = n-1;
  // the columns of AB are [A, rhs]
  for (i=0; i<n; i++) {
    rowmax = 0.0;
    for (j=0; j<n; j++) {
      rowmax = rowmax > fabs( A[i][j] ) ? rowmax : fabs(A[i][j]);
    }
    // if the vector is NULL (zero-filled), then there is no solution
    if ( rowmax == 0.0 ) {
      flag = -1;
    }
    
    // scale and assign values for AB...
    scale = 1.0 / rowmax;
    for (j=0; j<n; j++) {
      ab[i][j] = A[i][j] * scale;
    }
    ab[i][n] = rhs[i] * scale;
  }


  // begin the elimination loop where rows of AB are interchanged to bring into 
  // alignment for pivotal positioning...
  for (k=0; k<n-1; k++) {
    big = 0.0;
    for (i=k; i<n; i++) {
      tempb = fabs( A[i][k] );
      if ( big >= tempb ) continue;
      big = tempb;
      bigindex = i;
    }

    // again if we have a null vector then bail out...
    if ( big == 0.0 ) {
      flag = -2;
    }

    // interchange rows (if we're not in the same row..)...
    if ( bigindex != k ) {
      for (i=k; i<n+1; i++) { 
        tempi = ab[k][i];
        ab[k][i] = ab[bigindex][i];
        ab[bigindex][i] = tempi;
      }
    }

    kp1 = k+1;
    
    // eliminate x(k) from equations 
    for (i=kp1; i<n; i++) {
      qout = ab[i][k] / ab[k][k];
      for (j=kp1; j<n+1; j++) {
        ab[i][j] = ab[i][j] - qout * ab[k][j];
      }
    }
  } // end of elimination loop

  //begin the calculation of the solution vector using back substitution...
  if ( ab[n-1][n-1] == 0.0 ) {
    // no solution 
    flag = -3;
  }

  x[n-1] = ab[n-1][n] / ab[n-1][n-1];
  for (i=1; i<=n; i++) {
    k = n - i;
    sum = 0.0;
    for (j=k+1; j<n; j++) {
      sum += ab[k][j] * x[j];
    }
    x[k] = ( ab[k][n] - sum ) / ab[k][k];
  }

  // calculate the maximum residual...
  residual = 0.0;
  for (i=0; i<n; i++) {
    sumr = 0.0;
    for (j=0; j<n; j++) {
      sumr += A[i][j] * x[j];
      residual = residual > fabs(sumr-rhs[i]) ? residual : fabs(sumr-rhs[i]);
      flag = 0;
    }
  }
  for (int i = 0; i < n; i++)
    delete [] ab[i];
  delete [] ab;
}

///////////////////////////////////////////////////////////////////////
//
// Perform Gauss-Jordan elimination with row-pivoting to obtain the solution to 
// the system of linear equations A x = B
// 
// Arguments:
//   lhs      - left-hand side of the equation, matrix A
//   rhs      - right-hand side of the equation, matrix B
//   nrows    - number of rows in the arrays lhs and rhs
//   ncolsrhs - number of columns in the array rhs
// 
// The function uses Gauss-Jordan elimination with pivoting.  The solution x to 
// the linear system winds up stored in the array rhs; create a copy to pass to
// the function if you wish to retain the original RHS array.
// 
// Passing the identity matrix as the rhs argument results in the inverse of 
// matrix A, if it exists.
// 
///////////////////////////////////////////////////////////////////////

void swaprows(double** arr, long row0, long row1)
{
  //  exchanges the contents of row0 and row1 in a 2d array
  double* temp;
  temp=arr[row0];
  arr[row0]=arr[row1];
  arr[row1]=temp;
}

void GaussJordanElim(
                double** lhs,
                double** rhs,
                long nrows,
                long ncolsrhs)
{
  // augment lhs array with rhs array and store in arr2
  double** arr2=new double*[nrows];
  for (long row=0; row<nrows; ++row)
    arr2[row]=new double[nrows+ncolsrhs];

  for (long row=0; row<nrows; ++row) {
    for (long col=0; col<nrows; ++col) {
      arr2[row][col]=lhs[row][col];
    }
    for (long col=nrows; col<nrows+ncolsrhs; ++col) {
      arr2[row][col]=rhs[row][col-nrows];
    }
  }

  // perform forward elimination to get arr2 in row-echelon form
  for (long dindex=0; dindex<nrows; ++dindex) {

    // run along diagonal, swapping rows to move zeros in working position 
    // (along the diagonal) downwards
    if ((dindex==(nrows-1)) && (arr2[dindex][dindex]==0)) {
      return; //  no solution
    } else if (arr2[dindex][dindex]==0) {
      swaprows(arr2, dindex, dindex+1);
    }

    // divide working row by value of working position to get 1 on the diagonal
    if (arr2[dindex][dindex] == 0.0) {
      return;
    } else {
      double tempval=arr2[dindex][dindex];
      for (long col=0; col<nrows+ncolsrhs; ++col) {
          arr2[dindex][col]/=tempval;
      }
    }

    // eliminate value below working position by subtracting a multiple of 
    // the current row
    for (long row=dindex+1; row<nrows; ++row) {
      double wval=arr2[row][dindex];
      for (long col=0; col<nrows+ncolsrhs; ++col) {
        arr2[row][col]-=wval*arr2[dindex][col];
      }
    }
  }

  // backward substitution steps
  for (long dindex=nrows-1; dindex>=0; --dindex) {

    // eliminate value above working position by subtracting a multiple of 
    // the current row
    for (long row=dindex-1; row>=0; --row) {
      double wval=arr2[row][dindex];
      for (long col=0; col<nrows+ncolsrhs; ++col) {
        arr2[row][col]-=wval*arr2[dindex][col];
      }
    }
  }

  // assign result to replace rhs
  for (long row=0; row<nrows; ++row) {
    for (long col=0; col<ncolsrhs; ++col) {
      rhs[row][col]=arr2[row][col+nrows];
    }
  }

  for (long row=0; row<nrows; ++row)
    delete[] arr2[row];
  delete[] arr2;
}

#define TINY 1.0e-20

///////////////////////////////////////////////////////////////////////
//
// Given a matrix a[0..n-1][0..n-1] this routine replaces it by the
// LU decomposition of a rowwise permutation of itself.  a and n are input.
// a is output, arranged for back substitution
// indx[0..n-1] is an output vector that records the row permutation
// effected by the partial pivoting.  d is output as +1/-1 depending on
// whether the number of row interchanges was even or odd.  This routine
// is used in combination with lubksb to solver linear equations or
// invert a matrix.
//
///////////////////////////////////////////////////////////////////////

void ludcmp(double **a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv = new double[n];   // store implicit scaling of each row

  *d=1.0;                       // no row interchanges yet

  // Loop over rows to get implicit scaling information
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++) {
      if ((temp=fabs(a[i][j])) > big)
        big=temp;
    }
    vv[i]=1.0/big;
  }

  // Loop over columns of Crout's method
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    // Initialize search for largest pivot element
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    // Do we need to interchange rows
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      // Change parity of number of interchanges, and exchange scale of rows
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0)
      a[j][j]=TINY;
    // Divide by the pivot element
    if (j != n-1) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++)
        a[i][j] *= dum;
    }
  }
  delete [] vv;
}
#undef TINY

///////////////////////////////////////////////////////////////////////
//
// Solves the set of n linear equations Ax = B.  a[0..n-1][0..n-1] is
// input, not as the matrix A but as its LU decomposition from ludcmp.
// indx[0..n-1] is input as the permutation vector.  b[0..n] is input
// as the right hand side vector B and returns with the solution vector X.
//
///////////////////////////////////////////////////////////////////////

void lubksb(double **a, int n, int *indx, double* b)
{
  int i,ii=0,ip,j;
  double sum;

  // Forward substitution
  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++)
        sum -= a[i][j]*b[j];
    else if (sum)
      ii=i;
    b[i]=sum;
  }

  // Back subtitution
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++)
      sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}
