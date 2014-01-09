#include "InitSolve.h"



void solve_tridiagonal_in_place_destructive(double x[],
		const size_t N, const double a[], const double b[], double c[]) {
    /* unsigned integer of same size as pointer */
    size_t in;


    /*
     solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
     note that contents of input vector c will be modified, making this a one-time-use function
     x[] - initially contains the input vector v, and returns the solution x. indexed from [0, ..., N - 1]
     N — number of equations
     a[] - subdiagonal (means it is the diagonal below the main diagonal) -- indexed from [1, ..., N - 1]
     b[] - the main diagonal, indexed from [0, ..., N - 1]
     c[] - superdiagonal (means it is the diagonal above the main diagonal) -- indexed from [0, ..., N - 2]
     */

    c[0] = c[0] / b[0];
    x[0] = x[0] / b[0];

    /* loop from 1 to N - 1 inclusive */
    for (in = 1; in < N; in++) {
        double m = 1.0 / (b[in] - a[in] * c[in - 1]);
        c[in] = c[in] * m;
        x[in] = (x[in] - a[in] * x[in - 1]) * m;
    }

    /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
    for (in = N - 1; in-- > 0; )
        x[in] = x[in] - c[in] * x[in + 1];
}



void solve_poisson_periodic(double x[],double dx,const int N)
{
	double* x1 = (double*)malloc(N*sizeof(double));
	double* x2 = (double*)malloc(N*sizeof(double));

	double* a = (double*)malloc(N*sizeof(double));
	double* b = (double*)malloc(N*sizeof(double));
	double* c = (double*)malloc(N*sizeof(double));

//	for(int i=0;i<N;i++)
//	{
//		a[i] = 1;
//		b[i] = -2;
//		c[i] = 1;
//
//		x1[i] = x[i];
//		x2[i] = 0;
//	}
//
//	for(int i=0;i<N;i++)
//		x[i] *= -dx*dx;
//
//	for(int i=0;i<N;i++)
//		x2[0] += i*x[i];
//
//	x2[0] = x2[0]/((double)N);
//	x2[1] = x[0] + 2.0*x2[0];
//
//	for(int i=2;i<N;i++)
//		x2[i] = x[i-1] + 2.0*x2[i-1] - x2[i-2];
//
//	double mean = 0.0;
//	for(int i=0;i<N;i++)
//		mean -= x2[i];
//
//	for(int i=0;i<N;i++)
//		x[i] = x2[i] + mean/((double)N);

// #pragma omp parallel for
// 	for(int i=0;i<N;i++)
// 	{
// 		b[i] = 1/(dx*dx);
// 		a[i] = -2/(dx*dx);
// 		c[i] = 1/(dx*dx);

// 		x1[i] = x[i];
// 		x2[i] = 0;
// 	}
// 	x2[0] = -a[0];
// 	x2[N-2] = -c[N-2];

// 	solve_tridiagonal_in_place_destructive(x1+1,N-1,b,a,c);
// 	solve_tridiagonal_in_place_destructive(x2+1,N-1,b,a,c);

// 	double xn = (x[N-1] - c[N-1]*x1[0] - b[N-1]*x1[N-2])/(a[N-1] + c[N-1]*x2[0] + b[N-1]*x2[N-2]);

// #pragma omp parallel for
// 	for(int i=0;i<N-1;i++)
// 	{
// 		x[i] = x1[i] + x2[i]*xn;
// 	}

// 	x[N-1] = xn;


// Thomas algorithm 
	
	for(int i=0;i<N;i++){
	  //rhs
	  x1[i] = -x[i]*dx*dx;
	  x[i]  = 0;
	}


	for(int i=0;i<N;i++){
	  x[0] = x[0] + (i+1)*x1[i];
	}
	x[0] = x[0]/((double)(N));
	x[1] = x1[0] + 2.0*x[0];
	for(int i=2;i<N;i++){
	  x[i] = x1[i-1] + 2.0*x[i-1] - x[i-2];
	}
	
	//ensure average zero
	double sum=0.0;
	for(int i=0;i<N;i++){
	  sum += x[i];
	}
	sum /= double(N);
	for(int i=0;i<N;i++){
	  x[i] = x[i] - sum;
	}


	free(x1);
	free(x2);
	free(a);
	free(b);
	free(c);

}

void InitSolve(PlasmaData* pdata,FieldDataCPU* fields,HOMomentsCPU* moments)
{
  //to check the initial conditions
  FILE * ftemp = fopen("Ini.txt","w");
  if (ftemp == NULL)
    {
      printf("Error opening file 'Ini.txt'!\n");
      exit(1);
    }

	double phi[pdata->nx];
	double Ay[pdata->nx];
	double Az[pdata->nx];

	double jmeany = 0;
	double jmeanz = 0;

	for(int i=0;i<pdata->nx;i++)
	{
		double current = 0;
		for(int l=0;l<pdata->nspecies;l++)
		{
			current += pdata->qspecies[l]*moments->get_val(i,0,0,l,HOMoments_currenty);
		}

		jmeany += current;

		current = 0;
		for(int l=0;l<pdata->nspecies;l++)
		{
			current += pdata->qspecies[l]*moments->get_val(i,0,0,l,HOMoments_currentz);
		}

		jmeanz += current;
	}

	jmeany /= (double)pdata->nx;
	jmeanz /= (double)pdata->nx;


	// Solve for Ex
#pragma omp parallel for
	for(int i=0;i<pdata->nx;i++)
	{
		double charge = 0;
		for(int l=0;l<pdata->nspecies;l++)
		{
			charge += pdata->qspecies[l]*moments->get_val(i,0,0,l,HOMoments_charge);
		}

		double currenty = 0;
		double currentz = 0;
		for(int l=0;l<pdata->nspecies;l++)
		{
			currenty += pdata->qspecies[l]*moments->get_val(i,0,0,l,HOMoments_currenty);
			currentz += pdata->qspecies[l]*moments->get_val(i,0,0,l,HOMoments_currentz);

		}


		phi[i] = charge/(pdata->xi*pdata->xi);

		Ay[i] = (currenty - jmeany)/pdata->bobzeta;
		Az[i] = (currentz - jmeanz)/pdata->bobzeta;



	}

	//print out the charge and current density
	for(int i=0;i<pdata->nx;i++)
	  {
	    fprintf(ftemp, " %e %e %e\n", phi[i],Ay[i],Az[i]);
	  }

	// Solve the system
	solve_poisson_periodic(phi,pdata->dxdi,pdata->nx);
	solve_poisson_periodic(Ay,pdata->dxdi,pdata->nx);
	solve_poisson_periodic(Az,pdata->dxdi,pdata->nx);

	double phi_avg = 0;
	double ay_avg = 0;
	double az_avg = 0;

	for(int i=0;i<pdata->nx;i++)
	{
		phi_avg += phi[i];
		ay_avg += Ay[i];
		az_avg += Az[i];
	}

	phi_avg /= (double)pdata->nx;
	ay_avg /= (double)pdata->nx;
	az_avg /= (double)pdata->nx;
#pragma omp parallel for
	for(int i=0;i<pdata->nx;i++)
	{
		phi[i] -= phi_avg;
		Ay[i] -= ay_avg;
		Az[i] -= az_avg;
	}

	// Convert Phi to Ex
#pragma omp parallel for
	for(int i=0;i<pdata->nx;i++)
	{
		int im = (((i-1)%(pdata->nx))+pdata->nx)%pdata->nx;
		fields->getE(i,0,0,0) = (phi[i] - phi[im])/pdata->dxdi;

//		fields->getE(i,0,0,1) = 2.0*Ay[i]/pdata->dt;
//		fields->getE(i,0,0,2) = 2.0*Az[i]/pdata->dt;


		fields->getA(i,0,0,1) = Ay[i];
		fields->getA(i,0,0,2) = Az[i];

		fields->getB(i,0,0,1) = -(Az[i] - Az[im])/pdata->dxdi;
		fields->getB(i,0,0,2) = (Ay[i] - Ay[im])/pdata->dxdi;
	}

	//close the ini.txt
	fclose(ftemp);
}


























