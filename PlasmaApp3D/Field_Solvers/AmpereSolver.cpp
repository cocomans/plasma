#include "AmpereSolver.h"
#include "../FieldData.h"
#include "../HOMoments.h"


AmpereSolver::~AmpereSolver()
{

}

void AmpereSolver::solve(PlasmaData* pdata,FieldData* fields,HOMoments* moments)
{
	int nx = pdata->nx;
	int ny = pdata->ny;
	int nz = pdata->nz;

	realkind jmean = 0.0;

	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int i=0;i<nx;i++)
			{
				for(int l=0;l<pdata->nspecies;l++)
				{
					jmean += moments->get_val(i,j,k,l,HOMoments_currentx)*pdata->qspecies[l];

				}

			}
		}
	}
//	printf("jmean = %f\n",jmean/(nx*ny*nz));

	jmean /= (double)(nx)*(ny)*(nz);




	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int i=0;i<nx;i++)
			{

				realkind jtotal = 0;



				for(int l=0;l<pdata->nspecies;l++)
				{
					jtotal += moments->get_val(i,j,k,l,HOMoments_currentx)*pdata->qspecies[l];
				}

				fields->getE(i,j,k,0) -= qe * pdata->dt/(epsilon_naught)
						* (jtotal - jmean);
				//fields->getE(i,j,k,1) -= qe*qe2me*pdata->dt/(pdata->nptcls*1.0*epsilon_naught)
				//		* moments->get_val(i,j,k,HOMoments_currenty)*0;
				//fields->getE(i,j,k,2) -= qe*qe2me*pdata->dt/(pdata->nptcls*1.0*epsilon_naught)
				//		* moments->get_val(i,j,k,HOMoments_currentz)*0;

			}
		}
	}
}


realkind AmpereSolver::calc_residual(PlasmaData* pdata,
		FieldData* fields_next,
		FieldData* fields_old,
		HOMoments* moments)
{

	int nx = pdata->nx;
	int ny = pdata->ny;
	int nz = pdata->nz;

	double jmean = 0.0;

	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int i=0;i<nx;i++)
			{
				for(int l=0;l<pdata->nspecies;l++)
				{
					jmean += moments->get_val(i,j,k,l,HOMoments_currentx)*pdata->qspecies[l];

				}

			}
		}
	}

	//printf("jmean = %f\n",jmean/(ny*nz));
	jmean /= (double)nx*ny*nz*1.0;




	double Fx,Fy,Fz;

	Fx = 0;
	Fy = 0;
	Fz = 0;

	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int i=0;i<nx;i++)
			{

					realkind jtotal = 0;



				for(int l=0;l<pdata->nspecies;l++)
				{
					jtotal += moments->get_val(i,j,k,l,HOMoments_currentx)*pdata->qspecies[l];
				}

				Fx += (fabs((fields_next->getE(i,j,k,0) - fields_old->getE(i,j,k,0))/pdata->dt
										+ qe/(epsilon_naught)
						* (jtotal-jmean)));

			}
		}
	}


	double residual = Fx/((double)nx*ny*nz);











	return residual;



}
