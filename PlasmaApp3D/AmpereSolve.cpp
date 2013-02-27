/*-------------------------------------------------------------------------*/
/**
  @file		AmpereSolve.cpp
*/
/*--------------------------------------------------------------------------*/
#include "AmpereSolve.h"
#include "FieldDataCPU.h"
#include "HOMoments.h"
#include "PlasmaData.h"




void AmpereSolve(PlasmaData* pdata,
		FieldDataCPU* fields_next,
		FieldDataCPU* fields_old,
		HOMoments* moments)
{
	int nx = pdata->nx;
	int ny = pdata->ny;
	int nz = pdata->nz;

	for(int l=0;l<pdata->nspecies;l++)
	{
		for(int k=0;k<nz;k++)
		{
			for(int j=0;j<ny;j++)
			{
				for(int i=0;i<nx;i++)
				{
					fields_next->getE(i,j,k,0) = fields_old->getE(i,j,k,0)
							- pdata->n0*qe*qe2me*pdata->dt/(pdata->nptcls*1.0f)
							* moments->get_val(i,j,k,l,HOMoments_currentx);
					fields_next->getE(i,j,k,1) -= fields_old->getE(i,j,k,1)
							- pdata->n0*qe*qe2me*pdata->dt/(pdata->nptcls*1.0f)
							* moments->get_val(i,j,k,l,HOMoments_currenty);
					fields_next->getE(i,j,k,2) -= fields_old->getE(i,j,k,2)
							- pdata->n0*qe*qe2me*pdata->dt/(pdata->nptcls*1.0f)
							* moments->get_val(i,j,k,l,HOMoments_currentz);
				}
			}
		}
	}
}

double Calc_Residual(PlasmaData* pdata,
		FieldDataCPU* fields_next,
		FieldDataCPU* fields_old,
		HOMoments* moments)
{
	int nx = pdata->nx;
	int ny = pdata->ny;
	int nz = pdata->nz;

	double Fx,Fy,Fz;

	Fx = 0;
	Fy = 0;
	Fz = 0;

	for(int l=0;l<pdata->nspecies;l++)
	{
		for(int k=0;k<nz;k++)
		{
			for(int j=0;j<ny;j++)
			{
				for(int i=0;i<nx;i++)
				{
					Fx += ((fields_next->getE(i,j,k,0) - fields_old->getE(i,j,k,0))/pdata->dt
											+ pdata->n0*qe*qe2me/(pdata->nptcls*1.0f)
							* moments->get_val(i,j,k,l,HOMoments_currentx));
					Fy += ((fields_next->getE(i,j,k,1) - fields_old->getE(i,j,k,1))/pdata->dt
											+ pdata->n0*qe*qe2me/(pdata->nptcls*1.0f)
							* moments->get_val(i,j,k,l,HOMoments_currenty));
					Fz += ((fields_next->getE(i,j,k,2) - fields_old->getE(i,j,k,2))/pdata->dt
											+ pdata->n0*qe*qe2me/(pdata->nptcls*1.0f)
							* moments->get_val(i,j,k,l,HOMoments_currentz));

				}
			}
		}
	}

	Fx /= (double)nx*ny*nz;
	Fy /= (double)nx*ny*nz;
	Fz /= (double)nx*ny*nz;

	double residual = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);

	return residual;

}
