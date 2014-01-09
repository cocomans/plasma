#include "AmpereSolver.h"
#include "../FieldDataCPU.h"
#include "../NodeHOMoments.h"


AmpereSolver::~AmpereSolver()
{

}

void AmpereSolver::init(PlasmaData* pdata,
		FieldDataCPU* _fields_old,
		NodeFieldData* fields_half,
		FieldDataCPU* _fields_next,
		NodeHOMoments* moments)
{
	fields_old = _fields_old;
	fields_next = _fields_next;
}

void AmpereSolver::solve(PlasmaData* pdata,FieldDataCPU* fields,NodeHOMoments* moments)
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
					jmean += moments->get_val(i,j,k,l,HOMoments_currentx,1)*pdata->qspecies[l];

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
					jtotal += moments->get_val(i,j,k,l,HOMoments_currentx,1)*pdata->qspecies[l];
				}

				fields_next->getE(i,j,k,0) = fields_old->getE(i,j,k,0) - pdata->dt
						* (jtotal - jmean)/pdata->xi;
				//fields->getE(i,j,k,1) -= qe*qe2me*pdata->dt/(pdata->nptcls*1.0*epsilon_naught)
				//		* moments->get_val(i,j,k,HOMoments_currenty)*0;
				//fields->getE(i,j,k,2) -= qe*qe2me*pdata->dt/(pdata->nptcls*1.0*epsilon_naught)
				//		* moments->get_val(i,j,k,HOMoments_currentz)*0;

			}
		}
	}

	fields->copy_from(fields_next);
}


realkind AmpereSolver::calc_residual(PlasmaData* pdata,
		FieldDataCPU* _fields_next,
		FieldDataCPU* _fields_old,
		NodeHOMoments* moments)
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
					jmean += moments->get_val(i,j,k,l,HOMoments_currentx,1)*pdata->qspecies[l];

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
					jtotal += moments->get_val(i,j,k,l,HOMoments_currentx,1)*pdata->qspecies[l];
				}

				Fx += (fabs((fields_next->getE(i,j,k,0) - fields_old->getE(i,j,k,0))/pdata->dt
										+  (jtotal-jmean)/pdata->xi));


			}
		}
	}


	double residual = Fx/((double)nx*ny*nz);











	return residual;



}
