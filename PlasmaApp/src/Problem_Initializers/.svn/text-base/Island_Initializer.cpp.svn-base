#include "Island_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../NodeHOMoments.h"
#include "../NodeParticleList.h"
#include "../FieldDataCPU.h"
#include "../PlasmaData.h"

#include "../ParallelInfo.h"

#include "../RunData.h"

Island_Initializer::Island_Initializer(PlasmaData* pdata_in)
	{
		pdata = pdata_in;
		kt = pdata->Lx/(2.0*3.14159265359);
		alpha = 0.001;
		title = "Island Coalescence";

		pdata->rdata->SimName = title;

		B0 = 0.3;

		lambda = 1.0/sqrt(pdata->mspecies[0]/pdata->mspecies[1]);

		cdf = new island_coalescence(0.25,
				lambda,
				pdata->Lx,
				pdata->Ly,
				0.0,
				pdata->ymin,
				pdata->nx,
				pdata->ny);



	}

void Island_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
		 int& ix, int& iy, int& iz,
		 realkind& vx, realkind& vy, realkind& vz,int ispecies, int iptcl,
		 int nptcls,int ioffset)
{
	// Set Position Values, ifloat = 0-2
//	raised_cos cdf(alpha,1.0/kt);


	double temp_x;


	temp_x = distribution_intrp(drandom(),
			0.0,pdata->nx*pdata->ny,*cdf);

	int itemp = floor(temp_x);
	iy = itemp/pdata->nx ;
	ix = itemp - iy*pdata->nx -1;

	px = drandom();
	py = drandom();
//	int iptcl_new = iptcl;
//

//	int idist = iptcl_new%(pdata->nptcls/2);
//	int ivel = iptcl_new/(pdata->nptcls/2);
//
//	double temp_x;
//
//	if(ispecies == 0)
//		temp_x = distribution_intrp((idist+0.5)/(0.5*pdata->nptcls),0.0,pdata->Lx,cdf);
//	else
//		px = iptcl/(1.0*pdata->nptcls) * pdata->Lx + pdata->xmin;
	//px = box_muller(0.5,0.1);

	///px = ((2.0*iptcl)/((double)pdata->nptcls))*pdata->Lx+pdata->xmin;
	//vx = 1.0*(2*(ivel%2)-1)*(exp(-(ivel/2)/100.0));

	//vx = sqrt(2.0*pdata->Te/(pdata->mspecies[0]*mass_e))
	//	*sqrt(-2.0*log((ivel2+0.5)/((realkind)nv2)))*cos(2.0*3.14159265359*(ivel1)/((realkind)nv1));

	//temp_x = drandom();

	//printf("random number = %f\n",px);

	// Set Velocity Values, ifloat = 3-5

	vx = box_muller(0.0,1.0/sqrt(pdata->mspecies[ispecies]));
	vy = box_muller(0.0,1.0/sqrt(pdata->mspecies[ispecies]));
	vz = B0/(4*pi_const*lambda)*pdata->qspecies[ispecies];


	pz = 0.0;

	iz = 0;

	pz = 0;

	//px = drandom();

	ix = ((ix%pdata->nx)+pdata->nx)%pdata->nx;
	iy = ((iy%pdata->ny)+pdata->ny)%pdata->ny;

//	printf("ix,iy,iz[%i] = %i %i %i\n",iptcl,ix,iy,iz);


}

void Island_Initializer::init_velocities(realkind& vx, realkind& vy, realkind& vz, int ispecies, int iptcl)
{
	// Set Velocity Values, ifloat = 3-5
	if(ispecies == 0)
		vx = 0.1*(1-(2*(rand()%2)));
	else
		vx = 0;
	vy = 0;

	vz = 0;



}

void Island_Initializer::initialize_particles(NodeParticleList* particles,
										  NodeHOMoments* moments)
{


	particles -> init(this,moments);



}

void Island_Initializer::initialize_fields(FieldDataCPU* fields,
										NodeHOMoments* moments)
{
	// Setup E-field
	for(int k=0;k<pdata->nz+1;k++)
	{
#pragma omp parallel for
		for(int j=0;j<pdata->ny+1;j++)
		{
			for(int i=0;i<pdata->nx+1;i++)
			{
				realkind x2 = i*pdata->dxdi+pdata->xmin;
				realkind y2 = (j+0.5)*pdata->dydi-0.5*pdata->Ly;

				realkind x1 = (i+0.5)*pdata->dxdi+pdata->xmin;
				realkind y1 = j*pdata->dydi-0.5*pdata->Ly;
				//realkind z = k*pdata->dzdi+pdata->zmin;




				fields -> getE(i,j,k,0) = 0;
				fields -> getE(i,j,k,1) = 0;
				fields -> getE(i,j,k,2) = 0;

				fields -> getB(i,j,k,0) = -B0*sinh(y1/(lambda))/(cosh(y1/(lambda))+0.25*cos(x2/lambda));
				fields -> getB(i,j,k,1) = -B0*0.25*sin(x2/lambda)/(cosh(y1/(lambda))+0.25*cos(x2/lambda));
				fields -> getB(i,j,k,2) = 0;

				fields->getA(i,j,k,0) = 0.0;
				fields->getA(i,j,k,1) = lambda*B0*log(cosh(y1/lambda)+0.25*cos(x1/lambda));
				fields->getA(i,j,k,2) = 0.0;

				fields->getPhi(i,j,k) = 0.0;
				fields->getChi(i,j,k) = 0.0;


			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

}

void Island_Initializer::check_step(NodeParticleList* particles_next,
		NodeHOMoments* moments,FieldDataCPU* fields_old,FieldDataCPU* fields_next)
{
	if(pdata->node_info->rank_g == 0)
	{
		double errors = 0;
		double denom = 0;

		for(int m=0;m<4;m++)

		for(int j=0;j<pdata->ny;j++)
			for(int i=0;i<pdata->nx;i++)
			{
				enum HOMoments_moment mom = moments->moments_next->get_enum(m);
				double temp0 = 0;
				double temp1 = 0;
				for(int l=0;l<pdata->nspecies;l++)
				{
					temp0 += moments->get_val(i,j,0,l,mom,0)*pdata->qspecies[l];
					temp1 += moments->get_val(i,j,0,l,mom,1)*pdata->qspecies[l];
				}
				errors+= fabs(temp1-temp0)/(fabs(temp0)+fabs(temp1));
				denom++;
			}

		printf("Change in moments = %e\n",errors/(denom));
	}


}

void Island_Initializer::finish(NodeParticleList* particles,
		NodeHOMoments* moments,NodeFieldData* fields_half)
{




}

