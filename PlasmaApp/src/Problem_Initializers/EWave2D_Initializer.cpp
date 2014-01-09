#include "EWave2D_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../NodeHOMoments.h"
#include "../NodeParticleList.h"
#include "../FieldDataCPU.h"
#include "../PlasmaData.h"

#include "../ParallelInfo.h"
#include"../Util/rand.h"


#include "../RunData.h"

EWave2D_Initializer::EWave2D_Initializer(PlasmaData* pdata_in)
	{
		pdata = pdata_in;
		kt = pdata->Lx/(2.0*3.14159265359);
		alpha = 0.001;
		title = "TwoStream2D";

		Te0 = 0.000;

		pdata->rdata->SimName = title;

		np_y = 8*pdata->ny;
		np_x = (pdata->nptcls/2)/np_y;
		np_in = pdata->nptcls;

//		np_x = 8*pdata->nx;
//		np_y = (pdata->nptcls/2)/np_x;

		B0 = 0.3;

		lambda = 1.0/sqrt(pdata->mspecies[0]/pdata->mspecies[1]);

		cdf = new raised_cos(alpha,1.0/kt);

		px_dist[0] = (double*)malloc(np_in*sizeof(double));
		ix_dist[0] = (int*)malloc(np_in*sizeof(int));

		px_dist[1] = (double*)malloc(np_in*sizeof(double));
		ix_dist[1] = (int*)malloc(np_in*sizeof(int));

		py_dist[0] = (double*)malloc(np_in*sizeof(double));
		iy_dist[0] = (int*)malloc(np_in*sizeof(int));

		py_dist[1] = (double*)malloc(np_in*sizeof(double));
		iy_dist[1] = (int*)malloc(np_in*sizeof(int));

		vx_dist[0] = (double*)malloc(np_in*sizeof(double));
		vy_dist[0] = (double*)malloc(np_in*sizeof(double));
		vz_dist[0] = (double*)malloc(np_in*sizeof(double));

		vx_dist[1] = (double*)malloc(np_in*sizeof(double));
		vy_dist[1] = (double*)malloc(np_in*sizeof(double));
		vz_dist[1] = (double*)malloc(np_in*sizeof(double));




		for(int l=0;l<pdata->nspecies;l++)
			for(int i=0;i<np_in;i++)
			{


				init_velocities(vx_dist[l][i],vy_dist[l][i],vz_dist[l][i],l,i);

				int ixt = i%np_x;
				int ivel = (i)/(np_in/2);
				int iyt = (i-ixt)/np_x - np_y*ivel;

				int itemp = i%(np_x*np_y);

				double temp_x;
				double temp_y;
				if(l == 0)
					temp_x = distribution_intrp((ixt+0.5)/((double)np_x),0.0,pdata->Lx,*cdf);
				else
					temp_x = (ixt)/((double)np_x) * pdata->Lx + pdata->xmin;

				temp_y =  (iyt+0.5)/((double)np_y) * pdata->Ly + pdata->ymin;


//				double temp_x;
//				double temp_y;
//				if(l == 0)
//					temp_y = distribution_intrp((iyt+0.5)/((double)np_y),0.0,pdata->Ly,*cdf);
//				else
//					temp_y = (iyt)/((double)np_y) * pdata->Ly + pdata->ymin;
//
//				temp_x =  (ixt+0.5)/((double)np_x) * pdata->Lx + pdata->xmin;

				int ix = floor(temp_x*pdata->didx);
				int iy = floor(temp_y*pdata->didy);


				double px = temp_x*pdata->didx - ix;
				double py = temp_y*pdata->didy - iy;


				//px = drandom();

				ix = ((ix%pdata->nx)+pdata->nx)%pdata->nx;
				iy = ((iy%pdata->ny)+pdata->ny)%pdata->ny;

				px_dist[l][i] = px;
				ix_dist[l][i] = ix;
				py_dist[l][i] = py;
				iy_dist[l][i] = iy;


				if(ivel == 1)
				{
					vx_dist[l][i] = -vx_dist[l][itemp];
					vy_dist[l][i] = -vy_dist[l][itemp];
					vz_dist[l][i] = -vz_dist[l][itemp];
					px_dist[l][i] =  px_dist[l][itemp];
					ix_dist[l][i] =  ix_dist[l][itemp];
					py_dist[l][i] =  py_dist[l][itemp];
					iy_dist[l][i] =  iy_dist[l][itemp];

				}



			}



	}

void EWave2D_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
		 int& ix, int& iy, int& iz,
		 realkind& vx, realkind& vy, realkind& vz,int ispecies, int iptcl,
		 int nptcls,int ioffset)
{
	// Set Position Values, ifloat = 0-2
//	raised_cos cdf(alpha,1.0/kt);

	int itemp = iptcl%(np_in);


//	py = 0.5;
	pz = 0.0;

	iz = floor(pz*pdata->didz);

	pz = pz*pdata->didz - iz;


	px = px_dist[ispecies][itemp];
	py = py_dist[ispecies][itemp];
	ix = ix_dist[ispecies][itemp];
	iy = iy_dist[ispecies][itemp];
	vx = vx_dist[ispecies][itemp];
	vy = vy_dist[ispecies][itemp];
	vz = vz_dist[ispecies][itemp];

//	printf("ix,iy,iz[%i] = %i %i %i\n",iptcl,ix,iy,iz);


}

void EWave2D_Initializer::init_velocities(realkind& vx, realkind& vy, realkind& vz, int ispecies, int iptcl)
{
	randn(vx);
	randn(vy);
	randn(vz);

	vx*=Te0;
	vy*=Te0;
	vz*=Te0;

	if(ispecies == 0)
	{
		vx = 0.1;
		vy = 0.0;
		vz = 0.0;
	}




}

void EWave2D_Initializer::initialize_particles(NodeParticleList* particles,
										  NodeHOMoments* moments)
{


	particles -> init(this,moments);



}

void EWave2D_Initializer::initialize_fields(FieldDataCPU* fields,
										NodeHOMoments* moments)
{
	// Setup E-field
	for(int k=0;k<pdata->nz;k++)
	{
#pragma omp parallel for
		for(int j=0;j<pdata->ny;j++)
		{
			for(int i=0;i<pdata->nx;i++)
			{
				double xc = (i+0.5)*pdata->dxdi + pdata->xmin;
				double xf = (i)*pdata->dxdi + pdata->xmin;

				double yc = (j+0.5)*pdata->dydi + pdata->ymin;
				double yf = j*pdata->dydi + pdata->ymin;

				fields->getPhi(i,j,k) = -alpha*cos(xc/kt)*(kt*kt);
				fields->getE(i,j,k,0) = alpha*sin(xf/kt)*(kt);
				fields->getE(i,j,k,1) = 0;
				fields->getE(i,j,k,2) = 0;


				fields->getChi(i,j,k) = 0.0;
				fields->getA(i,j,k,0) = 0.0;
				fields->getA(i,j,k,1) = 0.0;
				fields->getA(i,j,k,2) = 0.0;
				fields->getB(i,j,k,0) = 0.0;
				fields->getB(i,j,k,1) = 0.0;
				fields->getB(i,j,k,2) = 0.0;
			}
		}
	}

}

void EWave2D_Initializer::check_step(NodeParticleList* particles_next,
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

void EWave2D_Initializer::finish(NodeParticleList* particles,
		NodeHOMoments* moments,NodeFieldData* fields_half)
{




}

