#include "IonAcoustic_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../NodeHOMoments.h"
#include "../NodeParticleList.h"
#include "../FieldDataCPU.h"
#include "../PlasmaData.h"

#include "../ParallelInfo.h"


IonAcoustic_Initializer::IonAcoustic_Initializer(PlasmaData* pdata_in)
{
		pdata = pdata_in;
		Ey = 2.0;
		Bz = 50;
		kt = 1.0*(pdata->Lx)/(2.0*pi_const);
		v_shift = sqrt(1.0/pdata->mspecies[1]*1.0*2.0e-4/(1.0+1.0/(kt*kt)) + 1.0*2.0e-4/pdata->mspecies[1]);
		v_shift = -2.23e-2 + (93.5-72)/3000.0;
		alpha = 0.4;

//		for(int i=0;i<pdata->nx;i++)
//		{
//			realkind ux = alpha*v_shift*sin((0.5+i)*pdata->dxdi/kt) - v_shift;
//			printf("u[%i] = %f\n",i,ux);
//		}


}


void IonAcoustic_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
									 int& ix, int& iy, int& iz,
									 realkind& vx, realkind& vy, realkind& vz,
									 int ispecies, int iptcl,int nptcls,int ioffset)
{
	// Set Position Values, ifloat = 0-2
	raised_cos cdf(alpha,1.0/kt);
	int nv1 = 4;

	int iptcl_new = floor((((double)iptcl*pdata->nptcls_device[0]))/((double)nptcls));

	iptcl_new =	iptcl_new*pdata->node_info->nTasks_g + (pdata->node_info->rank_g);

	iptcl_new = iptcl;

	int idist = iptcl_new%(nptcls/nv1);
	double temp_x;

	px = distribution_intrp((idist)/(1.0/nv1*nptcls),0.0,pdata->Lx,cdf)+0.5*pdata->Lx;
//	px = distribution_intrp(drandom(),0.0,pdata->Lx,cdf);

	//vx = 0.0;
//	vx = (2*(iptcl%2)-1)*sqrt(2.0*pdata->Te/(pdata->mspecies[0]*mass_e));
	py = 0.5;
	pz = 0.5;



	ix = floor(px*pdata->didx);
	iy = floor(py*pdata->didy);
	iz = floor(pz*pdata->didz);

	px = px*pdata->didx - ix;
	py = py*pdata->didy - iy;
	pz = pz*pdata->didz - iz;

	ix = ((ix%pdata->nx)+pdata->nx)%pdata->nx;

	// Set Velocity Values, ifloat = 3-5


//	realkind ux = alpha*v_shift*sin((px+ix)*pdata->dxdi/kt) + v_shift;
	realkind ux = alpha*sin((px+ix)*pdata->dxdi/kt)*9.7e-4*kt + v_shift;


	if(ispecies == 0)
		vx = box_muller(ux,sqrt(2.0*2.0e-4/(pdata->mspecies[ispecies]*mass_e)));
	else
		vx = box_muller(ux,sqrt(2.0*2.0e-4/(pdata->mspecies[ispecies]*mass_e)));


	vy = 0;

	vz = 0;





}


void IonAcoustic_Initializer::init_velocities(realkind& vx, realkind& vy, realkind& vz, int ispecies, int iptcl)
{
	// Set Velocity Values, ifloat = 3-5
	if(ispecies == 0)
		vx = 0.1*(1-(2*(rand()%2)));
	else
		vx = 0;
	vy = 0;

	vz = 0;



}

void IonAcoustic_Initializer::initialize_particles(NodeParticleList* particles,
										  NodeHOMoments* moments)
{

	particles -> init(this,moments);

}

void IonAcoustic_Initializer::initialize_fields(FieldDataCPU* fields, NodeHOMoments* moments)
{


	// Setup E-field
	for(int k=0;k<pdata->nz;k++)
	{
		for(int j=0;j<pdata->ny;j++)
		{
			for(int i=0;i<pdata->nx;i++)
			{
				realkind x = i*pdata->dxdi+pdata->xmin;
				//realkind y = j*pdata->dydi+pdata->ymin;
				//realkind z = k*pdata->dzdi+pdata->zmin;

				realkind Ex = -0.0*alpha*sin(x/kt)*(kt);



				fields -> getE(i,j,k,0) = Ex/epsilon_naught;
				fields -> getE(i,j,k,1) = 0;
				fields -> getE(i,j,k,2) = 0;

				fields -> getB(i,j,k,0) = 0;
				fields -> getB(i,j,k,1) = 0;
				fields -> getB(i,j,k,2) = 0;



			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}




}

