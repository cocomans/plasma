#include "Landau_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../HOMoments.h"
#include "../CurrentTally.h"
#include "../PlasmaData.h"
#include "../FieldData.h"
#include "../ParallelInfo.h"



Landau_Initializer::Landau_Initializer(PlasmaData* pdata_in)
{
		pdata = pdata_in;
		Ey = 2.0;
		Bz = 50;
		kt = 1.0*(pdata->Lx)/(2.0*3.14159265359);
		v_phase = 2.0*sqrt(2.0*pdata->Te/(pdata->mspecies[0]*mass_e));
		alpha = 0.01;
		printf("Phase Velocity = %f\n",v_phase);

}


void Landau_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
									 int& ix, int& iy, int& iz,
									 realkind& vx, realkind& vy, realkind& vz,int ispecies, int iptcl)
{
	// Set Position Values, ifloat = 0-2

	raised_cos cdf(alpha,1.0/kt);

	int nv1 = 100;
	int nv2 = 10;
	int nv3 = pdata->nptcls/(nv1*nv2);


	int ivel1 = iptcl/(nv2*nv3);
	int ivel2 = (iptcl - nv3*nv2*ivel1)/nv3;
	int idist = iptcl - nv3*(ivel2 + nv2*ivel1);

	if(ispecies == 0)
	{
		//px = distribution_intrp(idist/(1.0*nv3),0.0,pdata->Lx,cdf);

		px = distribution_intrp(iptcl/(1.0*pdata->nptcls),0.0,pdata->Lx,cdf);

	//px = box_muller(0.5,0.1);

	///px = ((2.0*iptcl)/((double)pdata->nptcls))*pdata->Lx+pdata->xmin;
	//vx = 1.0*(2*(ivel%2)-1)*(exp(-(ivel/2)/100.0));

	//vx = sqrt(2.0*pdata->Te/(pdata->mspecies[0]*mass_e))
	//	*sqrt(-2.0*log((ivel2+1.0)/((realkind)nv2)))*cos(2.0*3.14159265359*(ivel1)/((realkind)nv1));
	}
	else
	{
		px = drandom() * pdata->Lx + pdata->xmin;
		vx = box_muller(0.0,sqrt(2.0*pdata->Te/(pdata->mspecies[ispecies]*mass_e)));
	}

	vx = box_muller(0.0,sqrt(2.0*pdata->Te/(pdata->mspecies[ispecies]*mass_e)));
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
	vy = 0;

	vz = 0;





}


void Landau_Initializer::init_velocities(realkind& vx, realkind& vy, realkind& vz, int ispecies, int iptcl)
{
	// Set Velocity Values, ifloat = 3-5
	if(ispecies == 0)
		vx = 0.1*(1-(2*(rand()%2)));
	else
		vx = 0;
	vy = 0;

	vz = 0;



}

void Landau_Initializer::initialize_particles(ParticleList* particles,
										  HOMoments* moments,
										  ParallelInfo* myinfo)
{

	particles -> ispecies = myinfo->myid_mpi;

	particles -> init(this,moments);

	particles->plot_particles(pdata);



}

void Landau_Initializer::initialize_fields(FieldData** fields,
										HOMoments** moments,
										ParallelInfo* myinfo)
{

	int tid = myinfo->tid;

	// Setup E-field
	for(int k=0;k<pdata->nz;k++)
	{
		for(int j=0;j<pdata->ny;j++)
		{
			for(int i=0;i<pdata->nx;i++)
			{
				realkind x = (i)*pdata->dxdi+pdata->xmin;
				//realkind y = j*pdata->dydi+pdata->ymin;
				//realkind z = k*pdata->dzdi+pdata->zmin;

				realkind Ex = -alpha*sin(x/kt)*(kt);



				fields[tid+1] -> getE(i,j,k,0) = Ex/epsilon_naught;
				fields[tid+1] -> getE(i,j,k,1) = 0;
				fields[tid+1] -> getE(i,j,k,2) = 0;

				fields[tid+1] -> getB(i,j,k,0) = 0;
				fields[tid+1] -> getB(i,j,k,1) = 0;
				fields[tid+1] -> getB(i,j,k,2) = 0;

				fields[0] -> getE(i,j,k,0) = Ex/epsilon_naught;
				fields[0] -> getE(i,j,k,1) = 0;
				fields[0] -> getE(i,j,k,2) = 0;

				fields[0] -> getB(i,j,k,0) = 0;
				fields[0] -> getB(i,j,k,1) = 0;
				fields[0] -> getB(i,j,k,2) = 0;



			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

	if(myinfo->myid_mpi == 0)
	{
		realkind gamma;

		gamma = pi_const * powf(pdata->omega_pe,3.0f);

		plot = gnuplot_init();
		kenergy_plot = gnuplot_init();
		Tenergy_plot = gnuplot_init();
		Energy_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
		kEnergy_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
		TEnergy_array = (float*)malloc(2*pdata->nsteps*sizeof(float));

		max_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
		max_t_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
		nplot = 0;
		nmax_array = 0;

		kE0 = moments[0] -> evaluate_energy();
		kEnergy_array[0] = kE0;
		E0 = fields[0] -> evaluate_energy();
		Energy_array[0] = (E0);

		TEnergy_array[0] = E0+kE0;
		nplot++;

		gnuplot_cmd(plot,"set log y");
		//gnuplot_cmd(kenergy_plot,"set log y");


	}

}

void Landau_Initializer::check_step(ParticleList* particles,
		HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo)
{

	particles->plot_particles(pdata);

	// Print out the sum of the x-component of the current
	if(myinfo->myid_mpi == 0)
	{
		Energy_array[nplot] = ((fields[0] -> evaluate_energy()));
		kEnergy_array[nplot] = moments[0] -> evaluate_energy();

		TEnergy_array[nplot] = kEnergy_array[nplot] + Energy_array[nplot];
		printf("Energy Error= %e\n",(E0+kE0 - (Energy_array[nplot]+kEnergy_array[nplot]))/(E0+kE0));

		if(nplot == 1)
		{
			max_array[0] = log(Energy_array[0]);
			max_t_array[nmax_array] = pdata->dt;
			nmax_array++;
		}
		else if(nplot > 1)
		{
			if(Energy_array[nplot-1] > fmax(Energy_array[nplot-2],Energy_array[nplot]))
			{
				max_array[nmax_array] = log(Energy_array[nplot-1]);
				max_t_array[nmax_array] = (nplot-1);

				float period = (max_t_array[nmax_array]-max_t_array[nmax_array-1])*pdata->dt;
				float vphase = 2.0*pi_const/period* 1.0*pdata->Lx/(2.0*pi_const);
				printf("Energy fluctuation period = %f\n",period);
				printf("Phase velocity = %f\n",vphase);
				printf("Energy decay rate = %f\n",(max_array[nmax_array]-max_array[nmax_array-1])/
						(pdata->dt*(max_t_array[nmax_array]-max_t_array[nmax_array-1])));

				nmax_array++;
			}
		}

		float charge_cons = moments[0] -> check_charge(moments_old);

		printf("Charge conservation = %e\n",charge_cons);



		if(nplot > 0)
		{
			gnuplot_resetplot(plot);
			gnuplot_resetplot(kenergy_plot);
			gnuplot_resetplot(Tenergy_plot);
		}

		nplot++;
		gnuplot_plot_x(plot,Energy_array,nplot,"Field Energy");
		gnuplot_plot_x(kenergy_plot,kEnergy_array,nplot,"Particle Energy");
		gnuplot_plot_x(Tenergy_plot,TEnergy_array,nplot,"Total Energy");
		//gnuplot_setstyle(plot,"lines");
		//gnuplot_plot_xy(plot,max_t_array,max_array,nmax_array,NULL);
		//gnuplot_setstyle(plot,"points");


	}
















}

