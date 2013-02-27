#include "Island_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../HOMoments.h"
#include "../CurrentTally.h"
#include "../PlasmaData.h"
#include "../FieldData.h"
#include "../ParallelInfo.h"

void Island_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
		 int& ix, int& iy, int& iz,
		 realkind& vx, realkind& vy, realkind& vz,int ispecies, int iptcl)
{
	// Set Position Values, ifloat = 0-2
//	raised_cos cdf(alpha,1.0/kt);

	int iptcl_new = floor((((double)iptcl*pdata->nptcls_cpu))/((double)pdata->my_nptcls));

		iptcl_new =	iptcl_new*pdata->num_nodes + (pdata->mynode);

	int idist = iptcl_new%(pdata->nptcls_cpu*pdata->num_nodes/2);
	int ivel = iptcl_new/(pdata->nptcls_cpu*pdata->num_nodes/2);
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

	iz = floor(pz*pdata->didz);

	pz = 0;

	//px = drandom();

	ix = ((ix%pdata->nx)+pdata->nx)%pdata->nx;
	iy = ((iy%pdata->ny)+pdata->ny)%pdata->ny;


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

void Island_Initializer::initialize_particles(ParticleList* particles,
										  HOMoments* moments,
										  ParallelInfo* myinfo)
{

	particles -> ispecies = pdata->my_species;

	//printf("my species = %i\n",particles->ispecies);

	particles -> init(this,moments);
	MPI_Barrier(MPI_COMM_WORLD);
//	kE0 = particles->evaluate_energy(pdata);

//	if(myinfo->myid_mpi == 0)
//		particles->plot_particles(pdata);

}

void Island_Initializer::initialize_fields(FieldData** fields,
										HOMoments** moments,
										ParallelInfo* myinfo)
{

	int tid = myinfo->tid;

	omp_set_num_threads(pdata->num_cores);
	// Setup E-field
	for(int k=0;k<pdata->nz;k++)
	{
#pragma omp for
		for(int j=0;j<pdata->ny;j++)
		{
			for(int i=0;i<pdata->nx;i++)
			{
				realkind x = i*pdata->dxdi+pdata->xmin;
				realkind y = j*pdata->dydi+pdata->ymin;
				//realkind z = k*pdata->dzdi+pdata->zmin;

				realkind Ex = -alpha*sin(x/kt)*(kt);



				fields[tid+1] -> getE(i,j,k,0) = 0;
				fields[tid+1] -> getE(i,j,k,1) = 0;
				fields[tid+1] -> getE(i,j,k,2) = 0;

				fields[tid+1] -> getB(i,j,k,0) = B0*sinh(y/lambda)/(cosh(y/lambda)+0.25*cos(x/lambda));
				fields[tid+1] -> getB(i,j,k,1) = B0*0.25*sin(x/lambda)/(cosh(y/lambda)+0.25*cos(x/lambda));
				fields[tid+1] -> getB(i,j,k,2) = 0;

				fields[0] -> getE(i,j,k,0) = 0;
				fields[0] -> getE(i,j,k,1) = 0;
				fields[0] -> getE(i,j,k,2) = 0;

				fields[0] -> getB(i,j,k,0) = B0*sinh(y/lambda)/(cosh(y/lambda)+0.25*cos(x/lambda));
				fields[0] -> getB(i,j,k,1) = B0*0.25*sin(x/lambda)/(cosh(y/lambda)+0.25*cos(x/lambda));
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

		if(pdata->plot_flag)
		{
			plot = gnuplot_init();
			kenergy_plot = gnuplot_init();
			Tenergy_plot = gnuplot_init();
		}
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

		TEnergy_array[0] =0;
		max_t_array[0] = 0;
		nplot++;

		if(pdata->plot_flag) gnuplot_cmd(plot,"set log y");
		//gnuplot_cmd(kenergy_plot,"set log y");

	}

}

void Island_Initializer::check_step(ParticleList* particles,
		HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo)
{


//	// Print out the sum of the x-component of the current
//
//	if(myinfo->myid_mpi == 0)
//	{
//		particles->plot_particles(pdata);
//
//
//		double pEnergy = ((fields[0] -> evaluate_energy()));
//		double kEnergy = moments[0] -> evaluate_energy();
//
//		Energy_array[nplot] = pEnergy;
//		kEnergy_array[nplot] = kEnergy;
//
//		TEnergy_array[nplot] = (E0+kE0 - (kEnergy + pEnergy))/(E0+kE0);
//		//printf("Energy Error= %e\n",(E0+kE0 - (Energy_array[nplot]+kEnergy_array[nplot]))/(E0+kE0));
//		max_t_array[nplot] = (nplot-1)*pdata->dt*2;
//		if(nplot == 1)
//		{
//			max_array[0] = log(Energy_array[0]);
//			max_t_array[nmax_array] = pdata->dt;
//			nmax_array++;
//		}
//		else if(nplot > 1)
//		{
//			if(Energy_array[nplot-1] > fmax(Energy_array[nplot-2],Energy_array[nplot]))
//			{
//				//max_array[nmax_array] = log(Energy_array[nplot-1]);
//				//max_t_array[nmax_array] = (nplot-1);
//
//				//realkind period = (max_t_array[nmax_array]-max_t_array[nmax_array-1])*pdata->dt;
//				//realkind vphase = 2.0*pi_const/period* 1.0*pdata->Lx/(2.0*pi_const);
//				//printf("Energy fluctuation period = %f\n",period);
//				//printf("Phase velocity = %f\n",vphase);
//				//printf("Energy decay rate = %f\n",(max_array[nmax_array]-max_array[nmax_array-1])/
//				//		(pdata->dt*(max_t_array[nmax_array]-max_t_array[nmax_array-1])));
//
//				nmax_array++;
//			}
//		}
//
//		realkind charge_cons = moments[0] -> check_charge(moments_old);
//
//		//printf("Charge conservation = %e\n",charge_cons);
//
//
//
//		if((nplot > 0)&&(nplot%10 == 9))
//		{
//			if(pdata->plot_flag){
//			gnuplot_resetplot(plot);
//			gnuplot_resetplot(kenergy_plot);
//			gnuplot_resetplot(Tenergy_plot);
//			}
//		}
//
//		nplot++;
//
//
//		if((pdata->plot_flag)&&(nplot%10 == 9))
//		{
//			gnuplot_plot_xy(plot,max_t_array+1,Energy_array+1,nplot-1,"Field Energy");
//			gnuplot_cmd(plot,"replot \"E_vec.csv\" title \"Matlab Code\" with points");
//			gnuplot_plot_x(kenergy_plot,kEnergy_array,nplot,"Particle Energy");
//			gnuplot_plot_x(Tenergy_plot,TEnergy_array,nplot,"Total Energy");
//		}
//
//
//	}

}

void Island_Initializer::finish(ParticleList* particles,
		HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo)
{


	// Print out the sum of the x-component of the current
	if(myinfo->myid_mpi == 0)
	{
		if(pdata->plot_flag)
		{
			gnuplot_resetplot(plot);
			gnuplot_cmd(plot,"set xlabel \"time\"");
			gnuplot_cmd(plot,"set ylabel \"Field Energy\"");
			gnuplot_plot_xy(plot,max_t_array,Energy_array,nplot,"PlasmaApp3D");
			gnuplot_cmd(plot,"replot \"E_vec.csv\" title \"Matlab Code\" with points");
			//getchar();
			gnuplot_save_pdf(plot,"TS_field_energy");
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);

}

