#include "ExB_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../HOMoments.h"
#include "../CurrentTally.h"
#include "../FieldData.h"
#include "../ParallelInfo.h"


ExB_Initializer::ExB_Initializer(PlasmaData* pdata_in)
	{
		pdata = pdata_in;
		Ex0 = -0.2;
		Ex1 = 0.2;
		Bz = 1;

		Ex0 = 0;
		Ex1 = 0;
//		Ex0 = 0.0;
//		Ex1 = 0.0;

		dEdx = (Ex1-Ex0)/pdata->Lx;

		particles_temp = new ParticleListCPU();
		particles_temp -> allocate(pdata,pdata->my_nptcls);
		particles_temp -> ispecies = pdata->my_species;

		int ispecies = pdata->my_species;
		realkind q2m = pdata->qspecies[ispecies]/pdata->mspecies[ispecies]*-1.0;
		omegac = Bz*q2m;

		Omega = -sqrt(fabs(omegac*omegac - q2m*dEdx));

		nsteps = pdata->nsteps+1;

		vx0 = (realkind*)malloc(pdata->my_nptcls*sizeof(realkind));
		vy0 = (realkind*)malloc(pdata->my_nptcls*sizeof(realkind));

		x0 = (realkind*)malloc(pdata->my_nptcls*sizeof(realkind));
		y0 = (realkind*)malloc(pdata->my_nptcls*sizeof(realkind));

		xt = (float*)malloc(pdata->my_nptcls*nsteps*sizeof(float));
		yt = (float*)malloc(pdata->my_nptcls*nsteps*sizeof(float));
		xrt = (float*)malloc(pdata->my_nptcls*nsteps*sizeof(float));
		yrt = (float*)malloc(pdata->my_nptcls*nsteps*sizeof(float));


		curtime = 0;

		iter = 0;


	}

void ExB_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
		 int& ix, int& iy, int& iz,
		 realkind& vx, realkind& vy, realkind& vz,int ispecies, int iptcl)
{
	// Set Position Values, ifloat = 0-2
	px = (rand()%10000/10000.0);

	if(pdata->ndimensions > 1)
		py = (rand()%10000/10000.0);
	else
		py = 0;

	pz = 0.0;

	vx = 0.4*(((rand()%10000))/10000.0 -0.5);
	vy = 0.4*(((rand()%10000))/10000.0 -0.5);
	vz = 0;
	// Set Position Index Values, iint = 0-2
	ix = pdata->nx/2;
	if(pdata->ndimensions > 1)
		iy = pdata->ny/2;
	else
		iy = 0;

	iz = 0;
}

void ExB_Initializer::init_velocities(realkind& vx, realkind& vy, realkind& vz, int ispecies, int iptcl)
{
	// Set Velocity Values, ifloat = 3-5
	vx = 1.5*(2*(rand()%10000))/10000.0 + 0.5;
	vy = 1.5*(2*(rand()%10000))/10000.0 + 0.5;
	vz = 0.0;
}

void ExB_Initializer::initialize_particles(ParticleList* particles,
										  HOMoments* moments,
										  ParallelInfo* myinfo)
{

	particles -> ispecies = 0;

	particles -> init(this,moments);

	particles_temp -> copy_from(particles);

//	particles->plot_particles(pdata);

	for(int i=0;i<pdata->my_nptcls;i++)
	{
		realkind x,y,xl,yl,vx,vy;
		int ix,iy;

		xl = particles_temp->get_fvalue(i,0);
		yl = particles_temp->get_fvalue(i,1);

		vx = particles_temp->get_fvalue(i,3);
		vy = particles_temp->get_fvalue(i,4);

		ix = particles_temp->get_ivalue(i,0);
		iy = particles_temp->get_ivalue(i,1);

		x = pdata->dxdi * (xl+ix) + pdata->xmin;
		y = pdata->dydi * (yl+iy) + pdata->ymin;

		if(pdata->ndimensions == 1)
		{
			y = curtime;
		}

		vx0[i] = vx;
		vy0[i] = vy;
		x0[i] = x;
		y0[i] = y;

		xt[nsteps*i] = x;
		yt[nsteps*i] = y;
		xrt[nsteps*i] = x;
		yrt[nsteps*i] = y;
	}


}

void ExB_Initializer::initialize_fields(FieldData** fields,
										HOMoments** moments,
										ParallelInfo* myinfo)
{

	int tid = myinfo->tid;

	// Setup E-field
	for(int k=0;k<pdata->nz+1;k++)
	{
		for(int j=0;j<pdata->ny+1;j++)
		{
			for(int i=0;i<pdata->nx+1;i++)
			{
				int it = ((i%pdata->nx)+pdata->nx)%pdata->nx;
				realkind x = it*pdata->dxdi+pdata->xmin;
				//realkind y = j*pdata->dydi+pdata->ymin;
				//realkind z = k*pdata->dzdi+pdata->zmin;

				realkind Ex = dEdx*(x-pdata->xmin)+Ex0;

				fields[tid+1] -> getE(i,j,k,0) = Ex;
				fields[tid+1] -> getE(i,j,k,1) = 0;
				fields[tid+1] -> getE(i,j,k,2) = 0;

				fields[tid+1] -> getB(i,j,k,0) = 0;
				fields[tid+1] -> getB(i,j,k,1) = 0;
				fields[tid+1] -> getB(i,j,k,2) = Bz;



				fields[0] -> getE(i,j,k,0) = Ex;
				fields[0] -> getE(i,j,k,1) = 0;
				fields[0] -> getE(i,j,k,2) = 0;

				fields[0] -> getB(i,j,k,0) = 0;
				fields[0] -> getB(i,j,k,1) = 0;
				fields[0] -> getB(i,j,k,2) = Bz;


			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

}

void ExB_Initializer::check_step(ParticleList* particles,
		HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo)
{



	iter += 1;
	curtime += pdata->dt;
	double total_error = 0;
	realkind t_error;
	realkind vx_error,vy_error;
	realkind x_error,y_error;



	// Copy in the pushed particles
	particles_temp -> copy_from(particles);

	//particles_temp->plot_particles(pdata);


	// iterate over particles, compare position to analytic solution
	for(int i=0;i<pdata->my_nptcls;i++)
	{
		realkind x,y,xl,yl,vx,vy;
		int ix,iy;
		realkind xr,yr,vxr,vyr;

		xl = particles_temp->get_fvalue(i,0);
		yl = particles_temp->get_fvalue(i,1);

		vx = particles_temp->get_fvalue(i,3);
		vy = particles_temp->get_fvalue(i,4);

		ix = particles_temp->get_ivalue(i,0);
		iy = particles_temp->get_ivalue(i,1);

		x = pdata->dxdi * (xl+ix) + pdata->xmin;
		y = pdata->dydi * (yl+iy) + pdata->ymin;



		// Get the analytic solution
		get_analytic_pos(xr,yr,vxr,vyr,i);

		if(pdata->ndimensions == 1)
		{
			y = curtime;
			yr = curtime;
		}

		xt[nsteps*i+iter] = x;
		yt[nsteps*i+iter] = y;
		xrt[nsteps*i+iter] = xr;
		yrt[nsteps*i+iter] = yr;

		// Compare the analytic and pusher solutions
		if(pdata->ndimensions == 1)
		{	// X position only
			y_error = 0;
			x_error = fabs(x-xr)/fabs(xr);
			vx_error = fabs(vx-vxr)/fabs(vxr);
			vy_error = fabs(vy-vyr)/fabs(vyr);

			t_error = sqrt(x_error*x_error+vx_error*vx_error+vy_error*vy_error);
		}
		else
		{
			x_error = fabs(x-xr)/fabs(xr);
			y_error = fabs(y-yr);
			vx_error = fabs(vx-vxr)/fabs(vxr);
			vy_error = fabs(vy-vyr)/fabs(vyr);

			t_error = sqrt(x_error*x_error+y_error*y_error
					+vx_error*vx_error+vy_error*vy_error);
		}

		if(t_error > 1.0e-6)
			printf("Errors[%i] = %e, %e, %e, %e\n",i,x_error,y_error,vx_error,vy_error);

		total_error += t_error;


	}

	printf("Average Error for step %i = %e\n",iter,total_error/pdata->my_nptcls);





}

void ExB_Initializer::get_analytic_pos(realkind &xout,realkind &yout,
									   realkind &vxout,realkind &vyout,
									   int iptcl)
{
	int ispecies = pdata->my_species;
	realkind q2m = pdata->qspecies[ispecies]/pdata->mspecies[ispecies]*-1.0;

	realkind dadx = omegac*omegac - q2m * dEdx;
	realkind a0 = omegac*(vy0[iptcl] + Ex0/Bz + omegac*x0[iptcl]);
	realkind x00;
	realkind sinodt,cosodt, A, B;
	x00 = a0/dadx;

	A = (x0[iptcl] - x00);
	B = vx0[iptcl]/Omega;

	if(dadx == 0)
	{
		xout = x0[iptcl] + vx0[iptcl]*curtime + 0.5*a0*curtime*curtime;
		vxout = vx0[iptcl] + a0*curtime;
	}
	else if(dadx > 0)
	{
		sinodt = sin(Omega*curtime);
		cosodt = cos(Omega*curtime);
		xout = A*cosodt + B*sinodt + x00;
		vxout = -A*Omega*sinodt + B*Omega*cosodt;
	}
	else
	{
		sinodt = sinh(Omega*curtime);
		cosodt = cosh(Omega*curtime);
		xout = A*cosodt + B*sinodt + x00;
		vxout = -A*Omega*sinodt - B*Omega*cosodt;
	}

	vyout = -omegac*(xout-x0[iptcl]) + vy0[iptcl];
	yout = -omegac*((A*sinodt - B*cosodt)/Omega + (x00 - x0[iptcl])*curtime)
			+ vy0[iptcl]*curtime + y0[iptcl];
	yout = yrt[nsteps*iptcl+iter-1]+vyout*pdata->dt;

}




void ExB_Initializer::finish(ParticleList* particles,
		HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo)
{


	printf("finishing");

	gnuplot_ctrl* plot = gnuplot_init();



	for(int i=0;i<(int)(fmin(2,pdata->my_nptcls));i++)
	{
		if(pdata->ndimensions == 1)
		{gnuplot_setstyle(plot,"points");
			gnuplot_plot_xy(plot,yt+nsteps*i,xt+nsteps*i,iter,"Tracks");
			gnuplot_setstyle(plot,"lines");
			gnuplot_plot_xy(plot,yrt+nsteps*i,xrt+nsteps*i,iter,"Tracks2");
		}
		else
		{gnuplot_setstyle(plot,"points");
			gnuplot_plot_xy(plot,xt+nsteps*i,yt+nsteps*i,iter,"Tracks");
			gnuplot_setstyle(plot,"lines");
			gnuplot_plot_xy(plot,xrt+nsteps*i,yrt+nsteps*i,iter,"Tracks2");
		}
	}

	getchar();



}











