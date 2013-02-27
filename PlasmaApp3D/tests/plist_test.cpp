#include "../FieldDataCPU.h"
#include "../ParticleListCPU.h"
#include "../PlasmaData.h"
#include "../CPUTimer.h"
#include "../HOMoments.h"
#include <omp.h>
#include <gnuplot_i.h>
#include <time.h>



int main(int argc,char* argv[])
{

	PlasmaData pdata(argc,argv);
	gnuplot_ctrl* plot;
	gnuplot_ctrl* plot_anim;


	plot = gnuplot_init();
	plot_anim = gnuplot_init();
	gnuplot_setstyle(plot,"lines");
	gnuplot_setstyle(plot_anim,"points");

	gnuplot_cmd(plot_anim,"set term gif animate nooptimize size 1280,1280 xffffffff");
	gnuplot_cmd(plot_anim,"set output \"particles.gif\"");

	gnuplot_cmd(plot_anim,"set xrange [-1:1]");
	gnuplot_cmd(plot_anim,"set yrange [-1:1]");

	float xmin = 0;
	float ymin = 0;
	float zmin = 0;

	float Lx = 5.0;
	float Ly = 5.0;
	float Lz = 5.0;

	int nx = 64;
	int ny = 64;
	int nz = 64;

	int nspecies = 1;

	const float dt = 0.01;

	const float dtau0 = 0.1;

	const int nptcls = 500;
	const int steps = 200;

	int iptcl[nptcls];

	float Ey = 5.0;
	float Bz = 100.0;




	pdata.nx = nx;
	pdata.ny = ny;
	pdata.nz = nz;

	pdata.Lx = Lx;
	pdata.Ly = Ly;
	pdata.Lz = Lz;

	pdata.xmin = xmin;
	pdata.ymin = ymin;
	pdata.zmin = zmin;
	pdata.epsilon_a = 1.0e-4;
	pdata.epsilon_r = 1.0e-10;

	pdata.dt = dt;

	pdata.niter_max = 20;

	pdata.nSubcycle_max = 1000;

	pdata.Bmag_avg = 1.0;
	pdata.ndimensions = 3;

	pdata.setup();

	FieldDataCPU fields;
	ParticleListCPU particles;
	HOMoments* moments;

	int numprocs = omp_get_num_procs();

	moments = (HOMoments*)malloc(numprocs*sizeof(HOMoments));

	for(int i=0;i<numprocs;i++)
	{
		moments[i] = *new HOMoments(&pdata);
	}
	float x_plot[nptcls][steps];
	float y_plot[nptcls][steps];
	float gx_plot[nptcls][steps];
	float gy_plot[nptcls][steps];

	float error_array[nptcls];


	//float x_plot_a[nptcls];
	//float y_plot_a[nptcls];


	fields.allocate(&pdata);
	particles.allocate(nptcls);

	fields.dx = pdata.dxdi;
	fields.dy = pdata.dydi;
	fields.dz = pdata.dzdi;

	particles.ispecies = 0;






	for(int i=0;i<nptcls;i++)
	{
		iptcl[i] = i;

		particles.px[i] = rand()%10000/10000.0;
		particles.py[i] = rand()%10000/10000.0;
		particles.pz[i] = 0.5;

		particles.ix[i] = nx/2;
		particles.iy[i] = ny/2;
		particles.iz[i] = nz/2;

		particles.vx[i] = 0.5*(2*(rand()%10000))/10000.0 + 0.5;
		particles.vy[i] = 0.5*(2*(rand()%10000))/10000.0 + 0.5;
		particles.vz[i] = 0.0* (rand()%50000 / 50000.0f - 0.5);

		error_array[i] = 0;


	}



	// Setup E-field
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int k=0;k<nz;k++)
			{
				float x = i*pdata.dxdi+xmin;
				float y = j*pdata.dydi+ymin;
				float z = k*pdata.dzdi+zmin;

				float Ex = -1.0*x;


				fields.getE(i,j,k,0) = 0;
				fields.getE(i,j,k,1) = Ey;
				fields.getE(i,j,k,2) = 0;

				fields.getB(i,j,k,0) = 0;
				fields.getB(i,j,k,1) = 0;
				fields.getB(i,j,k,2) = Bz;


			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

	fields.q2m[0] = 1.0;

	printf("Efield setup complete\n");

	float time;
	double avg_error = 0.0;
	int n_error = 0;

	CPUTimer timer;


	moments->init_plot();

	timer.start();
	for(int i=0;i<steps;i++)
	{
		//time = dtau0*(i);


		//moments.set_vals(0);
		particles.push(&pdata,&fields,moments);
		printf("finished step %i\n",i);


		for(int j=0;j<nptcls;j++)
		{

			float px,py,gx,gy;
			float rl;
			float vx,vy,vxy,vz,vxyz;

			float vgx,vgy;
			float verror;

			px = (particles.px[j] + particles.ix[j])*pdata.dxdi + pdata.xmin;
			py = (particles.py[j] + particles.iy[j])*pdata.dydi + pdata.ymin;

			vx = particles.vx[j];
			vy = particles.vy[j];
			vz = particles.vz[j];
			vxy = sqrt(vx*vx+vy*vy);

			vxyz = sqrt(vxy*vxy + vz*vz);

			rl = vxy/Bz;

			gx = vy*Bz/sqrt(vx*Bz*vx*Bz + vy*Bz*vy*Bz)*rl + px;
			gy = -vx*Bz/sqrt(vx*Bz*vx*Bz + vy*Bz*vy*Bz)*rl + py;

			x_plot[j][i] = px;
			y_plot[j][i] = py;

			gx_plot[j][i] = gx;
			gy_plot[j][i] = gy;

			if(i >= 1)
			{
				vgx = (gx_plot[j][i] - gx_plot[j][0])/(dt*(i));
				vgy = (gy_plot[j][i] - gy_plot[j][0])/(dt*(i));

				verror = fabs(Ey/Bz - vgx)/(Ey/Bz);

				error_array[j] = fmax(error_array[j],verror);

				avg_error += verror;
				n_error ++;

			//	printf("true[%i] v = %e, %e actual v = %e, %e, error = %e\n",
			//			j,Ey/Bz,0.0f,vgx,vgy,verror);
			}

		}




		//if((i+1)%64 == 0)
		//gnuplot_resetplot(plot_anim);
/*
		float diff_avg = 0.0;
		for(int j=0;j<nptcls;j++)
		{

			x_plot[j][i] = (particles.px[j] + particles.ix[j])*pdata.dxdi + pdata.xmin;
			y_plot[j][i] = (particles.py[j] + particles.iy[j])*pdata.dydi + pdata.ymin;

			//printf("particle %i with position %f, %f\n",j,x_plot[j][i],y_plot[j][i]);

		//	x_plot_a[j] = x_plot[j][i];
		//	y_plot_a[j] = y_plot[j][i];

		}
*/

		//avg_error += diff_avg / steps;


		//gnuplot_plot_xy(plot_anim,x_plot_a,y_plot_a,nptcls,NULL);


	}
	timer.stop();
	printf("average error = %e \n",avg_error/((float)n_error));
	printf("Run did %f particles per second\n",nptcls*steps/(timer.diff()*1.0e-3));

	for(int j=0;j<nptcls;j++)
	{
		if(error_array[j] >= 1.0e-2)
			gnuplot_plot_xy(plot,x_plot[j],y_plot[j],steps,NULL);


	}


	//moments->plot(nz/2,0,HOMoments_currentx);


	printf("Press 'Enter' to continue\n");
		getchar();

	moments->close_plot();



	gnuplot_close(plot);

	gnuplot_close(plot_anim);

}
