#include "../FieldDataCPU.h"
#include "../ParticleObjN.h"
#include "../PlasmaData.h"
#include "../CPUTimer.h"
#include <gnuplot_i.h>
#include <time.h>

int main(int argc,char* argv[])
{
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

	float xmin = 0.0;
	float ymin = 0.0;
	float zmin = 0.0;

	float Lx = 1.0;
	float Ly = 1.0;
	float Lz = 1.0;

	int nx = 32;
	int ny = 32;
	int nz = 32;

	int nspecies = 1;

	const float time_total = 10.0;

	const float dtau0 = 0.01;

	const int nptcls = 32;
	const int steps = 1000;

	int iptcl[nptcls];

	float Ey = 0.5;
	float Bz = 100.0;


	PlasmaData pdata(argc,argv);

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
	pdata.epsilon_r = 1.0e-15;

	pdata.niter_max = 200;

	pdata.Bmag_avg = 1.0;
	pdata.ndimensions = 3;

	FieldDataCPU fields;
	ParticleObjN<nptcls> particles(iptcl);
	ParticleObjN<nptcls> particles0(iptcl);
	typevecN<float,nptcls> dtau;

	float x_plot[nptcls][steps];
	float y_plot[nptcls][steps];
	float gx_plot[nptcls][steps];
	float gy_plot[nptcls][steps];

	float x_plot_a[nptcls];
	float y_plot_a[nptcls];

	pdata.setup();
	fields.allocate(&pdata);

	fields.dx = pdata.dxdi;
	fields.dy = pdata.dydi;
	fields.dz = pdata.dzdi;

	particles.species = 0;



	for(int i=0;i<nptcls;i++)
	{
		iptcl[i] = i;

		particles.px(i) = rand()%99/100.0;
		particles.py(i) = rand()%99/100.0;
		particles.pz(i) = 0.5;

		particles.ix(i) = nx/2 + 2*rand()%2 -1;
		particles.iy(i) = ny/2 + 2*rand()%2 -1;
		particles.iz(i) = nz/2;

		particles.vx(i) = 0.1*(2*(rand()%10000))/10000.0+0.1;
		particles.vy(i) = 0.1*(2*(rand()%10000))/10000.0+0.1;
		particles.vz(i) = 0.1* (rand()%50 / 50.0f - 0.5);

		dtau(i) = 1.0/sqrt(pow(particles.vx(i)/pdata.dxdi,2.0f)
								+pow(particles.vy(i)/pdata.dydi,2.0f)
								+pow(particles.vz(i)/pdata.dzdi,2.0f));

		particles.dt_finished(i) = 0;

		dtau(i) = dtau0;

		printf("dtau = %f\n",dtau(i));
	}

	particles0 = particles;

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
				//float Ey = -1.1;


				fields.getE(i,j,k,0) = 0;
				fields.getE(i,j,k,1) = Ey;
				fields.getE(i,j,k,2) = 0.0;

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
	CPUTimer timer;

	float avg_error = 0.0;
	int n_error = 0;

	timer.start();
	for(int i=0;i<steps;i++)
	{
		//time = dtau0*(i);



		particles.PicardIterate(&pdata,&fields,dtau);



		//if((i+1)%64 == 0)
		//gnuplot_resetplot(plot_anim);

		float diff_avg = 0.0;

		for(int j=0;j<nptcls;j++)
		{

			float px,py,gx,gy;
			float rl;
			float vx,vy,vxy,vz,vxyz;

			float vgx,vgy;
			float verror;

			px = (particles.px(j) + particles.ix(j))*pdata.dxdi + pdata.xmin;
			py = (particles.py(j) + particles.iy(j))*pdata.dydi + pdata.ymin;

			vx = particles.vx(j);
			vy = particles.vy(j);
			vz = particles.vz(j);
			vxy = sqrt(vx*vx+vy*vy );

			vxyz = sqrt(vxy*vxy + vz*vz);

			rl = vxy/Bz;

			gx = vy*Bz/sqrt(vx*Bz*vx*Bz + vy*Bz*vy*Bz)*rl + px;
			gy = -vx*Bz/sqrt(vx*Bz*vx*Bz + vy*Bz*vy*Bz)*rl + py;

			x_plot[j][i] = px;
			y_plot[j][i] = py;

			gx_plot[j][i] = gx;
			gy_plot[j][i] = gy;

			if(i > 0)
			{
				vgx = (gx_plot[j][i] - gx_plot[j][0])/particles.dt_finished(j);
				vgy = (gy_plot[j][i] - gy_plot[j][0])/particles.dt_finished(j);

				verror = fabs(Ey/Bz - vgx)/(Ey/Bz);

				avg_error += verror;
				n_error ++;

				//printf("true v = %e, %e actual v = %e, %e, error = %e\n",Ey/Bz,0.0f,vgx,vgy,verror);
			}


			x_plot_a[j] = x_plot[j][i];
			y_plot_a[j] = y_plot[j][i];
		}

		particles.dt_finished += dtau;




		//gnuplot_plot_xy(plot_anim,x_plot_a,y_plot_a,nptcls,NULL);


	}
	timer.stop();

	printf("average error = %e \n",avg_error/((float)n_error));
	printf("Run did %f particles per second\n",nptcls*steps/(timer.diff()*1.0e-3));

	for(int j=0;j<nptcls;j++)
	{
		gnuplot_plot_xy(plot,x_plot[j],y_plot[j],steps,NULL);


	}
	printf("Press 'Enter' to continue\n");
		getchar();



	gnuplot_close(plot);

	gnuplot_close(plot_anim);

}
