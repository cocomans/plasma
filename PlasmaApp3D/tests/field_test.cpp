#include "../FieldDataCPU.h"
#include "../ParticleObjN.h"
#include "../PlasmaData.h"
#include "../ShapeFunctions.h"
#include <gnuplot_i.h>
#include <time.h>



float Ex_function(float x, float y, float z)
{
	return sin(2.0*pi_const*x);//+cos(2.0*pi_const*y)+cos(2.0*pi_const*z);
}

float Ey_function(float x, float y, float z)
{
	return cos(2.0*pi_const*x);//+sin(2.0*pi_const*y)+cos(2.0*pi_const*z);
}

float Ez_function(float x, float y, float z)
{
	return cos(2.0*pi_const*x);//+cos(2.0*pi_const*y)+sin(2.0*pi_const*z);
}

int main(int argc,char* argv[])
{
	PlasmaData pdata(argc,argv);


	FieldDataCPU fields;

	pdata.Lx = 1.0;
	pdata.Ly = 1.0;
	pdata.Lz = 1.0;

	pdata.setup();
	fields.allocate(&pdata);

	int nx = pdata.nx;
	int ny = pdata.ny;
	int nz = pdata.nz;
	int nshape = 100;

	float shape1[nshape];
	float shape2[nshape];
	float shapex[nshape];

	for(int i=0;i<nshape;i++)
	{
		float x = i*4.0/nshape - 2.0;

		shape1[i] = S1_shape(x);
		shape2[i] = S2_shape(x);
		shapex[i] = x;
	}


	gnuplot_ctrl* shape_plot = gnuplot_init();

	gnuplot_plot_xy(shape_plot,shapex,shape1,nshape,"Shape 1");
	gnuplot_plot_xy(shape_plot,shapex,shape2,nshape,"Shape 2");



	// Setup E-field
	for(int i=0;i<pdata.nx;i++)
	{
		for(int j=0;j<pdata.ny;j++)
		{
			for(int k=0;k<pdata.nz;k++)
			{
				float x = i*pdata.dxdi+pdata.xmin;
				float y = j*pdata.dydi+pdata.ymin;
				float z = k*pdata.dzdi+pdata.zmin;


				fields.getE(i,j,k,0) = Ex_function(x,y+0.5*pdata.dydi,z+0.5*pdata.dzdi);
				fields.getE(i,j,k,1) = Ey_function(x+0.5*pdata.dxdi,y,z+0.5*pdata.dzdi);
				fields.getE(i,j,k,2) = Ez_function(x+0.5*pdata.dxdi,y+0.5*pdata.dydi,z);

				fields.getB(i,j,k,0) = Ex_function(x,y,z);
				fields.getB(i,j,k,1) = Ey_function(x,y,z);
				fields.getB(i,j,k,2) = Ez_function(x,y,z);

				//fields.getB(i,j,k,0) = 0.0;
				//fields.getB(i,j,k,1) = 0.0f;
				//fields.getB(i,j,k,2) = 0.0f;

			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

	fields.init_plot();

	fields.plot(&pdata,pdata.nz/2,0,0,0);

	getchar();


	double Ex_total = 0;
	double Ey_total = 0;
	double Ez_total = 0;

	double Bx_total = 0;
	double By_total = 0;
	double Bz_total = 0;


	// Setup E-field
	for(int i=0;i<pdata.nx;i++)
	{
		for(int j=0;j<pdata.ny;j++)
		{
			for(int k=0;k<pdata.nz;k++)
			{
				float xfrac = 0.5;
				float yfrac = 0.5;
				float zfrac = 0.5;

				float x = (xfrac + i)*pdata.dxdi+pdata.xmin;
				float y = (yfrac + j)*pdata.dydi+pdata.ymin;
				float z = (zfrac + k)*pdata.dzdi+pdata.zmin;


				float Ex_intrp = fields.intrpE(xfrac,yfrac,zfrac,i,j,k,0,FieldData_deriv_f);
				float Ey_intrp = fields.intrpE(xfrac,yfrac,zfrac,i,j,k,1,FieldData_deriv_f);
				float Ez_intrp = fields.intrpE(xfrac,yfrac,zfrac,i,j,k,2,FieldData_deriv_f);

				float Ex_real = Ex_function(x,y,z);
				float Ey_real = Ey_function(x,y,z);
				float Ez_real = Ez_function(x,y,z);


				float Ex_err = fabs(Ex_real - Ex_intrp)/fabs(Ex_real);
				float Ey_err = fabs(Ey_real - Ey_intrp)/fabs(Ey_real);
				float Ez_err = fabs(Ez_real - Ez_intrp)/fabs(Ez_real);

				float Bx_intrp = fields.intrpB(xfrac,yfrac,zfrac,i,j,k,0,FieldData_deriv_f);
				float By_intrp = fields.intrpB(xfrac,yfrac,zfrac,i,j,k,1,FieldData_deriv_f);
				float Bz_intrp = fields.intrpB(xfrac,yfrac,zfrac,i,j,k,2,FieldData_deriv_f);

				float Bx_real = Ex_function(x,y,z);
				float By_real = Ey_function(x,y,z);
				float Bz_real = Ez_function(x,y,z);


				float Bx_err = fabs(Bx_real - Bx_intrp)/fabs(Bx_real);
				float By_err = fabs(By_real - By_intrp)/fabs(By_real);
				float Bz_err = fabs(Bz_real - Bz_intrp)/fabs(Bz_real);

				Ex_total += Ex_err;
				Ey_total += Ey_err;
				Ez_total += Ez_err;

				Bx_total += Bx_err;
				By_total += By_err;
				Bz_total += Bz_err;
				//printf("values: x: %f / %f, y: %f / %f, z: %f / %f \n",Ex_real,Ex_intrp,Ey_intrp,Ey_real,Ez_intrp,Ez_real);
				//printf("Errors: x = %e, y = %e, z = %e\n",Ex_err,Ey_err,Ez_err);

			}
		}
	}



	printf("Errors_avg: x = %e, y = %e, z = %e\n",Ex_total/(nx*ny*nz*1.0),Ey_total/(nx*ny*nz*1.0),Ez_total/(nx*ny*nz*1.0));
	printf("Errors_avg: x = %e, y = %e, z = %e\n",Bx_total/(nx*ny*nz*1.0),By_total/(nx*ny*nz*1.0),Bz_total/(nx*ny*nz*1.0));

	printf("Efield setup complete\n");


}
