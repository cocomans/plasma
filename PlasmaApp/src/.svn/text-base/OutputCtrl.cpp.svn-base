/*
 * OutputCtrl.cpp
 *
 *  Created on: May 1, 2013
 *      Author: payne
 */

#include "OutputCtrl.h"
#include "NodeFieldData.h"
#include "NodeHOMoments.h"
#include "NodeParticleList.h"
#include "HOMomentsCPU.h"
#include "FieldDataCPU.h"
#include "CPUTimer.h"
#include "RunData.h"
#include <mpi.h>
#include "Util/mkpath.h"
#include "ParallelInfo.h"
#include "Util/webui/sse-server.h"
#include "Util/webui/LogStream.h"
#include "Util/webui/sseDataStream.h"

#include "Util/IPICAdaptor/IPICAdaptor.h"
#include "Util/webui/LogStream.h"




const int& OutputCtrl::istep(void){return pdata->rstatus->istep;}

OutputCtrl::OutputCtrl(PlasmaData* _pdata,int argc, char* argv[])
{
	pdata = _pdata;
	status = _pdata->rstatus;

	iplot = pdata->rdata->plot_flag;
	iplot = true;
	iSave_plots = true;
	iWrite_data = true;
	iFinalOutput = true;
	iExitSave = true;

	if((pdata->node_info->rank_g == 0)&&(iplot))
		iplot = true;
	else
		iplot = false;

//	iplot = false;
//	iSave_plots = false;
//	iWrite_data = false;
//	iFinalOutput = false;
//	iExitSave = false;

	OutputInterval = 5;
	PlotInterval = 5;
	PlotSaveInterval = 10;

	CurPlot = 0;

	push_timer = new CPUTimer();
	HOSolve_timer = new CPUTimer();
	Comm_timer = new CPUTimer();
	total_timer = new CPUTimer();
	LOSolve_timer = new CPUTimer();
	step_timer = new CPUTimer();
	coprocess_timer = new CPUTimer();



	charge_cons = (float*)malloc(2*pdata->nsteps*sizeof(float));
	energy_cons = (float*)malloc(2*pdata->nsteps*sizeof(float));
	momentum_cons = (float*)malloc(2*pdata->nsteps*sizeof(float));
	field_energy = (float*)malloc(2*pdata->nsteps*sizeof(float));
	Efield_energy = (float*)malloc(2*pdata->nsteps*sizeof(float));
	Bfield_energy = (float*)malloc(2*pdata->nsteps*sizeof(float));

	species_energy = (float**)malloc(pdata->nspecies*sizeof(float*));

	charge_diffs = (float**)malloc(2*sizeof(float*));
	current_diffs = (float**)malloc(6*sizeof(float*));
	ncalls_lotoho = 0;
	lotoho_array = (float*)malloc((pdata->nsteps+1)*pdata->nOuterPicard_max*sizeof(float));


	for(int i=0;i<2;i++)
	{
		charge_diffs[i] = (float*)malloc((pdata->nsteps+1)*pdata->nOuterPicard_max*sizeof(float));

	}

	for(int j=0;j<6;j++)
	{
		current_diffs[j] = (float*)malloc((pdata->nsteps+1)*pdata->nOuterPicard_max*sizeof(float));
	}

	for(int i=0;i<pdata->nspecies;i++)
		species_energy[i] = (float*)malloc(2*pdata->nsteps*sizeof(float));

	particle_energy = (float*)malloc(2*pdata->nsteps*sizeof(float));
	time_array = (float*)malloc(2*pdata->nsteps*sizeof(float));

	EfieldFile = new std::string("Efield");
	BfieldFile = new std::string("Bfield");

	nStreams = 4;
	nCharts = 4;


	lStreams = (LogStream**)malloc(nStreams*sizeof(LogStream*));
	lStreams[0] = &LO_log;
	lStreams[1] = &HO_log;
	lStreams[2] = &Outer_log;
	lStreams[3] = &debug_log;

	charts = (googleChart**)malloc(nCharts*sizeof(googleChart*));
	charts[0] = new googleChart("Conservation Quantities","time");
	charts[1] = new googleChart("Field Energy","px");
	charts[2] = new googleChart("Particle Energy","px");
	charts[3] = new googleChart("Moment Diffs","outer picard iterations");




	exit_control = 0;

	stream = new sseDataStream(pdata,nStreams,nCharts);
	control = new WebControlSystem(pdata,pdata->parser);

	stream->charts = charts;


	char* ROOT;
	int PORTi;
	pdata->parser->GetParam("SERVER_PORT_NUMBER","-p",PORTi,32000);

	std::stringstream tmp;
	tmp << PORTi;
	char PORTc[6];

	sprintf(PORTc,"%s",tmp.str().c_str());

	server = new sseServer(ROOT,PORTc,stream,control,&exit_control);

	debug_log.print("Starting webserver\n");
	if(pdata->node_info->rank_g == 0)
		server->StartServer();



}

void OutputCtrl::InitialOutput(NodeFieldData* fields,FieldDataCPU* fields_old,
		FieldDataCPU* fields_next,
		NodeHOMoments* moments,NodeParticleList* particles)
{


	if(iplot)
	{
		printf("allocating plots\n");

		n_cons_plots = 7+pdata->nspecies;

		if((pdata->ndimensions == 2)&&(pdata->nVelocity == 3))
		{

			n_field_plots = 8;
			n_moment_plots = 3*pdata->nspecies;
			n_particle_plots = pdata->nspecies;
		}
		else
		{
			n_cons_plots = 7+pdata->nspecies;
			n_field_plots = 8;
			n_moment_plots = 5;
			n_particle_plots = pdata->nspecies;

			if(pdata->nVelocity == 3)
			{
				n_cons_plots = 7+pdata->nspecies;
				n_field_plots = 11;
				n_moment_plots = 5;
				n_particle_plots = pdata->nspecies;
			}
		}

		plot_titles = (std::string***)malloc(4*sizeof(std::string**));

		plot_titles[0] = (std::string**)malloc(n_cons_plots*sizeof(std::string*));
		plot_titles[1] = (std::string**)malloc(n_field_plots*sizeof(std::string*));
		plot_titles[2] = (std::string**)malloc(n_moment_plots*sizeof(std::string*));
		plot_titles[3] = (std::string**)malloc(n_particle_plots*sizeof(std::string*));




		printf("Allocating Cons Plots\n");

		cons_plots = (gnuplot_ctrl**)malloc(n_cons_plots*sizeof(gnuplot_ctrl*));

		for(int i=0;i<n_cons_plots;i++)
		{

			plot_titles[0][i] = new std::string();
			cons_plots[i] = gnuplot_init();

		}

		gnuplot_cmd(cons_plots[2],"set log y");
		printf("Allocating Field Plots\n");

		field_plots = (gnuplot_ctrl**)malloc(n_field_plots*sizeof(gnuplot_ctrl*));

		plot_titles[1][0] = new std::string();
		field_plots[0] = gnuplot_init();
		gnuplot_setstyle(field_plots[0],"lines");

		for(int i=1;i<n_field_plots;i++)
		{
			plot_titles[1][i] = plot_titles[1][0];
			field_plots[i] = field_plots[0];
		}

		printf("Allocating Moment Plots\n");

		moment_plots = (gnuplot_ctrl**)malloc(n_moment_plots*sizeof(gnuplot_ctrl*));


		for(int i=0;i<n_moment_plots;i++)
		{
			plot_titles[2][i] = new std::string();
			moment_plots[i] = gnuplot_init();
			gnuplot_setstyle(moment_plots[i],"lines");
		}
		printf("Allocating Particle Plots\n");
		particle_plot = (gnuplot_ctrl**)malloc(n_particle_plots*sizeof(gnuplot_ctrl*));


		for(int i=0;i<n_particle_plots;i++)
		{
			plot_titles[3][i] = new std::string();
			particle_plot[i] = gnuplot_init();

			if(i==0)
			plot_titles[3][i] = new std::string("Electron Dist");
			else if(i==1)
			plot_titles[3][i] = new std::string("Ion Dist");

		}


		plot_titles[0][0][0] = std::string("Charge Conservation");
		plot_titles[0][1][0] = std::string("Energy Error");
		plot_titles[0][2][0] = std::string("Field Energy");
		plot_titles[0][3][0] = std::string("Particle Energy");

		plot_titles[0][4][0] = std::string("E-Field Energy");
		plot_titles[0][5][0] = std::string("B-Field Energy");

		plot_titles[0][6][0] = std::string("Momentum Conservation");


		plot_titles[0][7][0] = std::string("Electron Energy");

		if(pdata->nspecies > 1)
			plot_titles[0][8][0] = std::string("Ion Energy");



		if((pdata->ndimensions == 2)&&(pdata->nVelocity == 3))
		{

			plot_titles[1][3] = new std::string("Electric Field");
			plot_titles[1][7] = new std::string("Magnetic Field");

			plot_titles[2][0] = new std::string("Charge");
			plot_titles[2][1] = new std::string("Momentum");
			plot_titles[2][2] = new std::string("Kinetic Energy Density");
		}
		else
		{



			plot_titles[1][0] = new std::string("Electric Field");

			plot_titles[2][0] = new std::string("Charge");
			plot_titles[2][1] = new std::string("Momentum");
			plot_titles[2][2] = new std::string("Stress xx");

		}


	}

	printf("Checking Output Path\n");

	pdata->rdata->outpath = (char*)malloc(128*sizeof(char));

	sprintf(pdata->rdata->outpath,"./output/%s/%s",pdata->rdata->SimName,pdata->rdata->output_id);

	// Make the output path if it doesn't exist

	mkpath(pdata->rdata->outpath,0777);

	std::string outpath2(pdata->rdata->outpath);
	outpath2 += "/data";
	mkpath(outpath2.c_str(),0777);

	printf("Initial kinetic energy eval\n");
	kE0 = moments->evaluate_energy();
	printf("Initial potential energy eval\n");
	pE0 = fields_old->evaluate_energy();
	tE0 = kE0+pE0;
	printf("saving potential energy eval\n");
	charge_cons[0] = 0;
	energy_cons[0] = 0;
	field_energy[0] = pE0;
	particle_energy[0] = kE0;
	time_array[0] = 0;
	momentum_cons[0] = 0;

	Efield_energy[0] = fields_old->EFieldEnergy();
	Bfield_energy[0] = fields_old->BFieldEnergy();

	for(int i=0;i<pdata->nspecies;i++)
		species_energy[i][0] = moments->evaluate_energy(i);

//	SaveMesh();

//	for(int iField=0;iField<2;iField++)
//	{
//		std::string fieldname;
//		if(iField == 0)
//			fieldname = *EfieldFile;
//		else if(iField == 1)
//			fieldname = *BfieldFile;
//		else if(iField == 2)
//			fieldname = *AfieldFile;
//
//		std::string fName;
//
//		fName = std::string(pdata->rdata->outpath) + "/data/" + (fieldname) + ".dat";
//
//		std::ofstream oStr(fName.c_str(), std::ofstream::out);
//
//		if (!oStr) {
//			std::cout << "Could not open input file " << fName << std::endl;
//			exit(-1);
//		}
//
//		oStr.close();
//	}

	if(iplot)
	{
		if((pdata->ndimensions == 2)&&(pdata->nVelocity == 3))
		{
			PlotFields2D(fields->cpu_fields,3);
			PlotMoments2D(moments->moments_next,0);

		}
		else
		{
			PlotFields1D(fields->cpu_fields,0);
			PlotMoments1D(moments->moments_next);
			PlotParticles1D(particles);

			if(pdata->nVelocity > 1)
			{
				PlotFields1D(fields->cpu_fields,1);
				PlotFields1D(fields->cpu_fields,2);
				PlotFields1D(fields->cpu_fields,4);
				PlotFields1D(fields->cpu_fields,5);
				PlotFields1D(fields->cpu_fields,6);
				PlotFields1D(fields->cpu_fields,8);
				PlotFields1D(fields->cpu_fields,9);
				PlotFields1D(fields->cpu_fields,10);
			}
		}
	}
//	getchar();

	CurPlot++;

//	InitialCoprocess(fields, moments, particles);

}


void OutputCtrl::StepOutput(NodeHOMoments* moments,
		NodeParticleList* particles,FieldDataCPU* fields_old,FieldDataCPU* fields_next)
{


	double pEt,kEt;

	kEt = moments->evaluate_energy();
	particle_energy[istep()] = kEt;
	pEt = fields_next->evaluate_energy();
	field_energy[istep()] = pEt;
	energy_cons[istep()] = fabs((kEt+pEt - tE0));
	charge_cons[istep()] = moments->check_charge();
	momentum_cons[istep()] = moments->check_momentum(fields_next,fields_old);

	Efield_energy[istep()] = fields_next->EFieldEnergy();
	Bfield_energy[istep()] = fields_next->BFieldEnergy();

	for(int i=0;i<pdata->nspecies;i++)
		species_energy[i][istep()] = moments->evaluate_energy(i);

	printf("Field Energy = %e\n",Bfield_energy[istep()]);

	time_array[istep()] = istep()*pdata->dt;
	if(iplot)
	{
		if((pdata->ndimensions == 2)&&(pdata->nVelocity == 3))
		{
			PlotFields2D(fields_next,3);
			PlotMoments2D(moments->moments_old,0);

		}
		else
		{
			PlotFields1D(fields_next,0);
			PlotMoments1D(moments->moments_next);
			PlotParticles1D(particles);

			if(pdata->nVelocity > 1)
			{
				PlotFields1D(fields_next,1);
				PlotFields1D(fields_next,2);
				PlotFields1D(fields_next,4);
				PlotFields1D(fields_next,5);
				PlotFields1D(fields_next,6);
				PlotFields1D(fields_next,8);
				PlotFields1D(fields_next,9);
				PlotFields1D(fields_next,10);


			}
		}

		int offset = std::max(0,ncalls_lotoho-50);

		charts[3]->plotLines(8,std::min(ncalls_lotoho,50),lotoho_array+offset,
				"Density_e",charge_diffs[0]+offset,
				"Density_i",charge_diffs[1]+offset,
				"Momentum-x_e",current_diffs[0]+offset,
				"Momentum-y_e",current_diffs[1]+offset,
				"Momentum-z_e",current_diffs[2]+offset,
				"Momentum-x_i",current_diffs[3]+offset,
				"Momentum-y_i",current_diffs[4]+offset,
				"Momentum-z_i",current_diffs[5]+offset);

		charts[0]->plotLines(2,istep(),time_array,
				"Energy Conservation",energy_cons,
				"Charge Conservation",charge_cons);
		charts[2]->plotLines(1,istep()-1,time_array+1,
				"Electron Energy",species_energy[0]+1);
		charts[1]->plotLines(3,istep()-1,time_array+1,
				"Total Energy",field_energy+1,
				"B-Field Energy",Bfield_energy+1,
				"E-Field Energy",Efield_energy+1);
	}
	if((CurPlot > 0)&&(CurPlot%PlotInterval == 0))
	{
		if(iplot){
			for(int i=0;i<n_cons_plots;i++)
			{
				gnuplot_resetplot(cons_plots[i]);
			}
		}
	}

	CurPlot++;


	if((iplot)&&(CurPlot%PlotInterval == 0))
	{
		ChartInfo energy_cons_c = {"Energy Conservation",energy_cons};
		ChartInfo charge_cons_c = {"Charge Conservation",charge_cons};





		gnuplot_plot_xy(cons_plots[0],time_array,charge_cons,istep(),plot_titles[0][0]->c_str());
		gnuplot_plot_xy(cons_plots[1],time_array,energy_cons,istep(),plot_titles[0][1]->c_str());
		gnuplot_plot_xy(cons_plots[2],time_array+1,field_energy+1,istep()-1,plot_titles[0][2]->c_str());
		gnuplot_plot_xy(cons_plots[2],time_array+1,Efield_energy+1,istep()-1,plot_titles[0][4]->c_str());
		gnuplot_plot_xy(cons_plots[2],time_array+1,Bfield_energy+1,istep()-1,plot_titles[0][5]->c_str());


		gnuplot_plot_xy(cons_plots[3],time_array,particle_energy,istep(),plot_titles[0][3]->c_str());
		gnuplot_plot_xy(cons_plots[4],time_array,momentum_cons,istep(),plot_titles[0][6]->c_str());

		for(int i=0;i<pdata->nspecies;i++)
			gnuplot_plot_xy(cons_plots[5+i],time_array,species_energy[i],istep(),plot_titles[0][7+i]->c_str());


		if((std::string(pdata->rdata->SimName) == "TwoStream")||(std::string(pdata->rdata->SimName) == "TwoStream2D"))
			gnuplot_cmd(cons_plots[2],"replot \"E_vec.csv\" title \"Matlab Code\" with points");

		if(std::string(pdata->rdata->SimName) == "Weibel")
		{
			char temp[256];

			char temp2[256];


			double coe = 0.0;
			int istep_max = std::min(istep(),(int)ceil(15/pdata->dt));
			int istep_min = (int)floor(10/pdata->dt);

			for(int l=istep_min;l<istep_max;l++)
				coe += Bfield_energy[l]/exp(time_array[l]*0.2196*2);

			coe /= istep_max-istep_min;

//			sprintf(temp,"replot %e*exp(x*0.2196150*2) title \"Linear Fit\" with lines",coe);

//			gnuplot_cmd(cons_plots[2],temp);

			gnuplot_cmd(cons_plots[2],"set yrange [1e-7:]");
			gnuplot_cmd(cons_plots[2],"replot \"vxp/tplots.txt\" u 1:11 w l lw 4 t \"chen\" ");
			gnuplot_cmd(cons_plots[5],"replot \"vxp/tplots.txt\" u 1:5 w l lw 4 t \"chen\" ");
			gnuplot_cmd(cons_plots[6],"replot \"vxp/tplots.txt\" u 1:6 w l lw 4 t \"chen\" ");

		}

		if(istep()%PlotSaveInterval == 1)
		for(int i=0;i<5;i++)
		{
			std::string title = *(plot_titles[0][i]);

			std::string output = plotpath(title,istep());
			gnuplot_save_pdf(cons_plots[i],title.c_str());
		}
	}

	if(iplot&&(istep()%PlotSaveInterval == 1))
	{
		save_plots();



		SaveField(fields_old,istep(),0);
		SaveField(fields_next,istep()+1,0);

		SaveField(fields_old,istep(),1);
		SaveField(fields_next,istep()+1,1);

	}



//	StepCoprocess(istep(), moments, particles, fields_old, fields_next);


}

void OutputCtrl::FinalOutput(NodeHOMoments* moments,
		NodeParticleList* particles,FieldDataCPU* fields_old,FieldDataCPU* fields_next)
{
	particles->cpu_particles->piccard_stats(pdata);
	particles->cpu_particles->subcycle_stats(pdata);

	save_timing(particles);

	if((iplot))
	{

		gnuplot_plot_xy(cons_plots[0],time_array,charge_cons,istep(),plot_titles[0][0]->c_str());
		gnuplot_save_pdf(cons_plots[0],plotpath(plot_titles[0][0][0],0).c_str());
		gnuplot_plot_xy(cons_plots[1],time_array,energy_cons,istep(),plot_titles[0][1]->c_str());
		gnuplot_save_pdf(cons_plots[1],plotpath(plot_titles[0][1][0],0).c_str());
		gnuplot_plot_xy(cons_plots[2],time_array+1,field_energy+1,istep()-1,plot_titles[0][2]->c_str());
		gnuplot_plot_xy(cons_plots[2],time_array+1,Efield_energy+1,istep()-1,plot_titles[0][4]->c_str());
		gnuplot_plot_xy(cons_plots[2],time_array+1,Bfield_energy+1,istep()-1,plot_titles[0][5]->c_str());

		for(int i=0;i<pdata->nspecies;i++)
		{
			gnuplot_plot_xy(cons_plots[5+i],time_array,species_energy[i],istep(),plot_titles[0][7+i]->c_str());

		}

		if((std::string(pdata->rdata->SimName) == "TwoStream")||(std::string(pdata->rdata->SimName) == "TwoStream2D"))
			gnuplot_cmd(cons_plots[2],"replot \"E_vec.csv\" using 1:2 title \"Matlab Code\" with points");

		if(std::string(pdata->rdata->SimName) == "Weibel")
		{
			char temp[256];

			char temp2[256];


			double coe = 0.0;
			int istep_max = std::min(istep(),(int)ceil(15/pdata->dt));
			int istep_min = (int)floor(10/pdata->dt);

			for(int l=istep_min;l<istep_max;l++)
				coe += Bfield_energy[l]/exp(time_array[l]*0.2196*2);

			coe /= istep_max-istep_min;

//			sprintf(temp,"replot %e*exp(x*0.2196150*2) title \"Linear Fit\" with lines",coe);

//			gnuplot_cmd(cons_plots[2],temp);

			gnuplot_cmd(cons_plots[2],"set yrange [1e-7:]");
			gnuplot_cmd(cons_plots[2],"replot \"vxp/tplots.txt\" u 1:11 w l lw 4 t \"chen\" ");
			gnuplot_cmd(cons_plots[5],"replot \"vxp/tplots.txt\" u 1:5 w l lw 4 t \"chen\" ");
			gnuplot_cmd(cons_plots[6],"replot \"vxp/tplots.txt\" u 1:6 w l lw 4 t \"chen\" ");

		}


		gnuplot_save_pdf(cons_plots[2],plotpath(plot_titles[0][2][0],0).c_str());
		for(int i=0;i<pdata->nspecies;i++)
		{
			gnuplot_save_pdf(cons_plots[4+i],plotpath(plot_titles[0][6+i][0],0).c_str());

		}

		gnuplot_plot_xy(cons_plots[3],time_array,particle_energy,istep(),plot_titles[0][3]->c_str());
		gnuplot_save_pdf(cons_plots[3],plotpath(plot_titles[0][3][0],0).c_str());

		for(int i=0;i<pdata->nspecies;i++)
		{
			gnuplot_plot_xy(cons_plots[4+i],time_array,species_energy[i],istep(),plot_titles[0][6+i]->c_str());
			gnuplot_save_pdf(cons_plots[4+i],plotpath(plot_titles[0][6+i][0],0).c_str());

		}


//		gnuplot_cmd(cons_plots[2],"replot \"E_vec.csv\" title \"Matlab Code\" with points");

//		FinalCoprocess(moments, particles, fields_old, fields_next);

	}
}

void OutputCtrl::InitialCoprocess(
		NodeFieldData* fields,
		NodeHOMoments* moments,
		NodeParticleList* particles)
{
//	coprocessorInitialize(*(pdata->pythonScript));
//
//	coprocessorCreateGrid(
//		pdata->ndimensions, pdata->nVelocity,
//		pdata->nsteps, pdata->dt,
//		pdata->nptcls, pdata->nspecies,
//		pdata->nx, pdata->ny, pdata->nz,
//		pdata->xmin, pdata->ymin, pdata->zmin,
//		pdata->dxdi, pdata->dydi, pdata->dzdi,
//		charge_cons,
//		energy_cons,
//		field_energy,
//		particle_energy);
}

void OutputCtrl::StepCoprocess(
		int istep,
		NodeHOMoments* moments,
		NodeParticleList* particles,
		FieldDataCPU* fields_old,
		FieldDataCPU* fields_next)
{
//	coprocess_timer->start();
//	coprocessorFillParticles(particles->cpu_particles);
//	coprocessorFillField(fields_next, moments->moments_next);
//	coprocessorProcess(istep);
//	coprocess_timer->stop();
}

void OutputCtrl::FinalCoprocess(NodeHOMoments* moments,
		NodeParticleList* particles,FieldDataCPU* fields_old,FieldDataCPU* fields_next)
{
//	coprocessorFinalize();
}

void OutputCtrl::RecordMomentDiffs(double* charge_dif, double* current_dif)
{
	lotoho_array[ncalls_lotoho] = ncalls_lotoho;
	for(int i=0;i<pdata->nspecies;i++)
	{
		charge_diffs[i][ncalls_lotoho] = fmax(charge_dif[i],1.0e-17);
		for(int j=0;j<3;j++)
		{
			current_diffs[3*i + j][ncalls_lotoho] = fmax(current_dif[3*i+j],1.0e-17);

		}
	}

	ncalls_lotoho++;

}


void OutputCtrl::PlotFields1D(FieldDataCPU* fields,int icomponent)
{
	/*
	 * plane = 0: xy plane
	 * plane = 1: xz plane
	 * plane = 2: yz plane
	 * Fplot = 0: Efield
	 * Fplot = 1: Bfield
	 *
	 */
	int nxp = fields->nx;
	int i,j,k;

	float dxp,dyp;
	float x0,y0;

	float* x_vals;
	float* y_vals;
	float* y_vals_f;


	float* vals_in;



	x_vals = (float*)malloc(nxp*sizeof(float));
	y_vals = (float*)malloc(nxp*sizeof(float));

	y_vals_f = (float*)malloc(nxp*sizeof(float));

#pragma omp parallel for
	for(i=0;i<pdata->nx+1;i++)
	{
		x_vals[i] = pdata->dxdi*i + pdata->xmin;

		if(icomponent < 3)
			y_vals[i] = 0.25*(fields->getE(i-1,0,0,icomponent)
					+2.0*fields->getE(i,0,0,icomponent)
					+fields->getE(i+1,0,0,icomponent));
		else if(icomponent == 3)
			y_vals[i] = sqrt(
					powf(fields->getE(i,0,0,0),2.0)
					+ powf(fields->getE(i,0,0,1),2.0)
					+ powf(fields->getE(i,0,0,2),2.0));
		else if(icomponent < 7)
			y_vals[i] = 0.25*(fields->getB(i-1,0,0,icomponent-4)
					+2.0*fields->getB(i,0,0,icomponent-4)
					+fields->getB(i+1,0,0,icomponent-4));
		else if(icomponent == 7)
			y_vals[i] = sqrt(
					powf(fields->getB(i,0,0,0),2.0)
					+ powf(fields->getB(i,0,0,1),2.0)
					+ powf(fields->getB(i,0,0,2),2.0));
		else if(icomponent < 11)
		{
			y_vals[i] = 0.25*(fields->getA(i-1,0,0,icomponent-8)
					+2.0*fields->getA(i,0,0,icomponent-8)
					+fields->getA(i+1,0,0,icomponent-8));

		}
		else if(icomponent == 11)
		{	y_vals[i] = sqrt(
					powf(fields->getA(i,j,0,0),2.0)
					+ powf(fields->getA(i,j,0,1),2.0)
					+ powf(fields->getA(i,j,0,2),2.0));
		}


	}


	std::string Title;
	switch(icomponent)
	{
	case 0: Title = "Ex"; break;
	case 1: Title = "Ey"; break;
	case 2: Title = "Ez"; break;
	case 3: Title = "Et"; break;
	case 4: Title = "Bx"; break;
	case 5: Title = "By"; break;
	case 6: Title = "Bz"; break;
	case 7: Title = "Bt"; break;
	case 8: Title = "Ax"; break;
	case 9: Title = "Ay"; break;
	case 10: Title = "Az"; break;
	case 11: Title = "At"; break;

	}


	if((CurPlot > 0)&&(icomponent==0))
		gnuplot_resetplot(field_plots[icomponent]);

	gnuplot_plot_xy(field_plots[icomponent],x_vals,y_vals,pdata->nx+1,Title.c_str());


	free(x_vals);
	free(y_vals);
}


void OutputCtrl::PlotMoments1D(HOMomentsCPU* moments)
{
	/*
	 * plane = 0: xy plane
	 * plane = 1: xz plane
	 * plane = 2: yz plane
	 * Fplot = 0: Efield
	 * Fplot = 1: Bfield
	 *
	 */
	int nxp = pdata->nx;
	int i,j,k;

	float dxp,dyp;
	float x0,y0;

	float* x_vals;

	float* vals_in;


	x_vals = (float*)malloc(nxp*sizeof(float));
	float* charge = (float*)malloc(nxp*sizeof(float));
	float* current = (float*)malloc(nxp*sizeof(float));
	float* S2xx = (float*)malloc(nxp*sizeof(float));

#pragma omp parallel for
	for(i=0;i<pdata->nx;i++)
	{
		x_vals[i] = pdata->dxdi*i + pdata->xmin;


		charge[i] 	= 	moments->get_val(i,0,0,0,HOMoments_charge);
		current[i]  = 	moments->get_val(i,0,0,0,HOMoments_currentx);
		S2xx[i] 	= 	moments->get_val(i,0,0,0,HOMoments_S2xx);



	}


	if(CurPlot > 0)
	{
		gnuplot_resetplot(moment_plots[0]);
		gnuplot_resetplot(moment_plots[1]);
		gnuplot_resetplot(moment_plots[2]);
	}

	gnuplot_plot_xy(moment_plots[0],x_vals,charge,pdata->nx,"Electron Density");
	gnuplot_plot_xy(moment_plots[1],x_vals,current,pdata->nx,"E current");
	gnuplot_plot_xy(moment_plots[2],x_vals,S2xx,pdata->nx,"S2xx");
	if(pdata->nspecies > 1){
	for(i=0;i<pdata->nx;i++)
	{
		x_vals[i] = pdata->dxdi*i + pdata->xmin;


		charge[i] 	= 	moments->get_val(i,0,0,1,HOMoments_charge);
		current[i]  = 	moments->get_val(i,0,0,1,HOMoments_currentx);
		S2xx[i] 	= 	moments->get_val(i,0,0,1,HOMoments_S2xx);



	}

	gnuplot_plot_xy(moment_plots[0],x_vals,charge,pdata->nx,"Ion Density");
	gnuplot_plot_xy(moment_plots[1],x_vals,current,pdata->nx,"Ion current");
	gnuplot_plot_xy(moment_plots[2],x_vals,S2xx,pdata->nx,"S2xx");
	}



	free(x_vals);
	free(charge);
	free(current);
	free(S2xx);
}

void OutputCtrl::PlotParticles1D(NodeParticleList* particles)
{

	for(int i=0;i<n_particle_plots;i++)
	{
		ParticleListCPU* particlelist = particles->cpu_particles + i;
		int nptcls = particlelist->nptcls;
		float* x_vals = (float*)malloc(nptcls*sizeof(float));
		float* y_vals = (float*)malloc(nptcls*sizeof(float));



#pragma omp parallel for
		for(int j=0;j<nptcls;j++)
		{


			x_vals[j] = pdata->dxdi*(particlelist->px[j]
					+ particlelist->ix[j]) + pdata->xmin;

			y_vals[j] = particlelist->vx[j];
		}

		if(CurPlot > 0)
		{
			gnuplot_resetplot(particle_plot[i]);
//			charts[1]->appendPoints(3,nptcls,x_vals,
//					"x velocity",y_vals,
//					"time",time_vals,
//					"species",ispecies);
		}




		gnuplot_plot_xy(particle_plot[i],x_vals,y_vals,nptcls,plot_titles[3][i]->c_str());

		free(x_vals);
		free(y_vals);
	}

}

void OutputCtrl::PlotFields2D(FieldDataCPU* fields,int icomponent)
{
	/*
	 * plane = 0: xy plane
	 * plane = 1: xz plane
	 * plane = 2: yz plane
	 * Fplot = 0: Efield
	 * Fplot = 1: Bfield
	 *
	 */
	int nxp = fields->nx+1;
	int nyp = fields->ny+1;
	int i,j,k;
	k = 0;


	float dxp,dyp;
	float x0,y0;

	float* x_vals;
	float* y_vals;
	float* z_vals;

	float* dx_vals;
	float* dy_vals;
	float* dz_vals;

	float* vals_in;



	x_vals = (float*)malloc(nxp*nyp*sizeof(float));
	y_vals = (float*)malloc(nxp*nyp*sizeof(float));
	z_vals = (float*)malloc(nxp*nyp*sizeof(float));

	dx_vals = (float*)malloc(nxp*nyp*sizeof(float));
	dy_vals = (float*)malloc(nxp*nyp*sizeof(float));
	dz_vals = (float*)malloc(nxp*nyp*sizeof(float));

	float scale = sqrt(pdata->Lx*pdata->Ly/(1.0*pdata->nx*pdata->ny))/2.0;

	for(j=0;j<=pdata->ny;j++)
	{

		for(i=0;i<=pdata->nx;i++)
		{
			y_vals[j] = pdata->dydi*j + pdata->ymin;
			x_vals[i] = pdata->dxdi*i + pdata->xmin;


			if(icomponent < 3)
				z_vals[i+nxp*j] = fields->getE(i,0,0,icomponent);
			else if(icomponent == 3)
			{	z_vals[i+nxp*j] = sqrt(
						powf(fields->getE(i,j,0,0),2.0)
						+ powf(fields->getE(i,j,0,1),2.0)
						+ powf(fields->getE(i,j,0,2),2.0));
				dx_vals[i+nxp*j] = fields->getE(i,j,k,0)*scale/z_vals[i+nxp*j];
				dy_vals[i+nxp*j] = fields->getE(i,j,k,1)*scale/z_vals[i+nxp*j];
				dz_vals[i+nxp*j] = fields->getE(i,j,k,2)*scale/z_vals[i+nxp*j];
			}
			else if(icomponent < 7)
			{
				z_vals[i+nxp*j] = fields->getB(i,0,0,icomponent-4);
			}
			else if(icomponent == 7)
			{	z_vals[i+nxp*j] = sqrt(
						powf(fields->getB(i,j,0,0),2.0)
						+ powf(fields->getB(i,j,0,1),2.0)
						+ powf(fields->getB(i,j,0,2),2.0));
				dx_vals[i+nxp*j] = fields->getB(i,j,k,0)*scale/z_vals[i+nxp*j];
				dy_vals[i+nxp*j] = fields->getB(i,j,k,1)*scale/z_vals[i+nxp*j];
				dz_vals[i+nxp*j] = fields->getB(i,j,k,2)*scale/z_vals[i+nxp*j];
			}
			else if(icomponent < 11)
			{
				z_vals[i+nxp*j] = fields->getA(i,0,0,icomponent-8);

			}
			else if(icomponent == 11)
			{	z_vals[i+nxp*j] = sqrt(
						powf(fields->getA(i,j,0,0),2.0)
						+ powf(fields->getA(i,j,0,1),2.0)
						+ powf(fields->getA(i,j,0,2),2.0));
				dx_vals[i+nxp*j] = fields->getA(i,j,k,0)*scale/z_vals[i+nxp*j];
				dy_vals[i+nxp*j] = fields->getA(i,j,k,1)*scale/z_vals[i+nxp*j];
				dz_vals[i+nxp*j] = fields->getA(i,j,k,2)*scale/z_vals[i+nxp*j];
			}


		}

	}
	if(CurPlot > 0)
		gnuplot_resetplot(field_plots[icomponent]);


	if((istep()%PlotSaveInterval == 1)||(istep()>=pdata->nsteps))
	{
		std::string title = *(plot_titles[1][icomponent]);
		std::string output = plotpath(title,istep());

		char line[128];

		gnuplot_cmd(field_plots[icomponent],"set term pdf");

		sprintf(line,"set output \"%s.pdf\"",output.c_str());

		gnuplot_cmd(field_plots[icomponent],line);

		gnuplot_plot_xyz(field_plots[icomponent],x_vals,y_vals,z_vals,nxp,nyp,"Fields");
		gnuplot_cmd(field_plots[icomponent],"unset hidden3d");
		gnuplot_cmd(field_plots[icomponent],"set view map");
		gnuplot_plot_vector(field_plots[icomponent],x_vals,y_vals,z_vals,dx_vals,dy_vals,nxp,nyp,"Fields");

		gnuplot_cmd(field_plots[icomponent],"set term pop");
		gnuplot_cmd(field_plots[icomponent],"set out");
	}


	gnuplot_plot_xyz(field_plots[icomponent],x_vals,y_vals,z_vals,nxp,nyp,"Fields");
	gnuplot_cmd(field_plots[icomponent],"unset hidden3d");
	gnuplot_cmd(field_plots[icomponent],"set view map");
	gnuplot_plot_vector(field_plots[icomponent],x_vals,y_vals,z_vals,dx_vals,dy_vals,nxp,nyp,"Fields");



	free(x_vals);
	free(y_vals);
	free(z_vals);

	free(dx_vals);
	free(dy_vals);
	free(dz_vals);
}


void OutputCtrl::PlotMoments2D(HOMomentsCPU* moments,int ispecies)
{
	int nxp = pdata->nx+1;
	int nyp = pdata->ny+1;
	int i,j,k;
	k = 0;


	float dxp,dyp;
	float x0,y0;

	float* x_vals;
	float* y_vals;
	float* z_vals[3];

	float* dx_vals[2];
	float* dy_vals[2];
	float* dz_vals[2];

	float* vals_in;

	enum HOMoments_moment imom[2][3] = {
			{HOMoments_currentx,HOMoments_currenty,HOMoments_currentz},
			{HOMoments_S2xx,HOMoments_S2yy,HOMoments_S2zz}};



	x_vals = (float*)malloc(nxp*nyp*sizeof(float));
	y_vals = (float*)malloc(nxp*nyp*sizeof(float));

	for(int l=0;l<3;l++)
	z_vals[l] = (float*)malloc(nxp*nyp*sizeof(float));


	for(int l=0;l<2;l++)
	{
		dx_vals[l] = (float*)malloc(nxp*nyp*sizeof(float));
		dy_vals[l] = (float*)malloc(nxp*nyp*sizeof(float));
		dz_vals[l] = (float*)malloc(nxp*nyp*sizeof(float));
	}
	double scale = 1.0*sqrt(pdata->Lx*pdata->Ly/(pdata->nx*pdata->ny));

	for(j=0;j<=pdata->ny;j++)
	{

		for(i=0;i<=pdata->nx;i++)
		{
			y_vals[j] = pdata->dydi*j + pdata->ymin;
			x_vals[i] = pdata->dxdi*i + pdata->xmin;


			z_vals[0][i+nxp*j] = moments->get_val(i,j,0,ispecies,HOMoments_charge);

			for(int l=0;l<2;l++)
			{
				double zt = sqrt(
						powf(moments->get_val(i,j,0,ispecies,imom[l][0]),2.0)
						+ powf(moments->get_val(i,j,0,ispecies,imom[l][1]),2.0)
						+ powf(moments->get_val(i,j,0,ispecies,imom[l][2]),2.0));

				z_vals[l+1][i+nxp*j] = zt;

				dx_vals[l][i+nxp*j] = moments->get_val(i,j,0,ispecies,imom[l][0])*scale/zt;
				dy_vals[l][i+nxp*j] = moments->get_val(i,j,0,ispecies,imom[l][1])*scale/zt;
				dz_vals[l][i+nxp*j] = moments->get_val(i,j,0,ispecies,imom[l][2])*scale/zt;

			}


		}

	}

	for(int l=0;l<2;l++)
	{
		int lplot = 3*ispecies+l;
		if(CurPlot > 0)
		{
			gnuplot_resetplot(moment_plots[lplot]);
		}

		if(l == 0)
		{
			// Density Plot
			if((istep()%PlotSaveInterval == 1)||(istep()>=pdata->nsteps))
			{
				std::string title = *(plot_titles[2][lplot]);
				std::string output = plotpath(title,istep());

				char line[128];

				gnuplot_cmd(moment_plots[lplot],"set term pdf");

				sprintf(line,"set output \"%s.pdf\"",output.c_str());

				gnuplot_cmd(moment_plots[lplot],line);
				gnuplot_cmd(moment_plots[lplot],"unset hidden3d");
				gnuplot_cmd(moment_plots[lplot],"set view map");
				gnuplot_plot_xyz(moment_plots[lplot],x_vals,y_vals,z_vals[l],nxp,nyp,"Density");


				gnuplot_cmd(moment_plots[lplot],"set term pop");
				gnuplot_cmd(moment_plots[lplot],"set out");
			}
			gnuplot_cmd(moment_plots[lplot],"unset hidden3d");
			gnuplot_cmd(moment_plots[lplot],"set view map");
			gnuplot_plot_xyz(moment_plots[lplot],x_vals,y_vals,z_vals[l],nxp,nyp,"Density");


		}
		else
		{
			if((istep()%PlotSaveInterval == 1)||(istep()>=pdata->nsteps))
			{
				std::string title = *(plot_titles[2][lplot]);
				std::string output = plotpath(title,istep());

				char line[128];

				gnuplot_cmd(moment_plots[lplot],"set term pdf");

				sprintf(line,"set output \"%s.pdf\"",output.c_str());

				gnuplot_cmd(moment_plots[lplot],line);

				gnuplot_plot_xyz(moment_plots[lplot],x_vals,y_vals,z_vals[l],nxp,nyp,"Current");
				gnuplot_cmd(moment_plots[lplot],"unset hidden3d");
				gnuplot_cmd(moment_plots[lplot],"set view map");
				gnuplot_plot_vector(moment_plots[lplot],x_vals,y_vals,
						z_vals[l],dx_vals[l-1],dy_vals[l-1],nxp,nyp,"Current");

				gnuplot_cmd(moment_plots[lplot],"set term pop");
				gnuplot_cmd(moment_plots[lplot],"set out");
			}


			gnuplot_plot_xyz(moment_plots[lplot],x_vals,y_vals,z_vals[l],nxp,nyp,"Currents");
			gnuplot_cmd(moment_plots[lplot],"unset hidden3d");
			gnuplot_cmd(moment_plots[lplot],"set view map");
			gnuplot_plot_vector(moment_plots[lplot],x_vals,y_vals,
					z_vals[l],dx_vals[l-1],dy_vals[l-1],nxp,nyp,"Currents");

		}

	}

	free(x_vals);
	free(y_vals);
	for(int l=0;l<3;l++)
		free(z_vals[l]);

	for(int l=0;l<2;l++)
	{
		free(dx_vals[l]);
		free(dy_vals[l]);
		free(dz_vals[l]);
	}

}

void OutputCtrl::PlotParticles2D(NodeParticleList* particles)
{

	for(int i=0;i<n_particle_plots;i++)
	{
		ParticleListCPU* particlelist = particles->cpu_particles + i;
		int nptcls = particlelist->nptcls;
		float* x_vals = (float*)malloc(nptcls*sizeof(float));
		float* y_vals = (float*)malloc(nptcls*sizeof(float));


		for(int j=0;j<nptcls;j++)
		{
			x_vals[j] = pdata->dxdi*(particlelist->px[j]
					+ particlelist->ix[j]) + pdata->xmin;

			y_vals[j] = particlelist->vx[j];
		}

		if(CurPlot > 0)
		{
			gnuplot_resetplot(particle_plot[i]);
		}

		gnuplot_plot_xy(particle_plot[i],x_vals,y_vals,nptcls,plot_titles[3][i]->c_str());

		free(x_vals);
		free(y_vals);
	}

}


void OutputCtrl::save_plots(void)
{
	printf("saving plots\n");

	gnuplot_ctrl** all_plots[4] = {cons_plots,field_plots,moment_plots,particle_plot};
	int n_plots[4] = {0,n_field_plots,n_moment_plots,n_particle_plots};

	for(int i=0;i<4;i++)
	{
		for(int j=0;j<n_plots[i];j++)
		{
			std::string title = *(plot_titles[i][j]);

			std::string output = plotpath(title,istep());

			gnuplot_save_pdf((all_plots[i])[j],output.c_str());

		}
	}

}

std::string OutputCtrl::plotpath(std::string plot_name,int plot_number)
{
	std::string result = pdata->rdata->outpath;

	result = result + "/" + plot_name  + std::to_string((long long int)plot_number);
	return result;
}





double2 calc_perf(CPUTimer* timer,double weight)
{
	double2 result; // result, std

	double temp = timer->get_cummulative();
	double temp2 = temp*temp;

	double t1,t2;

	MPI_Allreduce(&temp,&t1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	MPI_Allreduce(&temp2,&t2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	result.x = t1*weight;
	result.y = sqrt(fabs(t1*t1 - t2))*weight;

	return result;

}

double2 calc_particle_perf(PlasmaData* pdata,ParticleListCPU* plist,int itime,int idevice)
{
	double2 result; // result, std

	int num_devices;
	if(idevice == 1)
		num_devices = 1;
	else
		num_devices = pdata->node_info->nTasks_g;

	double temp;
	if(plist->device_type == idevice)
		temp = plist->get_cummulative_time(itime)/num_devices;
	else
		temp = 0;

	double temp2 = temp*temp;

	double t1,t2;

	MPI_Allreduce(&temp,&t1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	MPI_Allreduce(&temp2,&t2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	result.x = t1;
	result.y = sqrt(fabs(t1*t1 - t2));

	return result;

}

std::string OutputCtrl::calc_timing(double perf_weight,const char* legend)
{
	double2 HOSolve_perf,push_perf,communication_perf,LOSolve_perf,total_perf,step_perf;

	double temp;

	double weight = 1.0e6/(perf_weight);

	HOSolve_perf = calc_perf(HOSolve_timer,weight);
	push_perf = calc_perf(push_timer,weight);
	communication_perf = calc_perf(Comm_timer,weight);
	step_perf = calc_perf(step_timer,weight);

	//communication_perf = HOSolve_perf - push_perf;

	if(pdata->rdata->lo_all)
		LOSolve_perf = calc_perf(LOSolve_timer,weight);
	else{
		LOSolve_perf.x = weight*LOSolve_timer->get_cummulative();
		LOSolve_perf.y = 0;
	}
	total_perf.x = weight*total_timer->get_cummulative();
	total_perf.y = 0;

	char tempc[128];
	std::string result("");

	sprintf(tempc,"Performance in %s\n",legend);
	result += tempc;
	sprintf(tempc,"HOSolve step took: %e +- %e\n",HOSolve_perf.x,HOSolve_perf.y);
	result += tempc;
	sprintf(tempc,"push step took: %e +- %e\n",push_perf.x,push_perf.y);
	result += tempc;
	sprintf(tempc,"Communication took: %e +- %e\n",communication_perf,communication_perf.y);
	result += tempc;
	sprintf(tempc,"LOSolve step took: %e +- %e\n",LOSolve_perf.x,LOSolve_perf.y);
	result += tempc;
	sprintf(tempc,"Total step took: %e +- %e\n",total_perf.x,total_perf.y);
	result += tempc;


	return result;
}

double OutputCtrl::total_substeps()
{
	double result;

	MPI_Allreduce(&(status->nsteps_node),&result,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	return result;
}


void OutputCtrl::save_timing(NodeParticleList* particles)
{
	//double4 piccard_stats;
	//double4 subcycle_stats;
	//subcycle_stats = particles->subcycle_stats(pdata);
	//piccard_stats = particles_old[0] -> piccard_stats(pdata);

	double nsubsteps = total_substeps();

	double2 HOSolve_perf,push_perf,communication_perf,LOSolve_perf,total_perf,step_perf;

	double perf_weight = 1.0e6/(nsubsteps);

	double temp;


	HOSolve_perf = calc_perf(HOSolve_timer,perf_weight);
	push_perf = calc_perf(push_timer,perf_weight);
	communication_perf = calc_perf(Comm_timer,perf_weight);
	step_perf = calc_perf(step_timer,perf_weight);



	//communication_perf = HOSolve_perf - push_perf;

	if(pdata->rdata->lo_all)
			LOSolve_perf = calc_perf(LOSolve_timer,perf_weight);
		else{
			LOSolve_perf.x = perf_weight*LOSolve_timer->get_cummulative();
			LOSolve_perf.y = 0;
		}
		total_perf.x = perf_weight*total_timer->get_cummulative();
		total_perf.y = 0;


	if(pdata->node_info->rank_g == 0){
	printf("Sim Params:\n");
	printf("Num Cores: %i\n",pdata->node_info->nCPU);
	printf("Num Nodes: %i\n",pdata->node_info->nTasks_g);
	printf("CPU Vec Length: %i\n",pdata->node_info->cpu_info->vec_length);
	printf("Num Ptcls per Node: %i\n",pdata->nptcls);
	printf("Nx : %i\n",pdata->nx);
	printf("Total Subcycle Steps: %e\n", (double)nsubsteps);
	printf("Total number of outer piccard: %i\n",status->npiccard_outer);

	printf("Performance in ns/particle-substep\n");
	printf("HOSolve step took: %f +- %f\n",HOSolve_perf.x,HOSolve_perf.y);
	printf("push step took: %f +- %f\n",push_perf.x,push_perf.y);
	printf("Communication took: %f +- %f\n",communication_perf,communication_perf.y);
	printf("LOSolve step took: %f +- %f\n",LOSolve_perf.x,LOSolve_perf.y);
	printf("Total step took: %f +- %f\n",total_perf.x,total_perf.y);

	}

	std::string piccard_perf = calc_timing((status->npiccard_outer*pdata->node_info->nTasks_g),
			"ns/outer piccard");

	if(pdata->node_info->rank_g == 0)
		printf("%s\n",piccard_perf.c_str());

	MPI_Barrier(MPI_COMM_WORLD);
	//printf("Node %i took %e ms to push %e steps\n",pdata->mynode,push_timer.get_cummulative(),(double)nsteps_node);


	int num_gpus;
	int igpu = 0;

//	if(pdata->device_type == 1)
//		igpu = 1;

	MPI_Reduce(&igpu,&num_gpus,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

	double np = status->npiccard_outer/((double)pdata->nsteps);
	double np2 = status->npiccard_outer2/((double)pdata->nsteps);
	double npiccard_std = sqrt(fabs(np*np - np2));
//	MPI_Barrier(MPI_COMM_WORLD);
//	printf("Node %i took %e ms to HO\n",pdata->mynode,HOSolve_timer.get_cummulative());
//




//	char filename[60];
//	FILE* fp;
//
//	int nparams = 19;
//	int ntimes;
//	char* run_params[nparams];
//	char* time_names[ntimes];
//
//
//
//	run_params[0] = "dt";
//	run_params[1] = "nsteps";
//	run_params[2] = "nx";
//	run_params[3] = "ny";
//	run_params[4] = "nz";
//	run_params[5] = "n_electrons";
//	run_params[6] = "n_ions";
//	run_params[7] = "nptcls_gpu";
//	run_params[8] = "nptcls_cpu";
//	run_params[9] = "ndimensions";
//
//	run_params[10] = "Lx";
//	run_params[11] = "Ly";
//	run_params[12] = "Lz";
//	run_params[13] = "Problem_Type";
//	run_params[14] = "mass_ratio";
//	run_params[15] = "cpu_vec_length";
//	run_params[16] = "num_cores";
//	run_params[17] = "num_nodes";
//
//	run_params[18] = "Num_SubCycle_Steps";
//
//
//
//	time_names[0] = "HOSolve_time";
//	time_names[1] = "Push_time";
//	time_names[2] = "Communication_time";
//	time_names[3] = "LOSolve_time";
//	time_names[4] = "Step_time";
//	time_names[5] = "Total_time";

	int nparticle_times = 10;
	double2 particle_times[nparticle_times];
	for(int i=0;i<nparticle_times;i++)
	{
		particle_times[i] = calc_particle_perf(pdata,particles->cpu_particles,i,0);
	}

	if(pdata->node_info->rank_g == 0)
	{
		char filename[60];
		FILE* fp;

		int nparams = 14;
		int ntimes = 16;
		const char* run_params[nparams];
		const char* time_names[ntimes];

		double2 run_times[ntimes];
		double2 adjusted_times[ntimes];



		run_params[0] = "dt";
		run_params[1] = "nsteps";
		run_params[2] = "nsubsteps_total";
		run_params[3] = "npiccard_total";
		run_params[4] = "nptcls_total";
		run_params[5] = "ncells";
		run_params[6] = "npiccard_outer";
		run_params[7] = "npiccard_std";
		run_params[8] = "vector_length";
		run_params[9] = "num_cores";

		run_params[10] = "num_nodes";
		run_params[11] = "num_gpus";
		run_params[12] = "lo-all";
		run_params[13] = "OutputID";



		time_names[0] = "HOSolve_time";
		time_names[1] = "Push_time";
		time_names[2] = "Comm_time";
		time_names[3] = "LOSolve_time";
		time_names[4] = "Step_time";
		time_names[5] = "Total_time";
		time_names[6] = "Ppicard_time";
		time_names[7] = "Accel_time";
		time_names[8] = "Tally_time";
		time_names[9] = "Crossing_time";
		time_names[10] = "Dtau_est_time";
		time_names[11] = "ChargeS2Tally_time";
		time_names[12] = "PLoadStore_time";
		time_names[13] = "PPiccardOther_time";
		time_names[14] = "Push2_time";
		time_names[15] = "PushOther_time";



		adjusted_times[0] = HOSolve_perf;
		adjusted_times[1] = push_perf;
		adjusted_times[2] = communication_perf;
		adjusted_times[3] = LOSolve_perf;
		adjusted_times[4] = step_perf;
		adjusted_times[5] = total_perf;





		for(int i=0;i<ntimes;i++)
		{
			run_times[i].x = adjusted_times[i].x / (perf_weight) ;
			run_times[i].y = adjusted_times[i].y / (perf_weight) ;

		}

		for(int i=0;i<nparticle_times;i++)
		{
			run_times[i+6] = particle_times[i];
			adjusted_times[i+6].x = run_times[i+6].x*1.0e6/nsubsteps;
			adjusted_times[i+6].y = run_times[i+6].y*1.0e6/nsubsteps;
		}

		// Check the output directory
		mkpath("./benchmarks",0777);

		// Setup the filename
		sprintf(filename,"./benchmarks/benchmark%i.dat",pdata->rdata->runid);

		// Check to see if the file exists
		fp = fopen(filename,"r+");

		// If the file doesn't exist yet, create it and write the top line
		if(fp == NULL)
		{

			fp = fopen(filename,"w");
			char header[nparams*16+35*(ntimes)];
			for(int i=0;i<nparams;i++)
			{
				fprintf(fp,"%s,",run_params[i]);
			}

			for(int i=0;i<ntimes;i++)
			{
				fprintf(fp,"%s(ms), %s(ms), %s(ns), %s(ns),",time_names[i],"std","adjusted","std");
			}

			fprintf(fp,"\n");
		}

		fclose(fp);

		fp = fopen(filename,"a");

		char lineout[nparams*16+35*(ntimes)];

		std::string nameout("");
		nameout += pdata->rdata->SimName;
		nameout += "/";
		nameout += pdata->rdata->output_id;

		fprintf(fp,"%f,%i,%e,%e,"
				   "%i,%i,"
				   "%i,%e,"
				   "%i,%i,%i,"
				   "%i,%i,",
					pdata->dt,pdata->nsteps,(double)nsubsteps,nsubsteps*0,
					pdata->nptcls_species[0]+pdata->nptcls_species[1],pdata->nx,
					status->npiccard_outer,npiccard_std,
					pdata->node_info->cpu_info->vec_length,pdata->node_info->nCPU,pdata->node_info->nTasks_g,
					num_gpus,pdata->rdata->lo_all);

		fprintf(fp,"%s,",nameout.c_str());


		for(int i=0;i<ntimes;i++)
		{
			fprintf(fp,"%f,%f,%f,%f,",run_times[i].x,run_times[i].y,adjusted_times[i].x,adjusted_times[i].y);
		}

		fprintf(fp,"\n");


		fclose(fp);

		printf("\n");
		for(int i=0;i<ntimes;i++)
		{
			char temp[30];
			char temp2[20];

			sprintf(temp,"%s runtime: ",time_names[i]);
			sprintf(temp2,"%*f",(40-strlen(temp)),adjusted_times[i].x);
			printf("%s %s(ns)\n",temp,temp2);
		}

	}

	//output();
}

void OutputCtrl::SaveField(FieldDataCPU* fields,int istep,int iField)
{
	using namespace std;
	std::string* fieldname;
	if(iField == 0)
		fieldname = EfieldFile;
	else if(iField == 1)
		fieldname = BfieldFile;
	else if(iField == 2)
		fieldname = AfieldFile;

	std::string fName;

	fName = std::string(pdata->rdata->outpath) + "/data/" + (*fieldname) + ".dat";

	std::ofstream oStr(fName.c_str(), std::ofstream::app);

	if (!oStr) {
		cout << "Could not open input file " << fName << endl;
		exit(-1);
	}

	oStr << "############################################" << endl;
	oStr << "TIMESTEP: " << istep << endl;
	oStr << "############################################" << endl;


	for(int k=0;k<pdata->nz+1;k++)
	{
		for(int j=0;j<pdata->ny+1;j++)
		{
			std::string line("");
			for(int i=0;i<pdata->nx+1;i++)
			{
				char temp[69];
				realkind Ex,Ey,Ez;

				// Interpolate field values to cell vertices
				if(iField == 0)
				{
					Ex = fields->intrpE(0.0,0.0,0.0,i,j,k,0,FieldData_deriv_f);
					Ey = fields->intrpE(0.0,0.0,0.0,i,j,k,1,FieldData_deriv_f);
					Ez = fields->intrpE(0.0,0.0,0.0,i,j,k,2,FieldData_deriv_f);
				}
				else if(iField == 1)
				{
					Ex = fields->intrpB(0.0,0.0,0.0,i,j,k,0,FieldData_deriv_f);
					Ey = fields->intrpB(0.0,0.0,0.0,i,j,k,1,FieldData_deriv_f);
					Ez = fields->intrpB(0.0,0.0,0.0,i,j,k,2,FieldData_deriv_f);
				}
				else if(iField == 2)
				{
//					Ex = fields->intrpA(0.0,0.0,0.0,i,j,k,0,FieldData_deriv_f);
//					Ey = fields->intrpA(0.0,0.0,0.0,i,j,k,1,FieldData_deriv_f);
//					Ez = fields->intrpA(0.0,0.0,0.0,i,j,k,2,FieldData_deriv_f);
				}

				// Store the tuple in a string
				sprintf(temp,"(%20.15e  %20.15e %20.15e), ",Ex,Ey,Ez);

				// append the string to the line
				line += temp;
			}

			// Write Line to File
			oStr << line << endl;
		}
	}

	oStr.close();

}

void OutputCtrl::SaveMesh(void)
{
	using namespace std;


	std::string fName;

	fName = std::string(pdata->rdata->outpath) + "/data/MeshData.dat";

	std::ofstream oStr(fName.c_str(), std::ofstream::out);

	if (!oStr) {
		cout << "Could not open input file " << fName << endl;
		exit(-1);
	}

	oStr << "############################################" << endl;
	oStr << "(X,Y,Z) VERTEX VALUES" << endl;
	oStr << "############################################" << endl;


	for(int k=0;k<pdata->nz+1;k++)
	{
		for(int j=0;j<pdata->ny+1;j++)
		{
			std::string line;
			for(int i=0;i<pdata->nx+1;i++)
			{
				char temp[69];
				realkind x,y,z;

				x = i*pdata->dxdi + pdata->xmin;
				y = j*pdata->dydi + pdata->ymin;
				z = k*pdata->dzdi + pdata->zmin;

				// Store the tuple in a string
				sprintf(temp,"(%20.15e  %20.15e  %20.15e), ",x,y,z);

				// append the string to the line
				line += temp;
			}

			// Write Line to File
			oStr << line << endl;
		}
	}

	oStr.close();

}


//OutputPage ImplicitPIC::get_subcycle_dists()
//{
//	OutputPage subcycle_dists;
//	/* Subcycle output data */
//	int ndists = pdata->nspecies * pdata->ndevices;
//	double dists[ndists][NSUBCYCLE_BINS+2];
//	double* dist_temp;
//	double null_dist[NSUBCYCLE_BINS+2];
//	int nnodes[ndists];
//
//	double4 my_sub_stats = particles_old[0]->subcycle_stats(pdata);
//	double4 substats[pdata->num_nodes];
//	double tempstats[4] = {my_sub_stats.x,my_sub_stats.y,my_sub_stats.z,my_sub_stats.w};
//	double substats_t[4*pdata->num_nodes];
//	printf("substats: %f, %f, %f, %f\n",tempstats[0],tempstats[1],tempstats[2],tempstats[3]);
//
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	// Gather all of the subcycle statistics
//	MPI_Gather(tempstats,4,MPI_DOUBLE,substats_t,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
//	MPI_Barrier(MPI_COMM_WORLD);
//	if(pdata->mynode == 0)
//	for(int i=0;i<pdata->num_nodes;i++)
//	{
//		substats[i] = make_double4(substats_t[4*i],substats_t[4*i+1],substats_t[4*i+2],substats_t[4*i+3]);
//		printf("substats[%i]: %f, %f, %f, %f\n",i,substats[i].x,substats[i].y,substats[i].z,substats[i].w);
//
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	for(int i=0;i<NSUBCYCLE_BINS+2;i++)
//		null_dist[i] = 0;
//
//	// get the distributions
//	dist_temp = subcycle_dist(pdata,particles_old[0]);
//
//
//	// get reduce everything
//	printf("Reducing Subcycle Distributions\n");
//	for(int i = 0;i<pdata->nspecies;i++)
//		for(int j=0;j<pdata->ndevices;j++)
//		{
//			int node_count = 0;
//			double* dist_reduce = null_dist;
//			if(pdata->my_species == i && pdata->device_type == j)
//			{
//				dist_reduce = dist_temp;
//				node_count = 1;
//			}
//
//			MPI_Reduce(dist_reduce,dists[j + pdata->ndevices*i],NSUBCYCLE_BINS+2,
//					MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//			MPI_Reduce(&node_count,nnodes+j + pdata->ndevices*i,1,
//					MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
//		}
//
//	if(pdata->mynode == 0)
//	{
//		printf("Normalizing Subcycle Distributions\n");
//		// Normalize everything
//		for(int i = 0;i<ndists;i++)
//			for(int k=0;k<NSUBCYCLE_BINS+2;k++)
//				dists[i][k] /= (double)nnodes[i];
//
//		for(int i=0;i<pdata->nspecies;i++)
//		{
//			double mins = dists[pdata->ndevices*i][NSUBCYCLE_BINS];
//			double maxs = dists[pdata->ndevices*i][NSUBCYCLE_BINS+1];
//			for(int j=0;j<pdata->ndevices;j++)
//			{
//				mins = fmin(mins,dists[j+pdata->ndevices*i][NSUBCYCLE_BINS]);
//				maxs = fmax(maxs,dists[j+pdata->ndevices*i][NSUBCYCLE_BINS+1]);
//			}
//
//			for(int j=0;j<pdata->ndevices;j++)
//			{
//				dists[j+pdata->ndevices*i][NSUBCYCLE_BINS] = mins;
//				dists[j+pdata->ndevices*i][NSUBCYCLE_BINS+1] = maxs;
//			}
//		}
//
//		if(pdata->plot_flag){
//		// Plot distributions to gnuplot.
//		gnuplot_ctrl* plots[ndists];
//		for(int i=0;i<ndists;i++)
//		{
//			float tempys[NSUBCYCLE_BINS];
//			float tempxs[NSUBCYCLE_BINS];
//			double mins = dists[i][NSUBCYCLE_BINS];
//			double maxs = dists[i][NSUBCYCLE_BINS+1];
//			double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);
//
//			if(!isnan(dists[i][0])){
//			for(int k=0;k<NSUBCYCLE_BINS;k++)
//			{
//				tempys[k] = dists[i][k];
//				tempxs[k] = k*dsdi+mins;
//			}
//
//			plots[i] = gnuplot_init();
//
//			gnuplot_plot_xy(plots[i],tempxs,tempys,NSUBCYCLE_BINS);
//			}
//
//
//		}
//		}
//
//		char temp[128];
//
//		// Do the Node stats
//		printf("Printing Node Subcycle Stats\n");
//		subcycle_dists.nextline() = "#Node, Average Subs, std, min, max";
//		for(int i=0;i<pdata->num_nodes;i++)
//		{
//			sprintf(temp,"#%i, %f, %f, %f, %f",
//					i,substats[i].x,substats[i].y,substats[i].z,substats[i].w);
//			subcycle_dists.nextline() = temp;
//		}
//
//		// Do the distributions
//		printf("Printing Subcycle Distributions\n");
//		for(int i=0;i<pdata->nspecies;i++)
//		{
//			std::string legend("#NSubcycles");
//			for(int j=0;j<pdata->ndevices;j++)
//			{
//				sprintf(temp,",device:%i",j);
//				legend += temp;
//			}
//			subcycle_dists.nextline() = legend;
//			sprintf(temp,"#ispecies: %i",i);
//			subcycle_dists.nextline() = temp;
//
//			// Do the arrays
//			double mins = dists[pdata->ndevices*i][NSUBCYCLE_BINS];
//			double maxs = dists[pdata->ndevices*i][NSUBCYCLE_BINS+1];
//			double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);
//
//
//			for(int k=0;k<NSUBCYCLE_BINS;k++)
//			{
//				std::string &line = subcycle_dists.nextline();
//				sprintf(temp,"%f",k*dsdi + mins);
//				line = temp;
//
//				for(int j=0;j<pdata->ndevices;j++)
//				{
//					double count = dists[j + pdata->ndevices*i][k];
//
//					sprintf(temp,", %e",count);
//		//			printf("line[%i] = %e\n",k,count);
//					line += temp;
//				}
//
//
//
//			}
//
//		}
//	}
//
//	return subcycle_dists;
//}
//
//OutputPage ImplicitPIC::get_piccard_dists()
//{
//	OutputPage piccard_dists;
//	/* Subcycle output data */
//	int ndists = pdata->nspecies * pdata->ndevices;
//	double dists[ndists][NSUBCYCLE_BINS+2];
//	double* dist_temp;
//	double null_dist[NSUBCYCLE_BINS+2];
//	int nnodes[ndists];
//
//	double4 my_sub_stats = particles_old[0]->piccard_stats(pdata);
//	double4 substats[pdata->num_nodes];
//	double tempstats[4] = {my_sub_stats.x,my_sub_stats.y,my_sub_stats.z,my_sub_stats.w};
//	double substats_t[4*pdata->num_nodes];
//	printf("piccardstats: %f, %f, %f, %f\n",tempstats[0],tempstats[1],tempstats[2],tempstats[3]);
//
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	// Gather all of the subcycle statistics
//	MPI_Gather(tempstats,4,MPI_DOUBLE,substats_t,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
//	MPI_Barrier(MPI_COMM_WORLD);
//	if(pdata->mynode == 0)
//	for(int i=0;i<pdata->num_nodes;i++)
//	{
//		substats[i] = make_double4(substats_t[4*i],substats_t[4*i+1],substats_t[4*i+2],substats_t[4*i+3]);
//		printf("substats[%i]: %f, %f, %f, %f\n",i,substats[i].x,substats[i].y,substats[i].z,substats[i].w);
//
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	for(int i=0;i<NSUBCYCLE_BINS+2;i++)
//		null_dist[i] = 0;
//
//	// get the distributions
//	dist_temp = piccard_dist(pdata,particles_old[0]);
//
//
//	// get reduce everything
//	if(pdata->mynode == 0)printf("Reducing Piccard Distributions\n");
//	for(int i = 0;i<pdata->nspecies;i++)
//		for(int j=0;j<pdata->ndevices;j++)
//		{
//			int node_count = 0;
//			double* dist_reduce = null_dist;
//			if(pdata->my_species == i && pdata->device_type == j)
//			{
//				dist_reduce = dist_temp;
//				node_count = 1;
//			}
//
//			MPI_Reduce(dist_reduce,dists[j + pdata->ndevices*i],NSUBCYCLE_BINS+2,
//					MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//			MPI_Reduce(&node_count,nnodes+j + pdata->ndevices*i,1,
//					MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
//		}
//
//	if(pdata->mynode == 0)
//	{
//		printf("Normalizing Piccard Distributions\n");
//		// Normalize everything
//		for(int i = 0;i<ndists;i++)
//			for(int k=0;k<NSUBCYCLE_BINS+2;k++)
//				dists[i][k] /= (double)nnodes[i];
//
//		for(int i=0;i<pdata->nspecies;i++)
//		{
//			double mins = dists[pdata->ndevices*i][NSUBCYCLE_BINS];
//			double maxs = dists[pdata->ndevices*i][NSUBCYCLE_BINS+1];
//			for(int j=0;j<pdata->ndevices;j++)
//			{
//				mins = fmin(mins,dists[j+pdata->ndevices*i][NSUBCYCLE_BINS]);
//				maxs = fmax(maxs,dists[j+pdata->ndevices*i][NSUBCYCLE_BINS+1]);
//			}
//
//			for(int j=0;j<pdata->ndevices;j++)
//			{
//				dists[j+pdata->ndevices*i][NSUBCYCLE_BINS] = mins;
//				dists[j+pdata->ndevices*i][NSUBCYCLE_BINS+1] = maxs;
//			}
//		}
//
//		// Plot distributions to gnuplot.
//		if(pdata->plot_flag){
//		gnuplot_ctrl* plots[ndists];
//		for(int i=0;i<ndists;i++)
//		{
//			float tempys[NSUBCYCLE_BINS];
//			float tempxs[NSUBCYCLE_BINS];
//			double mins = dists[i][NSUBCYCLE_BINS];
//			double maxs = dists[i][NSUBCYCLE_BINS+1];
//			double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);
//			if(!isnan(dists[i][0])){
//			for(int k=0;k<NSUBCYCLE_BINS;k++)
//			{
//				tempys[k] = dists[i][k];
//				tempxs[k] = k*dsdi+mins;
//			}
//
//			plots[i] = gnuplot_init();
//
//			gnuplot_plot_xy(plots[i],tempxs,tempys,NSUBCYCLE_BINS);
//
//			}
//
//		}
//		}
//
//		char temp[128];
//
//		// Do the Node stats
//		printf("Printing Node Piccard Stats\n");
//		piccard_dists.nextline() = "#Node, Average Iters, std, min, max";
//		for(int i=0;i<pdata->num_nodes;i++)
//		{
//			sprintf(temp,"#%i, %f, %f, %f, %f",
//					i,substats[i].x,substats[i].y,substats[i].z,substats[i].w);
//			piccard_dists.nextline() = temp;
//		}
//
//		// Do the distributions
//		printf("Printing Subcycle Distributions\n");
//		for(int i=0;i<pdata->nspecies;i++)
//		{
//			std::string legend("#NPiccard");
//			for(int j=0;j<pdata->ndevices;j++)
//			{
//				sprintf(temp,",device:%i",j);
//				legend += temp;
//			}
//			piccard_dists.nextline() = legend;
//			sprintf(temp,"#ispecies: %i",i);
//			piccard_dists.nextline() = temp;
//
//			// Do the arrays
//			double mins = dists[pdata->ndevices*i][NSUBCYCLE_BINS];
//			double maxs = dists[pdata->ndevices*i][NSUBCYCLE_BINS+1];
//			double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);
//
//
//			for(int k=0;k<NSUBCYCLE_BINS;k++)
//			{
//				std::string &line = piccard_dists.nextline();
//				sprintf(temp,"%f",k*dsdi + mins);
//				line = temp;
//				for(int j=0;j<pdata->ndevices;j++)
//				{
//					double count = dists[j + pdata->ndevices*i][k];
//
//					sprintf(temp,", %e",count);
//		//			printf("line[%i] = %e\n",k,count);
//					line += temp;
//				}
//
//
//
//			}
//
//		}
//	}
//
//	return piccard_dists;
//}
//
//
//
//
//void ImplicitPIC::output()
//{
//	OutputPage subcycle_page;
//	OutputPage piccard_page;
//
//	subcycle_page = get_subcycle_dists();
//	piccard_page = get_piccard_dists();
//
//
//	// Get subcycle and piccard distributions
//	if(pdata->mynode == 0)
//	{
//		// First setup the output path name
//		printf("Setting Output Path Name\n");
//		char output_path[128];
//		sprintf(output_path,"./output/%s/%s",pdata->SimName,pdata->output_id);
//
//		// Make the output path if it doesn't exist
//		printf("Checking Output Path\n");
//		mkpath(output_path,0777);
//
//		char filename[128];
//		sprintf(filename,"%s/SubcycleDist.dat",output_path);
//
//		printf("Writing Subcycle Distributions\n");
//		subcycle_page.writepage(filename);
//		sprintf(filename,"%s/PiccardDist.dat",output_path);
//		printf("Writing Subcycle Distributions\n");
//		piccard_page.writepage(filename);
//
//	}
//
//
//}
//
//
//
//
//double* subcycle_dist(PlasmaData* pdata,ParticleList* plist)
//{
//	// allocate distribution array
//	double* dist = (double*)malloc((NSUBCYCLE_BINS+2)*sizeof(double));
//	int* num_subcycles_temp = (int*)malloc(plist->nptcls*sizeof(int));
//
//	if(plist->device_type == 1){
//#ifndef NO_CUDA
//	CUDA_SAFE_CALL(cudaMemcpy(num_subcycles_temp,plist->num_subcycles,plist->nptcls*sizeof(int),
//			cudaMemcpyDeviceToHost));
//#endif
//	}else
//		memcpy(num_subcycles_temp,plist->num_subcycles,plist->nptcls*sizeof(int));
//
//	double4 stats = plist->subcycle_stats(pdata);
//	double mins_l = stats.z; // local min
//	double maxs_l = stats.w; // local max
//	double limits[2];
//	double mins,maxs;
//	// Get the global min and max for the species / device combo
//	for(int i = 0;i<pdata->nspecies;i++)
//	{
//			double limits_t[2] = {500,0};
//			if(plist->ispecies == i)
//			{
//				limits_t[0] = fmax(mins_l,0);
//				limits_t[1] = maxs_l;
//			}
//
//			MPI_Allreduce(limits_t,limits,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
//			MPI_Allreduce(limits_t+1,limits+1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
//
//			if(plist->ispecies == i)
//			{
//				mins = limits[0];
//				maxs = limits[1];
//			}
//	}
//
//
//
//	dist[NSUBCYCLE_BINS] = mins; // Minimum number of subcycles
//	dist[NSUBCYCLE_BINS+1] = maxs; // Maximum number of subcycles
//
//	for(int i=0;i<NSUBCYCLE_BINS;i++)
//		dist[i] = 0;
//
//	double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);
//	double dids = 1.0/dsdi;
//
//	double scale = 1.0/((double)plist->nptcls);
//	double scale2 = 1.0/pdata->npiccard_outer;
//
//	for(int i=0;i<plist->nptcls;i++)
//	{
//		double delta = ((num_subcycles_temp[i]*scale2-mins)*dids);
//		int isub = floor(delta);
//		delta = delta - isub;
//		for(int j=0;j<2;j++)
//		{
//			int i_out;
//			double xp;
//			double temp;
//
//			i_out = isub + j;
//			xp = j - delta;
//
//			dist[std::min(std::max(i_out,0),NSUBCYCLE_BINS-1)] += S1_shape(xp)*scale;
//
//		}
//
//	}
//
//	free(num_subcycles_temp);
//
//	return dist;
//
//}
//
//
//
//double* piccard_dist(PlasmaData* pdata,ParticleList* plist)
//{
//	// allocate distribution array
//	double* dist = (double*)malloc((NSUBCYCLE_BINS+2)*sizeof(double));
//	int* num_subcycles_temp = (int*)malloc(plist->nptcls*sizeof(int));
//	double* num_piccard_temp = (double*)malloc(plist->nptcls*sizeof(double));
//
//	if(plist->device_type == 1){
//#ifndef NO_CUDA
//	CUDA_SAFE_CALL(cudaMemcpy(num_subcycles_temp,plist->num_subcycles,plist->nptcls*sizeof(int),
//			cudaMemcpyDeviceToHost));
//	CUDA_SAFE_CALL(cudaMemcpy(num_piccard_temp,plist->num_piccard,plist->nptcls*sizeof(double),
//			cudaMemcpyDeviceToHost));
//#endif
//
//	}else{
//		memcpy(num_subcycles_temp,plist->num_subcycles,plist->nptcls*sizeof(int));
//		memcpy(num_piccard_temp,plist->num_piccard,plist->nptcls*sizeof(double));
//	}
//
//	double4 stats = plist->piccard_stats(pdata);
//	double mins_l = stats.z; // local min
//	double maxs_l = stats.w; // local max
//	double limits[2];
//	double mins,maxs;
//	// Get the global min and max for the species / device combo
//	for(int i = 0;i<pdata->nspecies;i++)
//	{
//			double limits_t[2] = {500,0};
//			if(plist->ispecies == i)
//			{
//				limits_t[0] = fmax(mins_l,0);
//				limits_t[1] = maxs_l;
//			}
//
//			MPI_Allreduce(limits_t,limits,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
//			MPI_Allreduce(limits_t+1,limits+1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
//
//			if(plist->ispecies == i)
//			{
//				mins = limits[0];
//				maxs = limits[1];
//			}
//	}
//
//
//
//	dist[NSUBCYCLE_BINS] = mins; // Minimum number of subcycles
//	dist[NSUBCYCLE_BINS+1] = maxs; // Maximum number of subcycles
//
//	for(int i=0;i<NSUBCYCLE_BINS;i++)
//		dist[i] = 0;
//
//	double dsdi = (maxs - mins)/((double)NSUBCYCLE_BINS);
//	double dids = 1.0/dsdi;
//
//	double scale = 1.0/((double)plist->nptcls);
//	double scale2 = 1.0/pdata->npiccard_outer;
//
//	for(int i=0;i<plist->nptcls;i++)
//	{
//		double delta = ((num_piccard_temp[i]/num_subcycles_temp[i]-mins)*dids);
//		int isub = floor(delta);
//		delta = delta - isub;
//		for(int j=0;j<2;j++)
//		{
//			int i_out;
//			double xp;
//			double temp;
//
//			i_out = isub + j;
//			xp = j - delta;
//
//			dist[std::min(std::max(i_out,0),NSUBCYCLE_BINS-1)] += S1_shape(xp)*scale;
//
//		}
//
//	}
//
//	free(num_subcycles_temp);
//	free(num_piccard_temp);
//
//	return dist;
//
//}
//
//
//
//
//
//
//






