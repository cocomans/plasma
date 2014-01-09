/*-------------------------------------------------------------------------*/
/**
	@file	OutputCtrl.h
	@author	J. Payne
	@date		04/22/2013
	@brief	Declares the RunData class. Also defines simulation precision,
	several macro calls, constants, and utility functions.


*/
/*--------------------------------------------------------------------------*/

#ifndef OUTPUTCTRL_H_
#define OUTPUTCTRL_H_

#include "Util/gnuplot_i.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

class CPUTimer;
class PlasmaData;
class RunStatus;
class NodeFieldData;
class NodeHOMoments;
class FieldDataCPU;
class NodeParticleList;
class HOMomentsCPU;

class LogStream;
class googleChart;
class sseDataStream;
class sseServer;
class WebControlSystem;
class googleChart;


class OutputCtrl
{
public:

	OutputCtrl(PlasmaData* _pdata,int argc, char* argv[]);

	void InitialOutput(NodeFieldData* fields,FieldDataCPU* fields_old,
			FieldDataCPU* fields_next,
			NodeHOMoments* moments,NodeParticleList* particles);

	void StepOutput(NodeHOMoments* moments,
			NodeParticleList* particles,FieldDataCPU* fields_old,FieldDataCPU* fields_next);

	void FinalOutput(NodeHOMoments* moments,
			NodeParticleList* particles,FieldDataCPU* fields_old,FieldDataCPU* fields_next);

	void InitialCoprocess(
			NodeFieldData* fields,
			NodeHOMoments* moments,
			NodeParticleList* particles);

	void StepCoprocess(int istep,
			NodeHOMoments* moments,
			NodeParticleList* particles,FieldDataCPU* fields_old,FieldDataCPU* fields_next);

	void FinalCoprocess(NodeHOMoments* moments,
			NodeParticleList* particles,FieldDataCPU* fields_old,FieldDataCPU* fields_next);

	void PlotFields1D(FieldDataCPU* fields,int icomponent);

	void PlotMoments1D(HOMomentsCPU* moments);

	void PlotParticles1D(NodeParticleList* particles);

	void PlotFields2D(FieldDataCPU* fields,int icomponent);

	void PlotMoments2D(HOMomentsCPU* moments,int ispecies);

	void PlotParticles2D(NodeParticleList* particles);

	void save_plots(void);

	std::string plotpath(std::string plot_name,int plot_number);

	void save_timing(NodeParticleList* particles);

	std::string calc_timing(double perf_weight,const char* legend);

	double total_substeps(void);

	const int& istep(void);

	void SaveField(FieldDataCPU* fields,int istep,int iField);

	void SaveMesh(void);

	void RecordMomentDiffs(double* charge_dif, double* current_dif);

	PlasmaData* pdata;
	RunStatus* status;

	bool iplot;
	bool iSave_plots;
	bool iWrite_data;
	bool iFinalOutput;

	bool iExitSave;

	int OutputInterval;
	int PlotInterval;
	int PlotSaveInterval;

	int CurPlot;


	float* charge_cons;
	float* energy_cons;
	float* field_energy;
	float* particle_energy;
	float* Bfield_energy;
	float* Efield_energy;
	float** species_energy;
	float* momentum_cons;

	float* time_array;

	float** charge_diffs;
	float** current_diffs;
	float* lotoho_array;
	int ncalls_lotoho;

	double kE0,pE0,tE0;

	std::string* EfieldFile;
	std::string* BfieldFile;
	std::string* AfieldFile;



	// Plot Handles
	gnuplot_ctrl** cons_plots;
	gnuplot_ctrl** field_plots;
	gnuplot_ctrl** moment_plots;
	gnuplot_ctrl** particle_plot;
	std::string*** plot_titles;
	std::string*** x_labels;
	std::string*** y_labels;
	int n_cons_plots;
	int n_field_plots;
	int n_moment_plots;
	int n_particle_plots;


	LogStream** lStreams;
	googleChart** charts;
	int nStreams;
	int nCharts;

	WebControlSystem* control;
	sseDataStream* stream;
	sseServer* server;
	int exit_control;

	/// Timer for particle push
	CPUTimer* push_timer;
	/// Timer for entire HO system
	CPUTimer* HOSolve_timer;
	/// Timer for MPI Communication
	CPUTimer* Comm_timer;
	/// Timer for total run time
	CPUTimer* total_timer;
	/// Timer for LO system solve
	CPUTimer* LOSolve_timer;
	/// Timer for single time step
	CPUTimer* step_timer;
	/// Timer for coprocessing
	CPUTimer* coprocess_timer;




};



#endif /* OUTPUTCTRL_H_ */
