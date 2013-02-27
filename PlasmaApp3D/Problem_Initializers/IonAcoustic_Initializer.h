#ifndef IONACOUSTIC_INITIALIZER_H
#define IONACOUSTIC_INITIALIZER_H

#include "../ProblemInitializer.h"

class PlasmaData;
class ParticleList;
class HOMoments;
class FieldData;


class IonAcoustic_Initializer : public ProblemInitializer
{
public:

	IonAcoustic_Initializer(PlasmaData* pdata_in);

	~IonAcoustic_Initializer(){};

	void initialize_particles(ParticleList* particles, HOMoments* moments,ParallelInfo* myinfo);

	void initialize_fields(FieldData** fields, HOMoments** moments,ParallelInfo* myinfo);

	void init_particle(realkind& px, realkind& py, realkind& pz, int& ix, int& iy, int& iz,
								realkind& vx,realkind& vy,realkind& vz,int ispecies,int iptcl);


	void init_velocities(realkind& vx, realkind& vy, realkind& vz,
						int ispecies,int iptcl);

	void check_step(ParticleList* particles,
			HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo);




	float Ey,Bz;
	double E0,kE0;
	float alpha;
	float kt;
	float v_shift;
	gnuplot_ctrl* plot;
	gnuplot_ctrl* kenergy_plot;
	gnuplot_ctrl* Tenergy_plot;
	gnuplot_ctrl* charge_cons_plot;
	float* Energy_array;
	float* kEnergy_array;
	float* TEnergy_array;
	float* max_array;
	float* max_t_array;
	float* charge_cons_array;
	int nplot;
	int nmax_array;

};

#endif /* IONACOUSTIC_INITIALIZER_H */
