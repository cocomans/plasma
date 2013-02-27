#ifndef ExB_INITIALIZER_H
#define ExB_INITIALIZER_H

#include "../ProblemInitializer.h"

class PlasmaData;
class ParticleList;
class HOMoments;
class FieldData;


class ExB_Initializer : public ProblemInitializer
{
/** @brief
 * This class is used to test the accuracy of the particle pusher
 * by comparing it to an analytic solution to a linear varying E-field
 * and constant B field
 */
public:

	ExB_Initializer(PlasmaData* pdata_in);

	~ExB_Initializer(){};

	void initialize_particles(ParticleList* particles, HOMoments* moments,ParallelInfo* myinfo);

	void initialize_fields(FieldData** fields, HOMoments** moments,ParallelInfo* myinfo);

	void init_particle(realkind& px, realkind& py, realkind& pz, int& ix, int& iy, int& iz,
								realkind& vx,realkind& vy,realkind& vz,int ispecies,int iptcl);


	void init_velocities(realkind& vx, realkind& vy, realkind& vz,
						int ispecies,int iptcl);

	void check_step(ParticleList* particles,
			HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo);

	void get_analytic_pos(realkind &xout,realkind &yout,
										   realkind &vxout,realkind &vyout,
										   int iptcl);

	void finish(ParticleList* particles,
			HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo);

	realkind Ex0,Ex1,dEdx;
	realkind Bz;
	ParticleList* particles_temp;
	realkind Omega,omegac;
	realkind* vx0;
	realkind* vy0;
	realkind* x0;
	realkind* y0;
	float* yt;
	float* xt;
	float* xrt;
	float* yrt;
	realkind curtime;
	int iter;
	int nsteps;

};

#endif /* ExB_INITIALIZER_H */
