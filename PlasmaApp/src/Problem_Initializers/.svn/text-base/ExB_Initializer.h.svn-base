#ifndef ExB_INITIALIZER_H
#define ExB_INITIALIZER_H

#include "../ProblemInitializer.h"

class PlasmaData;
class ParticleListCPU;
class HOMomentsCPU;
class FieldDataCPU;
class NodeHOMoments;
class NodeParticleList;
class NodeFieldData;


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

	void initialize_particles(NodeParticleList* particles, NodeHOMoments* moments);

	void initialize_fields(FieldDataCPU* fields, NodeHOMoments* moments);

	void init_particle(realkind& px, realkind& py, realkind& pz, int& ix, int& iy, int& iz,
								realkind& vx,realkind& vy,realkind& vz,int ispecies,int iptcl,
								int nptcls,int ioffset);

	void init_velocities(realkind& vx, realkind& vy, realkind& vz,
								int ispecies,int iptcl);

	void check_step(NodeParticleList* particles_next,
			NodeHOMoments* moments,FieldDataCPU* fields_old,FieldDataCPU* fields_next);

	void finish(NodeParticleList* particles,
			NodeHOMoments* moments,NodeFieldData* fields_half);

	void get_analytic_pos(realkind &xout,realkind &yout,
										   realkind &vxout,realkind &vyout,
										   int iptcl);

	realkind Ex0,Ex1,dEdx;
	realkind Bz;
	ParticleListCPU* particles_temp;
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
