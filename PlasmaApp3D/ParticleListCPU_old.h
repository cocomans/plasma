#ifndef PARTICLE_LIST_CPU_H
#define PARTICLE_LIST_CPU_H

#include "ParticleList.h"

class FieldData;
class HOMoments;
class ParticleObj;


//static int ParticleList_nfloats = 7;
//static int ParticleList_nints = 3;

class ParticleListCPU_old : public ParticleList
{
public:

	ParticleListCPU_old();

	~ParticleListCPU_old();

	void allocate(int nptcls_in);

	void init(ProblemInitializer* initializer, HOMoments* moments);

	realkind evaluate_energy(PlasmaData* pdata);

	long long int push(PlasmaData* pdata, FieldData* fields, HOMoments* moments);

	template<int VEC_LENGTH>
	long long int pushT(PlasmaData* pdata, FieldData* fields, HOMoments* moments);

	template<int VEC_LENGTH>
	long long int push_interface(PlasmaData* pdata, FieldData* fields, HOMoments* moments);

	void copy_from(const ParticleList* list_in);

	realkind& get_fvalue(int iptcl,int ifloat){return *((*get_float(ifloat)) + iptcl);};

	int& get_ivalue(int iptcl,int iint){return *((*get_int(iint)) + iptcl);};

	const realkind& get_fvalue(int iptcl,int ifloat)const{return *((*get_float(ifloat)) + iptcl);};

	const int& get_ivalue(int iptcl,int iint)const{return *((*get_int(iint)) + iptcl);};

	void CPUfree();




	//void accumulate_charge(PlasmaData* pdata, HOMoments* moments);

	// Diagnostic Methods
	void plot_particles(PlasmaData* pdata);





};





#endif /* PARTICLE_LIST_CPU_H */
