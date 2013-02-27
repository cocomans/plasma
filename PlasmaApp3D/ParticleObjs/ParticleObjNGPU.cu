
#define GPU_CODE // Directive to compile GPU specific optimizations
#include "../ParticleObjNT.inl"


// Dummy Function to force compilation of GPU code

__device__ __noinline__
void GPU_POBJ_DUMMY_FUNCTION(void)
{
	ParticleObjNT<1,1,1,0> particle1;
	ParticleObjNT<1,1,2,0> particle2;
	ParticleObjNT<1,1,3,0> particle3;
	ParticleObjNT<1,2,2,0> particle4;
	ParticleObjNT<1,2,3,0> particle5;
	ParticleObjNT<1,3,3,0> particle6;

	ParticleObjNT<1,1,1,1> particle7;
	ParticleObjNT<1,1,2,1> particle8;
	ParticleObjNT<1,1,3,1> particle9;
	ParticleObjNT<1,2,2,1> particle10;
	ParticleObjNT<1,2,3,1> particle11;
	ParticleObjNT<1,3,3,1> particle12;

	PlasmaData* pdata;
	FieldData* fields;
	CurrentTally* current;
	typevecN<int,1> iter;
	int nsubcycle_max;

	particle1.push(pdata,fields,current,iter,nsubcycle_max);
	particle2.push(pdata,fields,current,iter,nsubcycle_max);
	particle3.push(pdata,fields,current,iter,nsubcycle_max);
	particle4.push(pdata,fields,current,iter,nsubcycle_max);
	particle5.push(pdata,fields,current,iter,nsubcycle_max);
	particle6.push(pdata,fields,current,iter,nsubcycle_max);

	particle7.push(pdata,fields,current,iter,nsubcycle_max);
	particle8.push(pdata,fields,current,iter,nsubcycle_max);
	particle9.push(pdata,fields,current,iter,nsubcycle_max);
	particle10.push(pdata,fields,current,iter,nsubcycle_max);
	particle11.push(pdata,fields,current,iter,nsubcycle_max);
	particle12.push(pdata,fields,current,iter,nsubcycle_max);
}
