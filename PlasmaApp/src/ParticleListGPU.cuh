#ifndef PARTICLE_LIST_GPU_H
#define PARTICLE_LIST_GPU_H



// Basically this lets us use a CPU particle list
// in the cases where we want to build without CUDA support
// It also lets us build with different versions of ParticleListGPU without having to
// go to virtual functions

#ifndef NO_CUDA
#include "ParticleListGPUSorted.cuh"
typedef ParticleListGPUSorted ParticleListGPU;
#else
#include "ParticleListCPU.h"
typedef ParticleListCPU ParticleListGPU;
#endif

#endif /* PARTICLE_LIST_GPU_H */
