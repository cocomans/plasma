#ifndef __ImplicitPICAdaptor_h
#define __ImplicitPICAdaptor_h

class FieldDataCPU;
class ParticleListCPU;
class HOMomentsCPU;

#include <string>

#ifdef __cplusplus
extern "C" {
#endif

using namespace std;
        
// Coprocessor initialize (calls the Catalyst default)
void coprocessorInitialize(string pythonName);

// Create the VTK image data one time
void coprocessorCreateGrid(
                int ndim, int nvel, 
                int nsteps, double dt,
                int nparticle, int nspecies,
                int nx, int ny, int nz,
                double xmin, double ymin, double zmin,
                double dxdi, double dydi, double dzdi,
                float* chargeCons,
                float* energyCons,
                float* fieldEnergy,
                float* particleEnergy);

// Create the subblocks
void createParticleBlock(int specie);
void createFieldBlock();
void createTimeStepData(
                float* chargeCons,
                float* energyCons,
                float* fieldEnergy,
                float* particleEnergy);

// Fill in variable data per time step
void coprocessorFillParticles(
                ParticleListCPU* particles);
void coprocessorFillField(
                FieldDataCPU* field,
                HOMomentsCPU* moments);
void setRange(int var);

// Process the pipeline using the variable data
void coprocessorProcess(int timestep);

// Coprocessor finalize (calls the Catalyst default)
void coprocessorFinalize();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __ImplicitPICAdaptor_h */
