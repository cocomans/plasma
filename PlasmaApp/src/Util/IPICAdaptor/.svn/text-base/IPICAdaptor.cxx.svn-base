#include "IPICAdaptor.h"

#include "CAdaptorAPI.h"
#include "CPythonAdaptorAPI.h"

#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPAdaptorAPI.h"
#include "vtkCPPythonAdaptorAPI.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkImageData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMPI.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIController.h"
#include "vtkPointData.h"
#include "vtkPlot.h"
#include "vtkUnstructuredGrid.h"

#include <float.h>
#include <sstream>

#define NO_CUDA 1
#include "FieldDataCPU.h"
#include "ParticleListCPU.h"
#include "HOMomentsCPU.h"

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

namespace 
{
  // Types of data variables
  const int SCALAR = 1;
  const int VECTOR = 3;
  const int TENSOR = 9;

  // VTK Multiblock block types
  const int FIELD_BLOCK = 0;
  const int TIME_BLOCK = 1;
  const int PARTICLE_BLOCK = 2;

  // Pointers to actual VTK data storage
  float* eP;            // Field Electric (vector)
  float* bP;            // Field Magnetic (vector)
  float** chargeP;      // Moment Charge by species (scalar)
  float** currentP;     // Moment Current by species (vector)
  float** stressP;      // Moment Stress by species (tensor)

  int mpi_rank;
  int mpi_totalRank;

  int numberOfDim;
  int numberOfVel;
  int numberOfSteps;
  int numberOfParticles;
  int numberOfSpecies;
  int numberOfTuples;

  int numberOfParticlesMax;

  int dimension[VECTOR];
  double origin[VECTOR];
  double delta[VECTOR];
  double timeDelta;
}

using namespace std;

//////////////////////////////////////////////////////////////////////////
//
// Coprocessor initialize (calls default but saves extra declaration)
//
//////////////////////////////////////////////////////////////////////////

void coprocessorInitialize(string pythonName)
{
  vtkCPPythonAdaptorAPI::CoProcessorInitialize(pythonName.c_str());
}

//////////////////////////////////////////////////////////////////////////
//
// Create the coprocessing grid one time
//
//////////////////////////////////////////////////////////////////////////

void coprocessorCreateGrid(
                         int ndim, int nvel, 
                         int nsteps, double dt,
                         int nparticles, int nspecies,
                         int nx, int ny, int nz,
                         double xmin, double ymin, double zmin,
                         double dxdi, double dydi, double dzdi,
                         float* chargeCons,
                         float* energyCons,
                         float* fieldEnergy,
                         float* particleEnergy)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_totalRank);

  // Collect size parameters for grids based on type of problem
  numberOfDim = ndim;
  numberOfVel = nvel;
  numberOfSteps = nsteps;
  numberOfParticles = nparticles;
  numberOfSpecies = nspecies;
  timeDelta = dt;

  numberOfParticlesMax = 32000;

  origin[0] = xmin;  origin[1] = ymin;  origin[2] = zmin;
  delta[0] = dxdi;   delta[1] = dydi;   delta[2] = dzdi;

  if (numberOfDim == 1) {
    dimension[0] = nx;
    dimension[1] = 1;
    dimension[2] = 1;
  } else if (numberOfDim == 2) {
    dimension[0] = nx;
    dimension[1] = ny;
    dimension[2] = 1;
  }
  numberOfTuples = dimension[0] * dimension[1] * dimension[2];

  // Create enclosing multiblock data set and add to coprocessor
  vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::New();
  vtkCPAdaptorAPI::GetCoProcessorData()->
                   GetInputDescriptionByName("input")->SetGrid(grid);
  grid->SetNumberOfBlocks(numberOfSpecies + 2);

  // Create a block for field data which holds field and momentum
  createFieldBlock();

  // Create a block for time step data which is in FieldData of multiblock
  // but which also will be a vtkTable block for python pipeline
  createTimeStepData(chargeCons, energyCons, 
                     fieldEnergy, particleEnergy);

  // Create a block for each particle specie
  for (int specie = 0; specie < numberOfSpecies; specie++) {
    createParticleBlock(specie);
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Create the time plot arrays attached to the multiblock
//
/////////////////////////////////////////////////////////////////////////

void createTimeStepData(float* chargeCons,
                        float* energyCons,
                        float* fieldEnergy,
                        float* particleEnergy)
{
  vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::SafeDownCast(
      vtkCPAdaptorAPI::GetCoProcessorData()->
                       GetInputDescriptionByName ("input")->GetGrid ());

  vtkImageData* tgrid = vtkImageData::New();
  tgrid->SetDimensions(numberOfSteps, 1, 1);
  tgrid->SetOrigin(0, 0, 0);
  tgrid->SetSpacing(timeDelta, 1, 1);
  grid->SetBlock(TIME_BLOCK, tgrid);

  // Pass in the actual array from simulation so don't delete
  int saveArray = 1;

  vtkFloatArray* ChargeCons = vtkFloatArray::New();
  ChargeCons->SetName("Charge Conservation");
  ChargeCons->SetNumberOfComponents(1);
  ChargeCons->SetNumberOfTuples(numberOfSteps);
  ChargeCons->SetArray(chargeCons, numberOfSteps, saveArray);

  vtkFloatArray* EnergyCons = vtkFloatArray::New();
  EnergyCons->SetName("Energy Conservation");
  EnergyCons->SetNumberOfComponents(1);
  EnergyCons->SetNumberOfTuples(numberOfSteps);
  EnergyCons->SetArray(energyCons, numberOfSteps, saveArray);

  vtkFloatArray* FieldEnergy = vtkFloatArray::New();
  FieldEnergy->SetName("Field Energy");
  FieldEnergy->SetNumberOfComponents(1);
  FieldEnergy->SetNumberOfTuples(numberOfSteps);
  FieldEnergy->SetArray(fieldEnergy, numberOfSteps, saveArray);

  vtkFloatArray* ParticleEnergy = vtkFloatArray::New();
  ParticleEnergy->SetName("Particle Energy");
  ParticleEnergy->SetNumberOfComponents(1);
  ParticleEnergy->SetNumberOfTuples(numberOfSteps);
  ParticleEnergy->SetArray(particleEnergy, numberOfSteps, saveArray);

  tgrid->GetPointData()->AddArray(ChargeCons);
  tgrid->GetPointData()->AddArray(EnergyCons);
  tgrid->GetPointData()->AddArray(FieldEnergy);
  tgrid->GetPointData()->AddArray(ParticleEnergy);

  ChargeCons->Delete();
  EnergyCons->Delete();
  FieldEnergy->Delete();
  ParticleEnergy->Delete();
}

/////////////////////////////////////////////////////////////////////////
//
// Create unstructured grid which will contain Particles
//
/////////////////////////////////////////////////////////////////////////

void createParticleBlock(int specie)
{
  vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::SafeDownCast(
      vtkCPAdaptorAPI::GetCoProcessorData()->
                       GetInputDescriptionByName ("input")->GetGrid ());

  vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
  ugrid->Initialize();

  // Allocate storage for Particles
  vtkPoints* points = vtkPoints::New();
  points->SetDataTypeToFloat();
  ugrid->SetPoints(points);
  ugrid->Allocate(numberOfParticlesMax, numberOfParticlesMax);



  vtkFloatArray* vel = vtkFloatArray::New();
  std::stringstream varName;
  varName << "Velocity_" << specie;
  vel->SetName(varName.str().c_str());
  vel->SetNumberOfComponents(VECTOR);
  vel->Allocate(numberOfParticlesMax);

  float value[VECTOR];
  for (int i = 0; i < VECTOR; i++)
    value[i] = 0.0;

  // Points are always 3D
  for (int i = 0; i < numberOfParticlesMax; i++) {
    vtkIdType id = points->InsertNextPoint(value);
    ugrid->InsertNextCell(VTK_VERTEX, 1, &id);
    vel->InsertNextTuple(value);
  }

  ugrid->GetCellData()->AddArray(vel);
  grid->SetBlock(PARTICLE_BLOCK + specie, ugrid);

  points->Delete();
  vel->Delete();
  ugrid->Delete();
}
  
/////////////////////////////////////////////////////////////////////////
//
// Create ImageData (uniform rectilinear grid) for Field and Moment data
// Moment data is per species
//
/////////////////////////////////////////////////////////////////////////

void createFieldBlock()
{
  vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::SafeDownCast(
      vtkCPAdaptorAPI::GetCoProcessorData()->
                       GetInputDescriptionByName ("input")->GetGrid ());

  vtkImageData* sgrid = vtkImageData::New();
  sgrid->SetDimensions(dimension[0] + 1, dimension[1] + 1, dimension[2] + 1);
  sgrid->SetOrigin(origin);
  sgrid->SetSpacing(delta);
  grid->SetBlock(FIELD_BLOCK, sgrid);

  float svalue = 0;
  float vvalue[VECTOR] = {0, 0, 0};
  float tvalue[TENSOR] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  // Create arrays to hold Field data
  vtkFloatArray* eArray = vtkFloatArray::New();
  eArray->SetName("Electric");
  eArray->SetNumberOfComponents(VECTOR);
  eArray->Allocate(numberOfTuples);

  vtkFloatArray* bArray = vtkFloatArray::New();
  bArray->SetName("Magnetic");
  bArray->SetNumberOfComponents(VECTOR);
  bArray->Allocate(numberOfTuples);

  for (int i = 0; i < numberOfTuples; i++) {
    eArray->InsertNextTuple(vvalue);
    bArray->InsertNextTuple(vvalue);
  }

  eP = eArray->GetPointer(0);
  bP = bArray->GetPointer(0);
  sgrid->GetCellData()->AddArray(eArray);
  sgrid->GetCellData()->AddArray(bArray);
  eArray->Delete();
  bArray->Delete();

  // Create arrays to hold Moment data by species
  chargeP = new float*[numberOfSpecies]; 
  currentP = new float*[numberOfSpecies]; 
  stressP = new float*[numberOfSpecies]; 

  for (int specie = 0; specie < numberOfSpecies; specie++) {
    vtkFloatArray* chargeArray = vtkFloatArray::New();
    ostringstream chargeName;
    chargeName << "Charge" << specie;
    chargeArray->SetName(chargeName.str().c_str());
    chargeArray->SetNumberOfComponents(SCALAR);
    chargeArray->Allocate(numberOfTuples);

    vtkFloatArray* currentArray = vtkFloatArray::New();
    ostringstream currentName;
    currentName << "Current" << specie;
    currentArray->SetName(currentName.str().c_str());
    currentArray->SetNumberOfComponents(VECTOR);
    currentArray->Allocate(numberOfTuples);

    vtkFloatArray* stressArray = vtkFloatArray::New();
    ostringstream stressName;
    stressName << "Stress" << specie;
    stressArray->SetName(stressName.str().c_str());
    stressArray->SetNumberOfComponents(TENSOR);
    stressArray->Allocate(numberOfTuples);

    for (int i = 0; i < numberOfTuples; i++) {
      chargeArray->InsertNextValue(svalue);
      currentArray->InsertNextTuple(vvalue);
      stressArray->InsertNextTuple(tvalue);
    }

    chargeP[specie] = chargeArray->GetPointer(0);
    currentP[specie] = currentArray->GetPointer(0);
    stressP[specie] = stressArray->GetPointer(0);
    sgrid->GetCellData()->AddArray(chargeArray);
    sgrid->GetCellData()->AddArray(currentArray);
    sgrid->GetCellData()->AddArray(stressArray);
    chargeArray->Delete();
    currentArray->Delete();
    stressArray->Delete();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Fill the particle structure for this time step
//
//////////////////////////////////////////////////////////////////////////

void coprocessorFillParticles(ParticleListCPU* particleLists)
{
  if (!vtkCPAdaptorAPI::GetCoProcessor() || 
      !vtkCPAdaptorAPI::GetCoProcessorData()) {
    cerr << "CoProcessor has not been properly initialized" << endl;
    return;
  }

  vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::SafeDownCast(
      vtkCPAdaptorAPI::GetCoProcessorData()->
                       GetInputDescriptionByName ("input")->GetGrid ());

  float point[VECTOR], value[VECTOR];
  for (int specie = 0; specie < numberOfSpecies; specie++) {
    int indx = 0;
    ParticleListCPU* particles = particleLists + specie;

    vtkUnstructuredGrid* ugrid = 
      vtkUnstructuredGrid::SafeDownCast(grid->GetBlock(PARTICLE_BLOCK+specie));
    vtkPoints* points = ugrid->GetPoints();

    std::stringstream varName;
    varName << "Velocity_" << specie;
    vtkFloatArray* vel =  vtkFloatArray::SafeDownCast(
        ugrid->GetCellData()->GetArray(varName.str().c_str()));

    int stride = (numberOfParticles + numberOfParticlesMax -1)/numberOfParticlesMax;
    if (numberOfDim == 1) {
      for (int i = 0; i < numberOfParticles; i+=stride) {
        point[0] = origin[0] + delta[0] * (particles->px[i] + particles->ix[i]);
        point[1] = particles->vx[i];
        point[2] = 0.0;
        points->SetPoint(indx, point);

        value[0] = particles->vx[i];
        value[1] = particles->vy[i];
        value[2] = particles->vz[i];
        vel->SetTupleValue(indx, value);
        indx++;
      }
    }
    else {
      for (int i = 0; i < numberOfParticles; i+=stride) {
        point[0] = origin[0] + delta[0] * (particles->px[i] + particles->ix[i]);
        point[1] = origin[1] + delta[1] * (particles->py[i] + particles->iy[i]);
        point[2] = 0.0;
        points->SetPoint(indx, point);
  
        value[0] = particles->vx[i];
        value[1] = particles->vy[i];
        value[2] = particles->vz[i];
        vel->SetTupleValue(indx, value);
        indx++;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Fill the field data arrays for the grid for this time step
// Explicitly fill each variable because access is very specific
// If IPIC went to float* arrays this would be easier
// VTK data is stored with components of a tuple being contiguous
//
//////////////////////////////////////////////////////////////////////////

void coprocessorFillField(FieldDataCPU* fields, HOMomentsCPU* moments)
{
  if (!vtkCPAdaptorAPI::GetCoProcessor() || 
      !vtkCPAdaptorAPI::GetCoProcessorData()) {
    cerr << "CoProcessor has not been properly initialized" << endl;
    return;
  }

  vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::SafeDownCast(
      vtkCPAdaptorAPI::GetCoProcessorData()->
                       GetInputDescriptionByName ("input")->GetGrid ());

  // Iterate over data structure extracting data and inserting into VTK
  int indx = 0;
  for (int k = 0; k < dimension[2]; k++) {
    for (int j = 0; j < dimension[1]; j++) {
      for (int i = 0; i < dimension[0]; i++) {
        for (int comp = 0; comp < VECTOR; comp++) {
          eP[indx] = fields->getE(i, j, k, comp);
          bP[indx] = fields->getB(i, j, k, comp);
          indx++;
        }
      }
    }
  }

  // Iterate over data structure extracting data and inserting into VTK
  for (int specie = 0; specie < numberOfSpecies; specie++) {
    int charge = 0;
    int current = 0;
    int stress = 0;
    for (int k = 0; k < dimension[2]; k++) {
      for (int j = 0; j < dimension[1]; j++) {
        for (int i = 0; i < dimension[0]; i++) {

          chargeP[specie][charge++] = 
                  moments->get_val(i, j, k, specie, HOMoments_charge);

          currentP[specie][current++] = 
                  moments->get_val(i, j, k, specie, HOMoments_currentx);
          currentP[specie][current++] = 
                  moments->get_val(i, j, k, specie, HOMoments_currenty);
          currentP[specie][current++] = 
                  moments->get_val(i, j, k, specie, HOMoments_currentz);

          stressP[specie][stress++] = 
                  moments->get_val(i, j, k, specie, HOMoments_S2xx);
          stressP[specie][stress++] = 
                  moments->get_val(i, j, k, specie, HOMoments_S2xy);
          stressP[specie][stress++] = 
                  moments->get_val(i, j, k, specie, HOMoments_S2xz);
          stressP[specie][stress++] = 
                  moments->get_val(i, j, k, specie, HOMoments_S2xy);
          stressP[specie][stress++] = 
                  moments->get_val(i, j, k, specie, HOMoments_S2yy);
          stressP[specie][stress++] = 
                  moments->get_val(i, j, k, specie, HOMoments_S2yz);
          stressP[specie][stress++] = 
                  moments->get_val(i, j, k, specie, HOMoments_S2xz);
          stressP[specie][stress++] = 
                  moments->get_val(i, j, k, specie, HOMoments_S2yz);
          stressP[specie][stress++] = 
                  moments->get_val(i, j, k, specie, HOMoments_S2zz);
        }
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////
//
// Coprocessing for the timestep
// VTK structure has already been filled with data
// CoProcess() is called on the added pipeline which is hardcoded or python
//
//////////////////////////////////////////////////////////////////////////

void coprocessorProcess(int timestep)
{
  if (!vtkCPAdaptorAPI::GetCoProcessor() || 
      !vtkCPAdaptorAPI::GetCoProcessorData()) {
    cerr << "CoProcessor has not been properly initialized" << endl;
    return;
  }

  if (mpi_rank == 0) cout << "COPROCESS " << timestep << endl;

  vtkMultiBlockDataSet* grid = vtkMultiBlockDataSet::SafeDownCast(
      vtkCPAdaptorAPI::GetCoProcessorData()->
                       GetInputDescriptionByName ("input")->GetGrid ());

  double time = timestep * timeDelta;
  vtkCPAdaptorAPI::GetCoProcessorData()->SetTimeData (time, timestep);

  if (vtkCPAdaptorAPI::GetCoProcessor()->
         RequestDataDescription(vtkCPAdaptorAPI::GetCoProcessorData())) {
    vtkCPAdaptorAPI::GetCoProcessor()->CoProcess(
                                         vtkCPAdaptorAPI::GetCoProcessorData());
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Coprocessor finalize (calls default from Catalyst)
//
//////////////////////////////////////////////////////////////////////////

void coprocessorFinalize()
{
  vtkCPAdaptorAPI::CoProcessorFinalize();
}
