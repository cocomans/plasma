

CC=gcc#bf-gcc
CXX=mpicxx# -cxx=bf-g++ -bf-verbose -bf-thread-safe -bf-all-ops -bf-by-func 
NVCC=nvcc

#CC=icc
#CXX=mpicxx# -cxx=bf-g++ -bf-verbose -bf-thread-safe -bf-by-func -bf-all-ops

TRILINOS_DIR=/usr/local/trilinos


#BYFL_INCLUDE="_ZN12CurrentTally9tally1d1vEddid,_ZN12FieldDataCPU5getETILi0EEERdiii,_ZN12FieldDataCPU6intrpEEdddiiii15FieldData_deriv,_ZN12FieldDataCPU7intrpETILi0EL15FieldData_deriv0EEEddddiii,_ZN11ChargeTally7tally1dEdid,_ZN13ParticleObjNTILi1ELi1ELi1ELb0EE13PicardIterateEP10PlasmaDataP9FieldDataP12CurrentTallyR8typevecNIdLi1EE,_ZN13ParticleObjNTILi1ELi1ELi1ELb0EE13estimate_dtauEP10PlasmaDataP9FieldData,_ZN13ParticleObjNTILi1ELi1ELi1ELb0EE18accumulate_currentEP10PlasmaDataP12CurrentTally8typevecNIS5_IdLi1EELi1EES7_RKS6_,_ZN13ParticleObjNTILi1ELi1ELi1ELb0EE18time_till_crossingER8typevecNIdLi1EES3_S3_d,_ZN13ParticleObjNTILi1ELi1ELi1ELb0EE4pushEP10PlasmaDataP9FieldDataP12CurrentTallyR8typevecNIiLi1EEi,_ZN13ParticleObjNTILi1ELi1ELi1ELb0EE5get_EIL15FieldData_deriv0EEE8typevecNIS3_IdLi1EELi1EEP9FieldDataRKS5_RKS3_IS3_IiLi1EELi1EE,_ZN13ParticleObjNTILi1ELi1ELi1ELb0EE7copy_inERK12ParticleListi,_ZN13ParticleObjNTILi1ELi1ELi1ELb0EE9get_accelEP9FieldDataRK8typevecNIS3_IdLi1EELi1EES7_,_ZN15ParticleListCPU5pushTILi1ELi1ELi1ELb0EEExP10PlasmaDataP9FieldDataP9HOMoments._omp_fn.0,_ZN8typevecNIS_IdLi1EELi1EEmlIdEES1_RKT_,_ZNK8typevecNIS_IdLi1EELi1EEplIS0_EES1_RKS_IT_Li1EE"
BYFL_INCLUDE="_ZN15ParticleListCPU5pushTILi1ELi2ELi3ELb1EEExP10PlasmaDataP9FieldDataP9HOMoments._omp_fn.9,_ZN11ChargeTally5tallyEdddiiid,_ZN11ChargeTally7tally2dEddiid,_ZN11ChargeTallyC2EPd4int3dddi,_ZN11StressTally5tallyEddddddiiid,_ZN11StressTally9tally1d1vEddid,_ZN12CurrentTally5tallyEddddddiiid,_ZN12CurrentTally9tally2d3vEdddddiid,_ZN13ParticleObjNTILi1ELi2ELi3ELb1EE13PicardIterateEP10PlasmaDataP9FieldDataP12CurrentTallyR8typevecNIdLi1EE,_ZN13ParticleObjNTILi1ELi2ELi3ELb1EE13estimate_dtauEP10PlasmaDataP9FieldData,_ZN13ParticleObjNTILi1ELi2ELi3ELb1EE18accumulate_currentEP10PlasmaDataP12CurrentTally8typevecNIS5_IdLi1EELi2EES5_IS6_Li3EERKS6_,_ZN13ParticleObjNTILi1ELi2ELi3ELb1EE18time_till_crossingER8typevecNIdLi1EES3_S3_d,_ZN13ParticleObjNTILi1ELi2ELi3ELb1EE4pushEP10PlasmaDataP9FieldDataP12CurrentTallyR8typevecNIiLi1EEi,_ZN13ParticleObjNTILi1ELi2ELi3ELb1EE7copy_inERK12ParticleListi,_ZN13ParticleObjNTILi1ELi2ELi3ELb1EE9get_accelEP9FieldDataRK8typevecNIS3_IdLi1EELi2EERKS3_IS4_Li3EE,_ZN14FieldDataCPU2D10intrpAccelEddddddiiiRdS0_S0_,_ZN14FieldDataCPU2D4getBEiiii,_ZN14FieldDataCPU2D4getEEiiii,_ZN14FieldDataCPU2D6getBTzILi0EEERdiii,_ZN14FieldDataCPU2D6getBTzILi1EEERdiii,_ZN14FieldDataCPU2D6getBTzILi2EEERdiii,_ZN14FieldDataCPU2D6getETzILi0EEERdiii,_ZN14FieldDataCPU2D6getETzILi1EEERdiii,_ZN14FieldDataCPU2D6getETzILi2EEERdiii,_ZN15ParticleListCPU10get_fvalueEii,_ZN15ParticleListCPU10get_ivalueEii,_ZN15ParticleListCPU14push_interfaceILi10EEExP10PlasmaDataP9FieldDataP9HOMoments,_ZN15ParticleListCPU14push_interfaceILi11EEExP10PlasmaDataP9FieldDataP9HOMoments,_ZN15ParticleListCPU14push_interfaceILi12EEExP10PlasmaDataP9FieldDataP9HOMoments,_ZN15ParticleListCPU14push_interfaceILi13EEExP10PlasmaDataP9FieldDataP9HOMoments,_ZN15ParticleListCPU14push_interfaceILi14EEExP10PlasmaDataP9FieldDataP9HOMoments,_ZN15ParticleListCPU14push_interfaceILi15EEExP10PlasmaDataP9FieldDataP9HOMoments,_ZN15ParticleListCPU4pushEP10PlasmaDataP9FieldDataP9HOMoments,_ZN15ParticleListCPU5pushTILi1ELi2ELi3ELb1EEExP10PlasmaDataP9FieldDataP9HOMoments,_ZNK12ParticleList9get_floatEi,_ZNK15ParticleListCPU10get_fvalueEii,_ZNK15ParticleListCPU10get_ivalueEii"		

#include $(TRILINOS_DIR)/include/Makefile.export.NOX

#Includes
INCLUDE_PATHS:= -I/usr/local/cuda/include \
				-I/usr/local/cuda-5.0/include \
				-I./Util \
				-I/home/opt/cuda/cuda-5.0/include \
				-I/usr/include/mpich2-x86_64 \
				-I./Field_Solvers/ES1D \
				-I/projects/opt/mpi/openmpi/1.6.3_centos_gcc-4.4.6/include \
				-I/users/payne/cuda/include 
#				$(NOX_INCLUDE_DIRS)
#				-I/opt/cudatoolkit-4.1/include 
				
#Library Path
LIBRARY_PATHS=/usr/local/cuda-5.0/lib64:/usr/lib64/mpich2/lib:$(LD_LIBRARY_PATH):/usr/lib64/openmpi/lib:$(NOX_LIBRARY_DIRS)

CXX_FLAGS= -O3 -std=gnu++0x -fopenmp  -march=native  $(INCLUDE_PATHS) #-L${LIBRARY_PATHS}
#CXX_FLAGS+=-bf-include=$(BYFL_INCLUDE)


NVCC_FLAGS= -O3 -arch=compute_20 --use_fast_math -code=sm_20  -Xcompiler $(INCLUDE_PATHS) -L${LIBRARY_PATHS} -Xcompiler -pg
#CXX_FLAGS= -Wall -O0 -fopenmp -g -std=gnu++0x $(INCLUDE_PATHS) -L${LIBRARY_PATHS}

#CXX_FLAGS += -ftree-slp-vectorize -fkeep-static-consts -ftree-ch -floop-interchange -floop-parallelize-all -funroll-loops 

#Josh-Workstation Library Config
#LIBS = -lrt  -lmpichcxx -L/usr/local/cuda-5.0/lib64 -lcudart -lcudadevrt $(NOX_LIBRARY_DIRS) $(NOX_LIBRARIES)#-lmpi_cxx #-lmpichcxx  
#Darwin Library Config
#LIBS = -lrt -lmpi_cxx -L/projects/opt/cuda/cuda-5.0/lib64 -lcudart #-L/home/cuda/cuda4.2/lib64 -lcudart 
#Moonlight Library config
LIBS = -lrt -lmpi_cxx -L/opt/cudatoolkit-4.1/lib64  -lcudart #-L/users/payne/cuda/lib64 -lcudart #

VPATH = ./:./Problem_Initializers:./Field_Solvers:../gnuplot_c/src:./Field_Solvers/ES1D:./Util:/usr/local/cuda-5.0/lib64:/usr/local/cuda/include

A := $(shell echo abc)
$(info $(MACHINE))

ifndef MACHINE
	LIBS = -lrt  -lmpichcxx
	CUDA_LIB_DIR=/usr/local/cuda-5.0/lib64
else
ifeq "$(MACHINE)" "moonlight"
	LIBS = -lrt -lmpi_cxx
	CUDA_LIB_DIR=/opt/cudatoolkit-4.1/lib64
else
ifeq "$(MACHINE)" "darwin"
	LIBS = -lrt -lmpi_cxx
	CUDA_LIB_DIR=/projects/opt/cuda/cuda-5.0/lib64
endif
endif
endif


# command line input for no hand vectorization
ifdef NOHANDVEC
	CXX_FLAGS+=-DNO_HAND_VEC
endif

#command line input for using cuda
ifneq ($(origin USECUDA),undefined)
	LIBS += -L$(CUDA_LIB_DIR) -lcudart
	OBJ += GPUObjects.o $(GPUHOSTOBJ)
else
	CXX_FLAGS+= -DNO_CUDA
endif



#Object Director
OBJDIR:=./obj
ASMDIR:=./asm
BINDIR:=./bin
TESTDIR:=./tests

EXECS:= ExB_test \
	TwoStream_test \
	Landau_test \
	EnergyCons_test \
	IonAcoustic_test \
	Island_test

OBJ:=ParticleObjN.o \
	ParticleObjN1d1v.o \
	ParticleObjN1d2v.o \
	ParticleObjN1d3v.o \
	ParticleObjN2d2v.o \
	ParticleObjN2d3v.o \
	ParticleObjN3d3v.o \
	FieldDataCPU.o \
	gnuplot_i.o \
	ParticleListCPU.o \
	ParticleListCPUAoS.o \
	CurrentTally.o \
	ChargeTally.o \
	StressTally.o \
	AmpereSolve.o \
	HOMoments.o \
	PlasmaData.o \
	ExplicitPIC.o \
	ExB_Initializer.o \
	ConstantSolver.o \
	ConsistencyTermES.o \
	AmpereSolver.o \
	TwoStream_Initializer.o \
	Landau_Initializer.o \
	EnergyCons_Initializer.o \
	ImplicitPIC.o \
	IonAcoustic_Initializer.o \
	LOSolverSI.o \
	LOMomentsES.o \
	NLResidual.o \
	PlasmaUtility.o \
	Coupled1D_ES.o \
	Island_Initializer.o \
	FieldDataCPU2D.o \
	ParticleListCPUSorted.o
	
	
GPUOBJ:=ParticleObjNGPU.cubin \
		CurrentTally.cubin \
		ChargeTally.cubin \
		StressTally.cubin \
		CurrentTallyGPU.cubin \
		HOMoments.cubin \

		
GPUHOSTOBJ:= ParticleListGPUSimple.cu.cpp.ii \
			ParticleListGPU.cu.cpp.ii \
			FieldDataGPU.cu.cpp.ii 
		

		
	
ASMS:=$(addprefix $(ASMDIR)/, $(OBJ:.o=.s))
	
OBJ:=$(addprefix $(OBJDIR)/, $(OBJ))

GPUOBJ:=$(addprefix $(OBJDIR)/, $(GPUOBJ))

GPUHOSTOBJ:=$(addprefix $(OBJDIR)/, $(GPUHOSTOBJ))

EXECS:=$(addprefix $(BINDIR)/, $(EXECS))
		
all: checkdirs libs

checkdirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(ASMDIR)
	mkdir -p $(BINDIR)
	

libs: $(OBJ) 

asm: $(ASMS)

gpustuff: $(GPUOBJ)



tests: $(EXECS)


$(BINDIR)/push_test : $(TESTDIR)/push_test.cpp $(OBJ)
	$(CXX) $(CXX_FLAGS) ${OBJ} $(TESTDIR)/push_test.cpp  $(CXX_FLAGS) -L${LIBRARY_PATHS} $(LIBS) -o push_test
	mv push_test $(BINDIR)/push_test
	
$(BINDIR)/% : $(TESTDIR)/%.cpp $(OBJ)
	$(CXX) $(CXX_FLAGS) ${OBJ} $(TESTDIR)/$*.cpp   $(CXX_FLAGS) -L${LIBRARY_PATHS} $(LIBS) -o $*
	mv $* $(BINDIR)/$*



# Rules for making object files
$(OBJDIR)/GPUObjects.o : $(GPUOBJ) $(GPUHOSTOBJ)
	$(NVCC) $(NVCC_FLAGS) --device-link $(GPUOBJ) $(GPUHOSTOBJ) -o $@

$(OBJDIR)/gnuplot_i.o : Util/gnuplot_i.c
	$(CXX) $(CXX_FLAGS) -c Util/gnuplot_i.c -o $(OBJDIR)/gnuplot_i.o $(CXX_FLAGS)

$(OBJDIR)/%.o : %.cpp
	$(CXX) -c $< -o $@ $(CXX_FLAGS)
	
$(OBJDIR)/%.o : %.cu
	$(CXX) $(CXX_FLAGS) -x c++ -c $< -o $@ $(CXX_FLAGS)

	
$(OBJDIR)/%.o : %.cxx
	$(CXX) $(CXX_FLAGS) -c $< -o $@ $(CXX_FLAGS)
	
$(OBJDIR)/%.cubin : %.cu
	$(NVCC) $(NVCC_FLAGS) --cubin -D GPU_CODE -DCUDA_CODE --device-c $< -o $@

ifeq ($(origin USECUDA),undefined)
$(OBJDIR)/%.cu.cpp.ii : %.cu
	$(NVCC) $(NVCC_FLAGS) --cuda -DCUDA_CODE --device-c $< -o $@
else
$(OBJDIR)/%.cu.cpp.ii : %.cu
	$(CXX) $(CXX_FLAGS) $< -o $@
endif

$(OBJDIR)/%.cubin : %.cpp
	$(NVCC)  --cubin --device-c $< -o $@ $(NVCC_FLAGS)
		
$(ASMDIR)/%.s : %.cpp
	$(CXX) $(CXX_FLAGS) -S $< -o $@ $(CXX_FLAGS)

# Special Rule for making ParticleObj files
$(OBJDIR)/%.o : ./ParticleObjs/%.cpp ParticleObjNT.inl
	$(CXX) $(CXX_FLAGS) -c $< -o $@ $(CXX_FLAGS)
	
$(ASMDIR)/%.s : ./ParticleObjs/%.cpp ParticleObjNT.inl
	$(CXX) $(CXX_FLAGS) -S $< -o $@ $(CXX_FLAGS)
	
$(OBJDIR)/%.cubin : ./ParticleObjs/%.cu ParticleObjNT.inl
	$(NVCC) $(NVCC_FLAGS) --cubin --device-c $< -o $@
	
cleantemps: *.cu.cpp
	rm -f -v *.cu.cpp

clean:
	rm -f -v *.o
	rm -f -v $(OBJDIR)/*.o
	rm -f -v $(OBJDIR)/*.gpuo
	rm -f -v $(OBJDIR)/*.cubin
	rm -f -v $(OBJDIR)/*.cu.cpp.ii
	rm -f -v $(OBJDIR)/*.bc

cleans:
	rm -f -v $(ASMDIR)/*.s
