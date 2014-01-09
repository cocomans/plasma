

CC=gcc#bf-gcc
ifdef BYFL_BUILD
CXX=mpic++ -cxx=bf-g++ -bf-verbose -bf-thread-safe -bf-vectors -bf-every-bb -bf-merge-bb=1000000000#-bf-call-stack -bf-by-func 
else
CXX=mpic++
endif

NVCC=nvcc

include ../Make.common

ifdef TRILINOS

include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.Teuchos
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.Epetra
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.EpetraExt
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.AztecOO
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.NOX
endif
#Includes
INCLUDE_PATHS:= -I/usr/local/cuda/include \
				-I/usr/local/cuda-5.0/include \
				-I./Util \
				-I/projects/opt/cuda/5.5/include \
				-I./Field_Solvers/ES1D \
				-I/users/payne/cuda/include \
				-I/usr/include/mpich2-x86_64 \
				-I/opt/cudatoolkit-5.0/include 
#				$(NOX_INCLUDE_DIRS)
#				-I/opt/cudatoolkit-4.1/include 
	
				
#Library Path
#LIBRARY_PATHS=/usr/local/cuda-5.0/lib64:/usr/lib64/mpich2/lib:$(LD_LIBRARY_PATH):/usr/lib64/openmpi/lib:/usr/lib64
LIBRARY_PATHS=/usr/local/cuda-5.0/lib64:$(LD_LIBRARY_PATH):/usr/lib64
CXX_FLAGS_ALL = -O3  -std=c++0x  -fopenmp -pthread -march=native -DNO_MIC
CXX_FLAGS= $(CXX_FLAGS_ALL) $(INCLUDE_PATHS)  #-L${LIBRARY_PATHS}
#CXX_FLAGS+=-Dshitty_timer# -bf-include=$(BYFL_INCLUDE) 
#CXX_FLAGS= -std=gnu++0x -openmp -g -pthread $(INCLUDE_PATHS) -DNO_MIC  #-L${LIBRARY_PATHS}


NVCC_FLAGS= -m64 -O3 -rdc=true -arch=compute_20 -code=sm_20 --ptxas-options=-v -Xcompiler $(INCLUDE_PATHS) -L${LIBRARY_PATHS} -Xcompiler -fopenmp -Xcompiler -pthread
#NVCC_FLAGS= -O3 -arch=compute_30 -code=sm_30 --ptxas-options=-v -Xcompiler $(INCLUDE_PATHS) -L${LIBRARY_PATHS} -Xcompiler -fopenmp
#CXX_FLAGS= -Wall -O0 -fopenmp -g -std=gnu++0x $(INCLUDE_PATHS) -L${LIBRARY_PATHS}

#CXX_FLAGS += -ftree-slp-vectorize -fkeep-static-consts -ftree-ch -floop-interchange -floop-parallelize-all -funroll-loops 

#Josh-Workstation Library Config
#LIBS = -lrt  -lmpichcxx -L/usr/local/cuda-5.0/lib64 -lcudart -lcudadevrt $(NOX_LIBRARY_DIRS) $(NOX_LIBRARIES)#-lmpi_cxx #-lmpichcxx  
#Darwin Library Config
#LIBS = -lrt -lmpi_cxx -L/projects/opt/cuda/cuda-5.0/lib64 -lcudart #-L/home/cuda/cuda4.2/lib64 -lcudart 
#Moonlight Library config
LIBS = -lrt -lmpi_cxx -L/opt/cudatoolkit-4.1/lib64  -lcudart #-L/users/payne/cuda/lib64 -lcudart #

VPATH+= ./:./Problem_Initializers:./Field_Solvers:../gnuplot_c/src:./Field_Solvers/ES1D:./Util:/usr/local/cuda-5.0/lib64:/usr/local/cuda/include:./:./Field_Solvers/HoLoInterface:./Normalizations:./Util/webui:./Util/IPICAdaptor:./Field_Solvers/HoLoInterface2D


A := $(shell echo abc)
$(info $(MACHINE))

ifdef NO_COMP_VECTORIZE
CXX_FLAGS+=-O2  -finline-functions -funswitch-loops -fpredictive-commoning -fgcse-after-reload -fipa-cp-clone -fvect-cost-model 
else
CXX_FLAGS+=-O3
endif

ifdef SINGLE_PRECISION
CXX_FLAGS+=-DNO_DOUBLE_PREC 
NVCC_FLAGS+=-DNO_DOUBLE_PREC 
endif


ifndef MACHINE
	LIBS = -lrt  -lmpichcxx 
	CUDA_LIB_DIR=/usr/local/cuda-5.0/lib64
else
ifeq "$(MACHINE)" "moonlight"
	LIBS = -lrt -lmpi_cxx 
	CUDA_LIB_DIR=/opt/cudatoolkit-5.0/lib64
else
ifeq "$(MACHINE)" "darwin"
	LIBS = -lrt -lmpi_cxx -L/projects/opt/libraries/numactl-2.0.7/lib64 
	CUDA_LIB_DIR=/projects/opt/cuda/5.5/lib64
endif
endif
endif


# command line input for no hand vectorization
ifdef NOHANDVEC
	CXX_FLAGS+=-DNO_HAND_VEC
endif

#command line input for using cuda
ifneq ($(origin USECUDA),undefined)
	LIBS += -L$(CUDA_LIB_DIR) -lcudart -lcudadevrt
else
	CXX_FLAGS+= -DNO_CUDA
endif


ifdef TRILINOS
CXX_FLAGS+=$(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) -I. -I$(TRILINOS_INSTALL_DIR)/include
LIBS  += -L$(TRILINOS_INSTALL_DIR)/lib
LIBS += $(Teuchos_LIBRARIES)
LIBS += $(Epetra_LIBRARIES)
LIBS += $(EpetraExt_LIBRARIES)
LIBS += $(Epetra_TPL_LIBRARIES)
LIBS += $(AztecOO_LIBRARIES)
LIBS += $(NOX_LIBRARIES)
NVCC_FLAGS+=$(Trilinos_INCLUDE_DIRS) -I. -I$(TRILINOS_INSTALL_DIR)/include

# Define include directories
INCFLAGS = -I. -I$(TRILINOS_INSTALL_DIR)/include
else
CXX_FLAGS +=-I/usr/include/mpich2-x86_64
endif

#ParaView_DIR:= /usr/projects/views/paraview/ParaView-v4.0.1-build


##Object Director
#OBJDIR:=$(PWD)/obj
#ASMDIR:=$(PWD)/asm
#BINDIR:=$(PWD)/bin
#TESTDIR:=$(PWD)/src/tests
#LIBDIR:=$(PWD)/lib
#PKGDIR:=$(PWD)/pkgs
#
#LIBADAPTOR_BD:=$(PKGDIR)/Adaptor-build
#LIBADAPTOR_SRCD:=$(PWD)/src/Util/IPICAdaptor
#
#LLIBS:=  libIPICAdaptor.so
#
#
#EXECS:= TwoStream_test \
#		IonAcoustic_test \
#		ExB_test \
#		Island_test \
#		PushVal_test \
#		Weibel_test
#
#
#OBJ:=NodeHOMoments.o \
#	NodeParticleList.o \
#	NodeFieldData.o \
#	FieldDataCPU.o \
#	ParticleObjN1d1v.o \
#	ParticleObjN1d2v.o \
#	ParticleObjN1d3v.o \
#	ParticleObjN2d2v.o \
#	ParticleObjN2d3v.o \
#	ParticleObjN3d3v.o \
#	gnuplot_i.o \
#	ParticleListCPU.o \
#	ChargeTallyCPU.o \
#	StressTallyCPU.o \
#	HOMomentsCPU.o \
#	ParallelInfo.o \
#	PlasmaData.o \
#	AmpereSolver.o \
#	ConstantSolver.o \
#	TwoStream_Initializer.o \
#	IonAcoustic_Initializer.o \
#	Island_Initializer.o \
#	ExB_Initializer.o \
#	PushVal_Initializer.o \
#	ImplicitPIC.o \
#	PlasmaUtility.o \
#	CurrentTallyCPU.o \
#	OutputCtrl.o \
#	mkpath.o \
#	Weibel_Initializer.o \
#	NormElectronEM.o \
#	rand.o \
#	sse-server.o \
#	LogStream.o \
#	sseDataStream.o \
#	
#	
#ifdef TRILINOS
#OBJ+=LowOrderProblem.o \
#		SimParams.o \
#		DiscretePDE_1D1Ves.o \
#		DiscretePDE_1D3Vem.o \
#		PhysBasedPrec.o \
#		DiffMath.o \
#		JFNK.o \
#		StagMesh.o \
#		MapManager.o \
#		MLSolver.o \
#		HoLoInterface.o
#endif
#	
#	
#GPUOBJ:=ParticleObjNGPU.cubin \
#		ChargeTallyGPU.cubin \
#		StressTallyGPU.cubin \
#		CurrentTallyGPU.cubin 
#
#		
#GPUHOSTOBJ:= FieldDataGPUSimple.cu.cpp.ii \
#			ParticleListGPUSorted.cu.cpp.ii \
#			HOMomentsGPUs.cu.cpp.ii \
#			GPU_Utils.cu.cpp.ii
#			
#			
#			
#ifneq ($(origin USECUDA),undefined)
#OBJ+= GPUObjects.o
#OBJ+= $(GPUHOSTOBJ)
#endif

ifdef USE_MIC
CC=icc
CXX_FLAGS+=-offload-attribute-target=mic -DUSEMIC
endif

		
	
#ASMS:=$(addprefix $(ASMDIR)/, $(OBJ:.o=.s))
#	
#OBJ:=$(addprefix $(OBJDIR)/, $(OBJ))
#
#LLIBS:=$(addprefix $(LIBDIR)/, $(LLIBS))
#
#
#GPUOBJ:=$(addprefix $(OBJDIR)/, $(GPUOBJ))
#
#GPUHOSTOBJ:=$(addprefix $(OBJDIR)/, $(GPUHOSTOBJ))
#
#EXECS:=$(addprefix $(BINDIR)/, $(EXECS))


		
all: obj libs tests

checkdirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(ASMDIR)
	mkdir -p $(BINDIR)
	mkdir -p $(PKGDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(PKGDIR)/trilinos
	mkdir -p $(PKGDIR)/mesa
	mkdir -p $(PKGDIR)/paraview
	
	

	

obj: $(OBJDIR) $(OBJ)

asm: $(ASMDIR) $(ASMS)

gpustuff: $(GPUOBJ)

libs: $(LIBDIR) $(PKGDIR) $(LLIBS)

tests: obj libs $(BINDIR) $(EXECS)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(ASMDIR):
	mkdir -p $(ASMDIR)

$(BINDIR): 
	mkdir -p $(BINDIR)

$(PKGDIR):
	mkdir -p $(PKGDIR)


$(LIBDIR):
	mkdir -p $(LIBDIR)




$(BINDIR)/push_test : $(TESTDIR)/push_test.cpp $(OBJ) $(LLIBS)
	$(CXX) $(CXX_FLAGS) ${OBJ}  $(TESTDIR)/push_test.cpp  $(CXX_FLAGS) -L${LIBRARY_PATHS} $(LIBS) -L$(TOP_DIR)/lib $(LinkLIBS) -o push_test
	mv push_test $(BINDIR)/push_test
	
$(BINDIR)/% : $(TESTDIR)/%.cpp $(OBJ) $(LLIBS)
	$(CXX) $(CXX_FLAGS)  ${OBJ}  $(TESTDIR)/$*.cpp   $(CXX_FLAGS) -L${LIBRARY_PATHS} $(LIBS) -L$(TOP_DIR)/lib $(LinkLIBS) -o $*
	mv $* $(BINDIR)/$*



# Rules for making object files
$(OBJDIR)/GPUObjects.o : $(GPUOBJ) $(GPUHOSTOBJ)
	$(NVCC) $(NVCC_FLAGS) --device-link $(GPUOBJ) $(GPUHOSTOBJ) -o $@

$(OBJDIR)/%.o : %.c
	$(CXX) $(CXX_FLAGS) -c $< -o $@ $(CXX_FLAGS)
	
#$(OBJDIR)/mkpath.o : Util/mkpath.cpp
#	$(CXX) $(CXX_FLAGS) -c Util/mkpath.cpp -o $(OBJDIR)/mkpath.o $(CXX_FLAGS)

$(OBJDIR)/%.o : %.cpp
	$(CXX) -c $< -o $@ $(CXX_FLAGS)
	
$(OBJDIR)/%.o : %.cu
	$(CXX) $(CXX_FLAGS) -x c++ -c $< -o $@ $(CXX_FLAGS)

	
$(OBJDIR)/%.o : %.cxx
	$(CXX) $(CXX_FLAGS) -c $< -o $@ $(CXX_FLAGS)
	
$(OBJDIR)/%.cubin : %.cu
	$(NVCC) $(NVCC_FLAGS) --cubin -DGPU_CODE -DCUDA_CODE --device-c $< -o $@
	
$(OBJDIR)/%.cu.o : %.cu
	$(NVCC) $(NVCC_FLAGS) -DGPU_CODE -DCUDA_CODE --device-c $< -o $@
	
$(OBJDIR)/%.fatbin : %.cu
	$(NVCC) $(NVCC_FLAGS) --fatbin -DGPU_CODE -DCUDA_CODE --device-c $< -o $@

ifneq ($(origin USECUDA),undefined)
$(OBJDIR)/%.cu.cpp.ii : %.cu
	$(NVCC) $(NVCC_FLAGS) -DGPU_CODE -DCUDA_CODE --device-c $< -o $@
else
$(OBJDIR)/%.cu.cpp.ii : %.cu
	$(CXX) $(CXX_FLAGS) $< -o $@
endif

$(OBJDIR)/%.cubin : %.cpp
	$(NVCC)  --cubin --device-c $< -o $@ $(NVCC_FLAGS)
		
$(ASMDIR)/%.s : %.cpp
	$(CXX) $(CXX_FLAGS) -S $< -o $@ $(CXX_FLAGS)

# Special Rule for making ParticleObj files
$(OBJDIR)/%.o : ./ParticleObjs/%.cpp ./ParticleObjNT.inl
	$(CXX) $(CXX_FLAGS) -c $< -o $@ $(CXX_FLAGS)
	
$(ASMDIR)/%.s : ./ParticleObjs/%.cpp ./ParticleObjNT.inl
	$(CXX) $(CXX_FLAGS) -S $< -o $@ $(CXX_FLAGS)
	
$(OBJDIR)/%.cu.cpp.ii : ./ParticleObjs/%.cu ./ParticleObjNT.inl
	$(NVCC) $(NVCC_FLAGS) --cuda  $< -o $@
	
$(OBJDIR)/%.cu.o : ./ParticleObjs/%.cu ./ParticleObjNT.inl
	$(NVCC) $(NVCC_FLAGS) --device-c  $< -o $@
	
$(LIBADAPTOR_BD)/CMakeCache.txt: $(LIBADAPTOR_SRCD)/IPICAdaptor.cxx $(LIBADAPTOR_SRCD)/IPICAdaptor.h $(PKGDIR)
	cd $(PKGDIR); \
	mkdir -p $(LIBADAPTOR_BD); \
	cd $(LIBADAPTOR_BD); \
	cmake -D ParaView_DIR:PATH=$(ParaView_DIR) -D CMAKE_CXX_FLAGS="$(CXX_FLAGS_ALL)" $(LIBADAPTOR_SRCD) 
		
$(LIBADAPTOR_BD)/libIPICAdaptor.so: $(LIBADAPTOR_BD)/CMakeCache.txt $(LIBADAPTOR_SRCD)/IPICAdaptor.cxx $(LIBADAPTOR_SRCD)/IPICAdaptor.h $(PKGDIR)
	cd $(LIBADAPTOR_BD); \
	gmake
	
$(LIBDIR)/libIPICAdaptor.so: $(LIBADAPTOR_BD)/libIPICAdaptor.so $(LIBADAPTOR_SRCD)/IPICAdaptor.cxx $(LIBADAPTOR_SRCD)/IPICAdaptor.h $(PKGDIR)
	ln --force $(LIBADAPTOR_BD)/libIPICAdaptor.so $(LIBDIR)


		
	


cleantemps: *.cu.cpp
	rm -f -v *.cu.cpp
	


clean:
	rm -f -v *.o
	rm -f -v $(OBJDIR)/*.o
	rm -f -v $(OBJDIR)/*.gpuo
	rm -f -v $(OBJDIR)/*.cubin
	rm -f -v $(OBJDIR)/*.cu.cpp.ii
	rm -f -v $(LIBDIR)/*.so

deepclean: clean
	rm -rfv $(PKGDIR)/Adaptor-build
	rm -rfv $(BINDIR)
	rm -rfv $(LIBDIR)
	rm -rfv $(ASMDIR)
	rm -rfv $(OBJDIR)

cleans:
	rm -f -v $(ASMDIR)/*.s
