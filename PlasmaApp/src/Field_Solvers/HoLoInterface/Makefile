# Point to the location of Trilinos
TRILINOS_INSTALL_DIR=$(TRILINOS_DIR)

# Provide Makefile.export definitions
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.Teuchos
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.Epetra
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.EpetraExt
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.AztecOO
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.NOX

# Define libraries
LDFLAGS  = -L$(TRILINOS_INSTALL_DIR)/lib
LDFLAGS += $(Teuchos_LIBRARIES)
LDFLAGS += $(Epetra_LIBRARIES)
LDFLAGS += $(EpetraExt_LIBRARIES)
LDFLAGS += $(Epetra_TPL_LIBRARIES)
LDFLAGS += $(AztecOO_LIBRARIES)
LDFLAGS += $(NOX_LIBRARIES)

# Define include directories
INCFLAGS = -I. -I$(TRILINOS_INSTALL_DIR)/include

#default: main.cpp LowOrderProblem.o SimParams.o DiscretePDE.o PhysBasedPrec.o JFNK.o DiffMath.o StagMesh.o MapManager.o MLSolver.o HoLoInterface.o
#	mpic++ -g -O0 -Wall $(INCFLAGS) main.cpp LowOrderProblem.o SimParams.o DiscretePDE.o PhysBasedPrec.o DiffMath.o JFNK.o StagMesh.o MapManager.o MLSolver.o HoLoInterface.o $(LDFLAGS) -o run.x

LowOrderProblem.o: LowOrderProblem.cpp LowOrderProblem.h SimParams.h DiscretePDE.h PhysBasedPrec.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c LowOrderProblem.cpp -o LowOrderProblem.o

SimParams.o: SimParams.cpp SimParams.h MLSolver.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c SimParams.cpp -o SimParams.o

StagMesh.o: StagMesh.cpp StagMesh.h SimParams.h MapManager.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c StagMesh.cpp -o StagMesh.o

MapManager.o: MapManager.cpp MapManager.h SimParams.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c MapManager.cpp -o MapManager.o

DiscretePDE.o: DiscretePDE.cpp DiscretePDE.h SimParams.h DiffMath.h StagMesh.h MapManager.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c DiscretePDE.cpp -o DiscretePDE.o

PhysBasedPrec.o: PhysBasedPrec.cpp PhysBasedPrec.h SimParams.h DiffMath.h MapManager.h MLSolver.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c PhysBasedPrec.cpp -o PhysBasedPrec.o

DiffMath.o: DiffMath.cpp DiffMath.h MapManager.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c DiffMath.cpp -o DiffMath.o

JFNK.o: JFNK.cpp JFNK.h SimParams.h LowOrderProblem.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c JFNK.cpp -o JFNK.o

MLSolver.o: MLSolver.cpp MLSolver.h SimParams.h DiffMath.h MapManager.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c MLSolver.cpp -o MLSolver.o

HoLoInterface.o: HoLoInterface.cpp HoLoInterface.h SimParams.h JFNK.h LowOrderProblem.h
	mpic++ -g -O0 -Wall $(INCFLAGS) -c HoLoInterface.cpp -o HoLoInterface.o

clean:
	rm -rf *.o *.x *~
