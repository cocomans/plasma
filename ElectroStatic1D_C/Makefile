# $Id: Makefile,v 1.12 2012/05/22 21:34:07 pkf Exp $

# Override compiler chosen by system, for Mac or whatever
#CXX     = mpicxx

#include /home/pkf/software/rru/trunk/halo_finder/include.mk

PROGS    = \
	PlasmaDriver

CXX = g++

CXXFLAGS = -g

# Probably need to fix dependencies on some of the header files.
HDR= \
	PlasmaData.h \
	PlasmaUtility.h \
	Particles.h \
	ParticleInitialize.h \
	InitialLandau.h \
	InitialTSI.h \
	InitialIASW.h \
	ParticlePush.h \
	PushSimple.h \
	PushRayTrace.h \
	MomentSolution.h \
	ParticleSolution.h \
	NLResidual.h \
	LOSolver.h \
	LOSolverJFNK.h \
	LOSolverSI.h \
	DhatConsistency.h \
	PlasmaOut.h \
	OutGnuplot.h \
	AmpereVlasov.h

SRC= \
	PlasmaData.cxx \
	PlasmaUtility.cxx \
	Particles.cxx \
	ParticleInitialize.cxx \
	InitialLandau.cxx \
	InitialTSI.cxx \
	InitialIASW.cxx \
	PushSimple.cxx \
	PushRayTrace.cxx \
	MomentSolution.cxx \
	ParticleSolution.cxx \
	NLResidual.cxx \
	LOSolverJFNK.cxx \
	LOSolverSI.cxx \
	DhatConsistency.cxx \
	OutGnuplot.cxx \
	AmpereVlasov.cxx

OBJ= \
	PlasmaData.o \
	PlasmaUtility.o \
	Particles.o \
	ParticleInitialize.o \
	InitialLandau.o \
	InitialTSI.o \
	InitialIASW.o \
	PushSimple.o \
	PushRayTrace.o \
	MomentSolution.o \
	ParticleSolution.o \
	NLResidual.o \
	LOSolverJFNK.o \
	LOSolverSI.o \
	DhatConsistency.o \
	OutGnuplot.o \
	AmpereVlasov.o

PlasmaDriver: PlasmaDriver.cxx ${OBJ}
	${CXX} ${CXXFLAGS} -o $@ $^ ${CXXFLAGS} ${LDFLAGS}

PlasmaData.o: PlasmaData.cxx PlasmaData.h
	${CXX} ${CXXFLAGS} -c PlasmaData.cxx ${CXXFLAGS}

Particles.o: Particles.cxx Particles.h
	${CXX} ${CXXFLAGS} -c Particles.cxx ${CXXFLAGS}

ParticleInitialize.o: ParticleInitialize.cxx ParticleInitialize.h
	${CXX} ${CXXFLAGS} -c ParticleInitialize.cxx ${CXXFLAGS}

InitialLandau.o: InitialLandau.cxx InitialLandau.h ParticleInitialize.h
	${CXX} ${CXXFLAGS} -c InitialLandau.cxx ${CXXFLAGS}

InitialTSI.o: InitialTSI.cxx InitialTSI.h ParticleInitialize.h
	${CXX} ${CXXFLAGS} -c InitialTSI.cxx ${CXXFLAGS}

InitialIASW.o: InitialIASW.cxx InitialIASW.h ParticleInitialize.h
	${CXX} ${CXXFLAGS} -c InitialIASW.cxx ${CXXFLAGS}

PushSimple.o: PushSimple.cxx PushSimple.h ParticlePush.h
	${CXX} ${CXXFLAGS} -c PushSimple.cxx ${CXXFLAGS}

PushRayTrace.o: PushRayTrace.cxx PushRayTrace.h ParticlePush.h
	${CXX} ${CXXFLAGS} -c PushRayTrace.cxx ${CXXFLAGS}

PlasmaUtility.o: PlasmaUtility.cxx PlasmaUtility.h
	${CXX} ${CXXFLAGS} -c PlasmaUtility.cxx ${CXXFLAGS}

MomentSolution.o: MomentSolution.cxx MomentSolution.h
	${CXX} ${CXXFLAGS} -c MomentSolution.cxx ${CXXFLAGS}

ParticleSolution.o: ParticleSolution.cxx ParticleSolution.h
	${CXX} ${CXXFLAGS} -c ParticleSolution.cxx ${CXXFLAGS}

NLResidual.o: NLResidual.cxx NLResidual.h
	${CXX} ${CXXFLAGS} -c NLResidual.cxx ${CXXFLAGS}

LOSolverJFNK.o: LOSolverJFNK.cxx LOSolverJFNK.h LOSolver.h
	${CXX} ${CXXFLAGS} -c LOSolverJFNK.cxx ${CXXFLAGS}

LOSolverSI.o: LOSolverSI.cxx LOSolverSI.h LOSolver.h
	${CXX} ${CXXFLAGS} -c LOSolverSI.cxx ${CXXFLAGS}

DhatConsistency.o: DhatConsistency.cxx DhatConsistency.h
	${CXX} ${CXXFLAGS} -c DhatConsistency.cxx ${CXXFLAGS}

OutGnuplot.o: OutGnuplot.cxx OutGnuplot.h
	${CXX} ${CXXFLAGS} -c OutGnuplot.cxx ${CXXFLAGS}

AmpereVlasov.o: AmpereVlasov.cxx AmpereVlasov.h
	${CXX} ${CXXFLAGS} -c AmpereVlasov.cxx ${CXXFLAGS}

clean:
	rm *.o PlasmaDriver
