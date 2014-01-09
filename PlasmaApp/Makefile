

CC=gcc#bf-gcc
CXX=mpic++# -cxx=bf-g++ -bf-verbose -bf-thread-safe -bf-all-ops #-bf-by-func 
NVCC=nvcc

TOP_DIR=$(PWD)
include Make.local
include Make.common







TRILINOS_INCLUDES:= Makefile.export.Epetra \
					Makefile.export.EpetraExt \
					Makefile.export.AztecOO \
					Makefile.export.NOX
#					Makefile.export.Teuchos

TRILINOS_INCLUDES:=$(addprefix $(TRILINOS_INSTALL_DIR)/include/, $(TRILINOS_INCLUDES))


	

# Rules to build source code
all: 
	echo "Building the PlasmaApp Framework"
	cd $(SRCDIR); \
	gmake -j 48 tests TOP_DIR=$(TOP_DIR) $(MFLAGS)
	
packages: $(TRILINOS_INCLUDES) $(MESALIB_LIB) $(ParaViewInstallDir)/bin/paraview
	
	
# Rules to make trilinos include stuff
$(TRILINOS_INSTALL_DIR)/include/Makefile.export.Teuchos: $(SRCDIR)/Scripts/do-configure.sh $(SRCDIR)/Scripts/Makefile.trilinos
	echo "Trilinos is not detected, Configuring and Building Trilinos" 
	mkdir -p $(PKGDIR)/trilinos
	cp $(SRCDIR)/Scripts/do-configure.sh $(TRILINOS_PKG_DIR)/do-configure.sh
	cp $(SRCDIR)/Scripts/Makefile.trilinos $(TRILINOS_PKG_DIR)/Makefile
	cd $(TRILINOS_PKG_DIR); \
	gmake all TOP_DIR=$(TOP_DIR) $(MFLAGS)
	
$(TRILINOS_INCLUDES):$(TRILINOS_INSTALL_DIR)/include/Makefile.export.Teuchos

$(MESALIB_LIB): $(SRCDIR)/Scripts/Makefile.mesa
	mkdir -p $(MESALIB_DIR)
	cp $(SRCDIR)/Scripts/Makefile.mesa $(MESALIB_DIR)/Makefile
	echo "Building OSMesa"
	cd $(MESALIB_DIR); \
	gmake all TOP_DIR=$(TOP_DIR)

$(ParaViewInstallDir)/bin/paraview:$(MESA_LIB) $(SRCDIR)/Scripts/Makefile.paraview
	mkdir -p $(PKGDIR)/paraview
	cp $(SRCDIR)/Scripts/Makefile.paraview $(PKGDIR)/paraview/Makefile
	cp $(SRCDIR)/Scripts/pv-configure.sh $(PKGDIR)/paraview/pv-configure.sh
	cd $(PKGDIR)/paraview; \
	gmake all TOP_DIR=$(TOP_DIR)
	
	

checkdirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(ASMDIR)
	mkdir -p $(BINDIR)
	mkdir -p $(PKGDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(PKGDIR)/mesa
	mkdir -p $(PKGDIR)/paraview

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


$(TRILINOS_PKG_DIR):
	mkdir -p $(TRILINOS_PKG_DIR)	

cleantemps: *.cu.cpp
	rm -f -v *.cu.cpp
	
cudaclean:
	rm -f -v $(OBJDIR)/*.gpuo
	rm -f -v $(OBJDIR)/*.cubin
	rm -f -v $(OBJDIR)/*.cu.cpp.ii
	rm -f -v $(OBJDIR)/GPUObjects.o
	rm -f -v $(OBJDIR)/*.fatbin
	rm -fv $(OBJDIR)/*.cu.o

clean:
	rm -rfv $(BINDIR)
	rm -rfv $(LIBDIR)
	rm -rfv $(ASMDIR)
	rm -rfv $(OBJDIR)

deepclean: clean
	rm -rfv $(PKGDIR)
	rm -rfv $(BINDIR)
	rm -rfv $(LIBDIR)
	rm -rfv $(ASMDIR)
	rm -rfv $(OBJDIR)
	
pkgclean: clean
	cd $(PKGDIR)/paraview; \
	gmake realclean
	cd $(PKGDIR)/mesa; \
	gmake realclean
	cd $(PKGDIR)/trilinos; \
	gmake realclean

cleans:
	rm -f -v $(ASMDIR)/*.s
	
	





