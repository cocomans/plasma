PlasmaApp: A multi-architecture implicit particle-in-cell proxy
=========================================


### Building 

<dt><code>$>gmake</code></dt>
<dd>build all of the required libraries.

<dt><code>$>gmake tests</code></dt>
<dd>build all of the test routines

Be sure to use the parallel build option "-j N" where N is the number of
threads to use. 

Note: Double precision is toggled in the file PlasmaData.h via the preprocessor define
DOUBLE_PRECISION. to use single precision simply comment out this line of code.

Note 2: MPI libraries may be different on your machine. You may have to edit the makefile
to use the correct one.

Note 3: To build with cuda append USECUDA=1 to the make command


#############################
########## Running ##########

There are several test problems currently implemented.
1. Two Stream Instability
2. Ion Acoustic Shock


running the Two Stream Instability problem 
PlasmaApp3D $> mpirun -N $NUM_NODES -n $NUM_TASKS ./bin/TwoStream_test -np $NUM_PTCLS -nx 32 -Lx 1 -dt 0.5 -s 100

PlasmaApp3D $> mpirun -N $NUM_NODES -n $NUM_TASKS ./bin/IonAcoustic_test -np $NUM_PTCLS -nx 128 -Lx 144 -dt 0.5 -s 1000


##############################
### Command Line Arguments ###

-nx #, -ny #, -nz #
	Number of cells in the x, y, and z dimensions
	
-Lx #, -Ly #, -Lz #
	System Length in debye lengths
	
-x0 #, -y0 #, -z0 #
	System origin

--vec_length #
	Length of particle list object for cpu particle push

-dt #
	time step size in electron plasma frequency

-np #
	Number of particles per mpi task
	
-ns #
	Number of particle spiecies

--min-subcycles #
	Minimum number of subcycles to use during HO particle push

--num-cores #
	number of cpu cores to use for shared memory particle push.

-g
	Turns on plotting

--lo-all
	run the lo order solver on all nodes

--runid #
	Specify output file number. 