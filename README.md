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



### Running

There are several test problems currently implemented.
1. Two Stream Instability
2. Ion Acoustic Shock


running the Two Stream Instability problem 

<dt><code>$> mpirun -N $NUM_NODES -n $NUM_TASKS ./bin/TwoStream_test -np $NUM_PTCLS -nx 32 -Lx 1 -dt 0.5 -s 100</code></dt>

<dt><code>$> mpirun -N $NUM_NODES -n $NUM_TASKS ./bin/IonAcoustic_test -np $NUM_PTCLS -nx 128 -Lx 144 -dt 0.5 -s 1000</code></dt>



### Command Line Arguments

<dt><code>-nx #, -ny #, -nz #</dt></code>
<dd>Number of cells in the x, y, and z dimensions
	
<dt><code>-Lx #, -Ly #, -Lz #</dt></code>
<dd>System Length in debye lengths
	
<dt><code>-x0 #, -y0 #, -z0 #</dt></code>
<dd>System origin

<dt><code>--vec_length #</dt></code>
<dd>Length of particle list object for cpu particle push

<dt><code>-dt #</dt></code>
<dd>time step size in electron plasma frequency

<dt><code>-np #</dt></code>
<dd>Number of particles per mpi task
	
<dt><code>-ns #</dt></code>
<dd>Number of particle spiecies

<dt><code>--min-subcycles #</dt></code>
<dd>Minimum number of subcycles to use during HO particle push

<dt><code>--num-cores #</dt></code>
<dd>number of cpu cores to use for shared memory particle push.

<dt><code>-g</dt></code>
<dd>Turns on plotting

<dt><code>--lo-all</dt></code>
<dd>run the lo order solver on all nodes

<dt><code>--runid #</dt></code>
<dd>Specify output file number. 