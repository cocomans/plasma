PlasmaApp: A multi-architecture implicit particle-in-cell proxy
=========================================
PlasmaApp is a flexible implicit charge and energy conserving implicit PIC framework.
This codes aims to demonstrate the potential of using a fluid plasma model to accelerate a kinetic model 
through a High-Low order system coupling. The multi-granularity of this problem gives it the ability to 
map well to emerging heterogeneous architectures with multiple levels of parallelism. Additionally this problem 
maps very well to very fine parallel architectures, such as GPUs, since the vast majority of the work is 
encapsulated in the particle system, a trivially parallel problem. This approach also has applicability to 
very large scale systems, potential exascale, due to the large amount of particle work per communication.

The initial C++ implementation targets hybrid GPU + Multi-Core systems, but will preserve the flexibility to be easily implemented on other architectures.  

This flexibility will be accomplished by separating the physics algorithms from the underlying architecture considerations through the use of C++ templates and class inheritance.



### Building 

<dt><code>$>gmake</code></dt>
<dd>build all of the required libraries.

<dt><code>$>gmake tests</code></dt>
<dd>build all of the test routines

Be sure to use the parallel build option <code>-j N</code> where N is the number of
threads to use. 

Note: Double precision is toggled in the file PlasmaData.h via the preprocessor define
<code>DOUBLE_PRECISION</code>. to use single precision simply comment out this line of code.

Note 2: MPI libraries may be different on your machine. You may have to edit the makefile
to use the correct one.

#### Make arguments 
<dt><code>USECUDA=1</dt></code>
<dd> Enables and builds CUDA parts of the code (Requires CUDA 5.0 or later)

<dt><code>NOHANDVEC=1</dt></code>
<dd> Disables hand vectorization, This should be used if your machine does not support AVX




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

<dt><code>--epsilon #</dt></code>
<dd>Specify tolerance of particle picard loop

<dt><code>--nSpatial #</dt></code>
<dd>Number of spatial dimensions to use

<dt><code>-nVel #</dt></code>
<dd>Number of velocity dimensions

<dt><code>--plist-cpu #</dt></code>
<dd>Which CPU particle list optimization to use 0=default, 1=sorted

<dt><code>--min-subcycles #</dt></code>
<dd>Minimum number of subcycles to use during HO particle push

<dt><code>--num-cores #</dt></code>
<dd>number of cpu cores to use for shared memory particle push.

<dt><code>--gpu-mult #</dt></code>
<dd>Multiplier for number of particles to run on gpu vs multi-core for load balancing.

<dt><code>-g</dt></code>
<dd>Turns on plotting

<dt><code>--lo-all</dt></code>
<dd>run the lo order solver on all nodes

<dt><code>--runid #</dt></code>
<dd>Specify benchmark output file number. 

## Questions?
Email Joshua Payne payne@lanl.gov
