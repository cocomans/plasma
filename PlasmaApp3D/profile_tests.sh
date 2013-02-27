#!/bin/sh

echo TEST_ONE
#mpirun -np 1  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 1 --num_cores 1 -ns 1 --lo-all --min-subcycles 10

#mpirun -np 1  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 32 --num_cores 1 -ns 1 --lo-all --min-subcycles 10

#mpirun -np 1  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 32 --num_cores 48 -ns 1 --lo-all --min-subcycles 10

#mpirun -np 60  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 32 --num_cores 1 -ns 1 --lo-all --min-subcycles 10

#mpirun -np 60  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 32 --num_cores 48 -ns 1 --lo-all --min-subcycles 10

gmake clean
gmake tests -j 48 MACHINE=darwin NOHANDVEC=1

VEC_CONFIG='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32'
for i in $VEC_CONFIG
	do

		mpirun -np 1 --bind-to-none ./bin/TwoStream_test -np 500000 -nx 32 \
					-Lx 1 -s 50 -dt 1.0 --vec_length $i --num_cores 48 \
					-ns 1 --runid 21802 --lo-all
	done
	
mpirun -np 1 --bind-to-none ./bin/TwoStream_test -np 500000 -nx 32 -Lx 1 -s 50 -dt 1.0 --vec_length 32 --num_cores 48 -ns 1 --runid 21802 --lo-all


CORE_CONFIG='1 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48'
for i in $CORE_CONFIG
	do

		mpirun -np 1 --bind-to-none ./bin/TwoStream_test -np 500000 -nx 32 \
					-Lx 1 -s 50 -dt 1.0 --vec_length 32 --num_cores $i \
					-ns 1 --runid 21803
	done
	
	


NODE_CONFIG='1 2 3 4 5 6 7 8 9 10'


for i in $NODE_CONFIG
	do

		mpirun -np $i --bind-to-none ./bin/TwoStream_test -np 500000 -nx 32 \
					-Lx 1 -s 50 -dt 1.0 --vec_length 32 --num_cores 48 \
					-ns 1 --runid 21804
	done

GRID_SIZE='32 64 128 256 512'
	
for i in $GRID_SIZE
	do

		mpirun -np 10 --bind-to-none ./bin/TwoStream_test -np 500000 -nx $i \
					-Lx 1 -s 50 -dt 1.0 --vec_length 32 --num_cores 48 \
					-ns 1 --runid 21805
	done
	
DT_SWEEP='0.1 0.5 1.0 1.5 2.0'
for i in $GRID_SIZE
	do

		mpirun -np 10 --bind-to-none ./bin/TwoStream_test -np 500000 -nx 32 \
					-Lx 1 -s 50 -dt $i --vec_length 32 --num_cores 48 \
					-ns 1 --runid 21806
	done


	
