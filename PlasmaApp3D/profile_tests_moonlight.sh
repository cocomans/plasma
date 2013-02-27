#!/bin/sh

echo TEST_ONE
#mpirun -np 1  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 1 --num_cores 1 -ns 1 --lo-all --min-subcycles 10

#mpirun -np 1  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 32 --num_cores 1 -ns 1 --lo-all --min-subcycles 10

#mpirun -np 1  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 32 --num_cores 48 -ns 1 --lo-all --min-subcycles 10

#mpirun -np 60  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 32 --num_cores 1 -ns 1 --lo-all --min-subcycles 10

#mpirun -np 60  ./bin/TwoStream_test -np 1000000 -nx 128 -Lx 1 -s 25 -dt 0.5 --vec_length 32 --num_cores 48 -ns 1 --lo-all --min-subcycles 10


#GRID_CONFIG='32 64 128 256'
#for i in $GRID_CONFIG
#	do

#		mpirun -n 32 --bind-to-none  ./bin/Island_test -np 1024000 -nx $i -ny $i \
#					-s 10 -dt 10.0 --vec_length 32 --num_cores 32 \
#					-ns 1 --lo-all --runid 21826
#	done
	


#mpirun -n 40 ./bin/IonAcoustic_test -np 200000 -nx 64 -Lx 144 -s 10 -g -dt 3.0 -ns 2 --gpu-mult 1 --num_cores 16 --lo-all 





NODE_CONFIG='1 10'

#		mpirun -n 1  ./bin/TwoStream_test -np 1000000 -nx 32 \
#					-Lx 1 -s 50 -dt 1.0 --vec_length 32 --num_cores 1 \
#					-ns 1 --lo-all --runid 21

VEC_CONFIG='1 32'
for i in $VEC_CONFIG
	do

		mpirun -n $i --bind-to-none  ./bin/TwoStream_test -np 1000000 -nx 32 \
					-Lx 1 -s 25 -dt 2.0 --vec_length $i --num_cores 1 \
					-ns 1 --runid 21807
	done

for i in $NODE_CONFIG
	do

		mpirun -n $i --bind-to-none  ./bin/TwoStream_test -np 1000000 -nx 32 \
					-Lx 1 -s 25 -dt 2.0 --vec_length 32 --num_cores 16 \
					-ns 1 --lo-all --runid 21807
	done
	
#for i in $NODE_CONFIG
#	do

#		srun -N $i -n $i  ./bin/TwoStream_test -np 500000 -nx 32 \
#					-Lx 1 -s 20 -dt 0.5 --vec_length 32 --num_cores 32 \
#					-ns 1 --min-subcycles 10 --runid 2
#	done

#VEC_CONFIG='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32'
#for i in $VEC_CONFIG
#	do

#		srun -N 1 -n 1  ./bin/TwoStream_test -np 500000 -nx 128 \
#					-Lx 1 -s 50 -dt 0.5 --vec_length $i --num_cores 32 \
#					-ns 1 --min-subcycles 10 --runid 5
#	done
	
#CORE_CONFIG='1 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50'
#for i in $CORE_CONFIG
#	do

#		srun -N 1 -n 1  ./bin/TwoStream_test -np 1000000 -nx 32 \
#					-Lx 1 -s 50 -dt 0.5 --vec_length 32 --num_cores $i \
#					-ns 1 --min-subcycles 10 --runid 4
#	done

exit