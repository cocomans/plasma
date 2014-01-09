#ifndef RAND_GEN_H
#define RAND_GEN_H


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <tr1/random>

class RandGen
{
public:

	RandGen(int seed_in)
	{

		srand(seed_in);

		seeder.seed((rand()%2386157)+(rand()%98463));

		// combined engine
		xorgen = *(new std::tr1::xor_combine<std::tr1::ranlux64_base_01,7,std::tr1::mt19937,13>(gen1, gen2));

		// initialize the combined generator
		xorgen.seed(seeder);



		rand_unif = *(new std::tr1::uniform_int<int>(0,1000000-1));

	}

	double uniform()
	{

		return (rand_unif(xorgen)+1.0e-6)*1.0e-6;

	}
private:


	std::tr1::uniform_int<int> rand_unif;

	std::tr1::ranlux64_base_01 gen1;
	std::tr1::mt19937 gen2;
	std::tr1::linear_congruential<unsigned long, 33614, 0, 2147483647> seeder;
	std::tr1::ranlux64_base_01 generator;


	std::tr1::xor_combine<std::tr1::ranlux64_base_01,7,std::tr1::mt19937,13> xorgen;


};

#endif /* RAND_GEN_H */
