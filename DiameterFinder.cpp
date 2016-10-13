#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <new>
#include "SLS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;

/* This function finds the number of particles as well as their size distribution inside one
grid cell of the powder bed. The number of particles is stored in the first item of the 
particle diameter array so the size of it is (number of particles + 1) */

float* DiameterFinder(struct ParticleChar PC1, float grid_volume)
{

	/* First create the size distribution using the boost library */
	boost::mt19937 rng;

	boost::normal_distribution<> nd(PC1.avgrd, PC1.stddev);

	boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);

	/* Assigning diameter values based on the normal distribution */
	float diam[PC1.nmin];	// First guess of number of particles (more than the actual number)
	float d;	// Temporary saver of particle diameter
	int num;	// Total number of particles inside the grid cell
	num = 0;

	float volume = 0;	// Total volume of the grid cell occupied by the powder particles

	/* Adding particle volumes to reach the volume of the powder bed */
	for (int i = 0; i <= PC1.nmin; ++i)
	{
		d = var_nor();
		volume = volume + (4.0/3.0)*4.0*atan(1.0)*pow((d/2.0), 3.0);
		// If the volume is cell than the maximum packing volume, add more
		if (volume < grid_volume*PC1.packfrac)
		{
			diam[i] = d;
			num = num + 1;
		}
		else
		{
			break;
		}
	}

	/* Making a smaller array of diameters */
	float final_diam[num + 1];
	final_diam[0] = num;	// Adding the number of particles to the array to read it later
	for (int i = 1; i <= num; ++i)
	{
		final_diam[i] = diam[i];
	}
	/* We need to convert the array into a pointer to be able to return it to the main program. the following line does that */
	float *final_diam_pointer = final_diam;
	return final_diam_pointer;
}