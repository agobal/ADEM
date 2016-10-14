#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SLS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;

PowderBed PackingGenerator(struct ParticleChar PC, struct BedGeometry BG1, struct PowderBed PB)
{
	PowderBed Bed;

	// Add the particle size distribution to the new struct
	Bed.particle_count = PB.particle_count;
	for (int i = 0; i < PB.particle_count)
	{
		Bed.d_particles[i] = PB.d_particles[i];
	}

	// Assign random locations to particles
	int cell_x_num[BG1.num_grid], cell_y_num[BG1.num_grid], cell_z_num[BG1.num_grid];	// Number of the cell grid in x, y and z directions
	for (int i = 0; i < BG1.num_grid; ++i)
	{
		cell_x_num[i] = remainder(i - 1, BG1.num_grid_x) + 1;
		cell_y_num[i] = floor( remainder( (i - 1), (BG1.num_grid_x*BG1.num_grid_y) )/BG1.num_grid_x ) + 1;
		cell_z_num[i] = floor( (i - 1)/(BG1.num_grid_x*BG1.num_grid_y) ) + 1;
		for (int j = 0; j < PB.particle_count; ++j)
		{
			Bed.x_particles[i][j] = ((double) rand() / (RAND_MAX))*(BG1.grid_x - Bed.d_particles[i]) + Bed.d_particles[i] + (cell_x_num - 1)*BG1.grid_x;
			Bed.y_particles[i][j] = ((double) rand() / (RAND_MAX))*(BG1.grid_y - Bed.d_particles[i]) + Bed.d_particles[i] + (cell_y_num - 1)*BG1.grid_y;
			Bed.z_particles[i][j] = ((double) rand() / (RAND_MAX))*(BG1.grid_z - Bed.d_particles[i]) + Bed.d_particles[i] + (cell_z_num - 1)*BG1.grid_z;
		}
	}

	float q, x_particle_middle, y_particle_middle, z_particle_middle;

	// Relocate particles to reduce overlaps
	for (int c = 0; c < BG1.num_grid; ++c)
	{
		// Filling out the neighboring particle numbers for relocation
		float neighbor_particles[7*Bed.particle_count];
		for (int c1 = 0; c1 < Bed.particle_count; ++c1)
		{
			neighbor_particles[c1] = c1;
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for (int k = 0; k < 600; ++k)
		{
			for (int i = 0; i < Bed.particle_count; ++i)
			{
				q = 0.0;
				x_particle_middle = 0.0;
				y_particle_middle = 0.0;
				z_particle_middle = 0.0;
				
				for (int j = 0; j < Bed.particle_count; ++j)
				{
					if (j != i)
					{
						if ((Bed.d_particles[i]/2.0 + Bed.d_particles[j]/2.0) > sqrt(pow(Bed.x_particles[i] - Bed.x_particles[j], 2) + pow(Bed.y_particles[i] - Bed.y_particles[j], 2) + pow(Bed.z_particles[i] - Bed.z_particles[j], 2)))
						{
							q = q + 1;
							x_particle_middle = x_particle_middle + Bed.x_particles[j] + (Bed.x_particles[i] - Bed.x_particles[j])*(Bed.d_particles[i]/2.0 + Bed.d_particles[j]/2.0)/sqrt(pow((Bed.x_particles[i] - Bed.x_particles[j]),2) + pow((Bed.y_particles[i] - Bed.y_particles[j]),2) + pow((Bed.z_particles[i] - Bed.z_particles[j]),2));
							y_particle_middle = y_particle_middle + Bed.y_particles[j] + (Bed.y_particles[i] - Bed.y_particles[j])*(Bed.d_particles[i]/2.0 + Bed.d_particles[j]/2.0)/sqrt(pow((Bed.x_particles[i] - Bed.x_particles[j]),2) + pow((Bed.y_particles[i] - Bed.y_particles[j]),2) + pow((Bed.z_particles[i] - Bed.z_particles[j]),2));
							z_particle_middle = z_particle_middle + Bed.z_particles[j] + (Bed.z_particles[i] - Bed.z_particles[j])*(Bed.d_particles[i]/2.0 + Bed.d_particles[j]/2.0)/sqrt(pow((Bed.x_particles[i] - Bed.x_particles[j]),2) + pow((Bed.y_particles[i] - Bed.y_particles[j]),2) + pow((Bed.z_particles[i] - Bed.z_particles[j]),2));
						}
					}
				}

				if (q >= 1)
				{
					Bed.x_particles[i] = x_particle_middle/q;
					Bed.y_particles[i] = y_particle_middle/q;
					Bed.z_particles[i] = z_particle_middle/q;
				}

				//Don't accept if the new location is outside boundary
				if (Bed.xparticles[i] < Bed.rparticles[i])
					Bed.xparticles[i] = Bed.rparticles[i];
				if (Bed.xparticles[i] > (Bed.xel - Bed.rparticles[i]))
					Bed.xparticles[i] = Bed.xel - Bed.rparticles[i];
				if (Bed.yparticles[i] < Bed.rparticles[i])
					Bed.yparticles[i] = Bed.rparticles[i];
				if (Bed.yparticles[i] > (Bed.yel - Bed.rparticles[i]))
					Bed.yparticles[i] = Bed.yel - Bed.rparticles[i];
				if (Bed.zparticles[i] < Bed.rparticles[i])
					Bed.zparticles[i] = Bed.rparticles[i];
				if (Bed.zparticles[i] > (Bed.zel - Bed.rparticles[i]))
					Bed.zparticles[i] = Bed.zel - Bed.rparticles[i];
			}
		}
	}

	//Find the array of adjacent particles
/*	for (int i = 0; i < mystruct.nmin; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			Bed.neighbors[i][j] = -1;
		}
	}

	for (int i = 0; i < mystruct.nmin; ++i)
	{
		q = 0;
		for (int j = 0; j < mystruct.nmin; ++j)
		{
			if ((i != j) && ( abs( sqrt(pow((Bed.xparticles[i] - Bed.xparticles[j]),2) + pow((Bed.yparticles[i] - Bed.yparticles[j]),2) + pow((Bed.zparticles[i] - Bed.zparticles[j]),2)) - (Bed.rparticles[i] + Bed.rparticles[j]) ) < 0.000005 ))
			{
				Bed.neighbors[i][(int)q] = j;
				q = q + 1;
			}
		}
	}*/

	//typedef boost::tuple<std::string, int> animal;
	return Bed;
}