#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SLS.h"
#include <cstring>

using namespace std;

PowderBed PackingGenerator(struct ParticleChar PC, struct BedGeometry BG1, struct PowderBed PB)
{
	PowderBed Bed;

	// Add the particle size distribution to the new struct
	Bed.particle_count = PB.particle_count;
	for (int i = 0; i < PB.particle_count; ++i)
	{
		Bed.r_particles[i] = PB.r_particles[i];
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
			Bed.x_particles[i][j] = ((double) rand() / (RAND_MAX))*(BG1.grid_x - Bed.r_particles[i]) + Bed.r_particles[i] + float((cell_x_num[i] - 1))*BG1.grid_x;
			Bed.y_particles[i][j] = ((double) rand() / (RAND_MAX))*(BG1.grid_y - Bed.r_particles[i]) + Bed.r_particles[i] + float((cell_y_num[i] - 1))*BG1.grid_y;
			Bed.z_particles[i][j] = ((double) rand() / (RAND_MAX))*(BG1.grid_z - Bed.r_particles[i]) + Bed.r_particles[i] + float((cell_z_num[i] - 1))*BG1.grid_z;
		}
	}

	float q, x_particle_middle, y_particle_middle, z_particle_middle;
	int counter, neighbor_part, neighbor_cell;
	// Relocate particles to reduce overlaps
	for (int c = 0; c < BG1.num_grid; ++c)
	{
		cout<<c<<endl;
		// Filling out the neighboring particle numbers for relocation
		float neighbor_particles[4*Bed.particle_count];
		counter = 0;
		for (int ct = 0; ct < 4*Bed.particle_count; ++ct)
			neighbor_particles[ct] = 0;
		for (int c1 = 0; c1 < Bed.particle_count; ++c1)
		{
			neighbor_particles[counter] = c1;
			counter = counter + 1;
		}
		if ((BG1.num_grid_x != 1) && (cell_x_num[c] != 1))
		{
			for (int c1 = 0; c1 < Bed.particle_count; ++c1)
			{
				neighbor_particles[counter] = c1 + c*1000;
				counter = counter + 1;
			}
		}
		if ((BG1.num_grid_y != 1) && (cell_y_num[c] != 1))
		{
			for (int c1 = 0; c1 < Bed.particle_count; ++c1)
			{
				neighbor_particles[counter] = c1 + c*1000;
				counter = counter + 1;
			}
		}
		if ((BG1.num_grid_z != 1) && (cell_z_num[c] != 1))
		{
			for (int c1 = 0; c1 < Bed.particle_count; ++c1)
			{
				neighbor_particles[counter] = c1 + c*1000;
				counter = counter + 1;
			}
		}
		// Performing the relocation of particles to reach a stable position
		for (int k = 0; k < 600; ++k)
		{
			for (int i = 0; i < Bed.particle_count; ++i)
			{
				q = 0.0;
				x_particle_middle = 0.0;
				y_particle_middle = 0.0;
				z_particle_middle = 0.0;

				for (int j = 0; j < 4*Bed.particle_count; ++j)
				{
					if (j < Bed.particle_count)
					{
						if (j != i)
						{
							if ((Bed.r_particles[i] + Bed.r_particles[j]) > sqrt(pow(Bed.x_particles[c][i] - Bed.x_particles[c][j], 2) + pow(Bed.y_particles[c][i] - Bed.y_particles[c][j], 2) + pow(Bed.z_particles[c][i] - Bed.z_particles[c][j], 2)))
							{
								q = q + 1;
								x_particle_middle = x_particle_middle + Bed.x_particles[c][j] + (Bed.x_particles[c][i] - Bed.x_particles[c][j])*(Bed.r_particles[i] + Bed.r_particles[j])/sqrt(pow((Bed.x_particles[c][i] - Bed.x_particles[c][j]),2) + pow((Bed.y_particles[c][i] - Bed.y_particles[c][j]),2) + pow((Bed.z_particles[c][i] - Bed.z_particles[c][j]),2));
								y_particle_middle = y_particle_middle + Bed.y_particles[c][j] + (Bed.y_particles[c][i] - Bed.y_particles[c][j])*(Bed.r_particles[i] + Bed.r_particles[j])/sqrt(pow((Bed.x_particles[c][i] - Bed.x_particles[c][j]),2) + pow((Bed.y_particles[c][i] - Bed.y_particles[c][j]),2) + pow((Bed.z_particles[c][i] - Bed.z_particles[c][j]),2));
								z_particle_middle = z_particle_middle + Bed.z_particles[c][j] + (Bed.z_particles[c][i] - Bed.z_particles[c][j])*(Bed.r_particles[i] + Bed.r_particles[j])/sqrt(pow((Bed.x_particles[c][i] - Bed.x_particles[c][j]),2) + pow((Bed.y_particles[c][i] - Bed.y_particles[c][j]),2) + pow((Bed.z_particles[c][i] - Bed.z_particles[c][j]),2));
							}
						}
					}
					else
					{
						neighbor_cell = (neighbor_particles[j]/1000);
						neighbor_part = remainder(neighbor_particles[j], 1000);
						if ((Bed.r_particles[i] + Bed.r_particles[neighbor_part]) > sqrt(pow(Bed.x_particles[c][i] - Bed.x_particles[neighbor_cell][neighbor_part], 2) + pow(Bed.y_particles[c][i] - Bed.y_particles[neighbor_cell][neighbor_part], 2) + pow(Bed.z_particles[c][i] - Bed.z_particles[neighbor_cell][neighbor_part], 2)))
						{
							q = q + 1;
							x_particle_middle = x_particle_middle + Bed.x_particles[neighbor_cell][neighbor_part] + (Bed.x_particles[c][i] - Bed.x_particles[neighbor_cell][neighbor_part])*(Bed.r_particles[i] + Bed.r_particles[neighbor_part])/sqrt(pow((Bed.x_particles[c][i] - Bed.x_particles[neighbor_cell][neighbor_part]),2) + pow((Bed.y_particles[c][i] - Bed.y_particles[neighbor_cell][neighbor_part]),2) + pow((Bed.z_particles[c][i] - Bed.z_particles[neighbor_cell][neighbor_part]),2));
							y_particle_middle = y_particle_middle + Bed.y_particles[neighbor_cell][neighbor_part] + (Bed.y_particles[c][i] - Bed.y_particles[neighbor_cell][neighbor_part])*(Bed.r_particles[i] + Bed.r_particles[neighbor_part])/sqrt(pow((Bed.x_particles[c][i] - Bed.x_particles[neighbor_cell][neighbor_part]),2) + pow((Bed.y_particles[c][i] - Bed.y_particles[neighbor_cell][neighbor_part]),2) + pow((Bed.z_particles[c][i] - Bed.z_particles[neighbor_cell][neighbor_part]),2));
							z_particle_middle = z_particle_middle + Bed.z_particles[neighbor_cell][neighbor_part] + (Bed.z_particles[c][i] - Bed.z_particles[neighbor_cell][neighbor_part])*(Bed.r_particles[i] + Bed.r_particles[neighbor_part])/sqrt(pow((Bed.x_particles[c][i] - Bed.x_particles[neighbor_cell][neighbor_part]),2) + pow((Bed.y_particles[c][i] - Bed.y_particles[neighbor_cell][neighbor_part]),2) + pow((Bed.z_particles[c][i] - Bed.z_particles[neighbor_cell][neighbor_part]),2));
						}
					}
				}

				if (q >= 1)
				{
					Bed.x_particles[c][i] = x_particle_middle/q;
					Bed.y_particles[c][i] = y_particle_middle/q;
					Bed.z_particles[c][i] = z_particle_middle/q;

					// Particles can't go over the overall boundaries of the bed
					if ((cell_x_num[c] == BG1.num_grid_x) && (Bed.x_particles[c][i] >= (BG1.grid_x*cell_x_num[c] - Bed.r_particles[i])))
						Bed.x_particles[c][i] = (BG1.grid_x - Bed.r_particles[i]);
					if ((cell_y_num[c] == BG1.num_grid_y) && (Bed.y_particles[c][i] >= (BG1.grid_y*cell_y_num[c] - Bed.r_particles[i])))
						Bed.y_particles[c][i] = (BG1.grid_y - Bed.r_particles[i]);
					// if (cell_z_num[c] == BG1.num_grid_z) && (Bed.z_particles[c][i] >= (BG1.grid_z*cell_z_num[c] - Bed.r_particles[i]))
					// 	Bed.z_particles[c][i] = (BG1.grid_z - Bed.r_particles[i]);
					// Particles can't go below a certain amount of the previous cell
					if (Bed.x_particles[c][i] <= (0.9*BG1.grid_x*(cell_x_num[c] - 1.0) + Bed.r_particles[i]))
						Bed.x_particles[c][i] = (0.9*BG1.grid_x*(cell_x_num[c] - 1.0) + Bed.r_particles[i]);
					if (Bed.y_particles[c][i] <= (0.9*BG1.grid_y*(cell_y_num[c] - 1.0) + Bed.r_particles[i]))
						Bed.y_particles[c][i] = (0.9*BG1.grid_y*(cell_y_num[c] - 1.0) + Bed.r_particles[i]);
					if (Bed.z_particles[c][i] <= (0.9*BG1.grid_z*(cell_z_num[c] - 1.0) + Bed.r_particles[i]))
						Bed.z_particles[c][i] = (0.9*BG1.grid_z*(cell_z_num[c] - 1.0) + Bed.r_particles[i]);
					// Particles can't go below the general size of the bed
					if ((cell_x_num[c] == 1) && (Bed.x_particles[c][i] <= (BG1.grid_x*(cell_x_num[c] - 1) + Bed.r_particles[i])))
						Bed.x_particles[c][i] = Bed.r_particles[i];
					if ((cell_y_num[c] == 1) && (Bed.y_particles[c][i] <= (BG1.grid_y*(cell_y_num[c] - 1) + Bed.r_particles[i])))
						Bed.y_particles[c][i] = Bed.r_particles[i];
					if ((cell_z_num[c] == 1) && (Bed.z_particles[c][i] <= (BG1.grid_z*(cell_z_num[c] - 1) + Bed.r_particles[i])))
						Bed.z_particles[c][i] = Bed.r_particles[i];
				}
			}
		}
	}
	return Bed;
}