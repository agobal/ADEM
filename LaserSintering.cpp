#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SLS.h"
#include <cstring>

using namespace std;

TempProfile LaserSintering(struct PowderBed PB)
{
	TempProfile TP;
	// First step is to initialize the temperature and internal energy matrices
	float solid_heat_capacity = 40;
	float particle_heat_capacity;
	float preheat_temperature = 500;
	for (int i = 0; i < PB.cell_count; ++i)
	{
		for (int j = 0; j < PB.particle_count; ++j)
		{
			TP.T[i][j] = preheat_temperature;
			particle_heat_capacity = solid_heat_capacity*(4.0/3.0)*4.0*atan(1)*pow(PB.r_particles[j], 3);
			TP.E[i][j] = TP.T[i][j]*particle_heat_capacity;
		}
	}

	// Assign the simulation timestep (seconds)
	float delta_t = 0.00001;
	LaserPath LP;
	LP = GcodeReader(delta_t);
	
	// Finding laser path and total number of required time steps and laser location at each time by reading the gcode

	return TP;
}

