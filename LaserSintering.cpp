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
	// Initializing temperatures for all particles inside the powder bed
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
	// Finding laser path and total number of required time steps and laser location at each time by reading the gcode
	LP = GcodeReader(delta_t);

	// Solving the actual laser sintering equations
	float K, Q;	// Heat transfer coefficient between two particles in contact
	int neigh, cell2;	// Middle parameter for setting the neighboring particle
	for (int t = 0; t < LP.time_steps; ++t)
	{
		for (int cell = 0; cell < PB.cell_count; ++cell)
	    {
	    	for (int i = 0; i < PB.particle_count; ++i)
	    	{   
	    		Q = 0;
				for (int j = 0; j < 15; ++j)
				{
					if (PB.neighbors[cell][i][j] != 0)
					{
						if  (PB.neighbors[cell][i][j] > 1000)
						{
							neigh = PB.neighbors[cell][i][j];
							K = CondCoeff(PB.x_particles[cell][i], PB.y_particles[cell][i], PB.z_particles[cell][i], PB.r_particles[i], PB.x_particles[cell][neigh], PB.y_particles[cell][neigh], PB.z_particles[cell][neigh], PB.r_particles[neigh], TP.T[cell][i], TP.T[cell][neigh]);
							Q = Q + K*(TP.T[cell][i] - TP.T[cell][neigh]);
						}
						else
						{
							cell2 = (PB.neighbors[cell][i][j]/1000) - 1;
							if (cell2 != 0)
								neigh = (PB.neighbors[cell][i][j] % (cell2*1000));
							K = CondCoeff(PB.x_particles[cell][i], PB.y_particles[cell][i], PB.z_particles[cell][i], PB.r_particles[i], PB.x_particles[cell2][neigh], PB.y_particles[cell2][neigh], PB.z_particles[cell2][neigh], PB.r_particles[neigh], TP.T[cell][i], TP.T[cell2][neigh]);
							Q = Q + K*(TP.T[cell][i] - TP.T[cell2][neigh]);
						}
					}
					// if ( sqrt( (x(cell, i) - V_c*t)^2 + (y(cell, i) - 0.0005)^2 ) < r_b )
					// 	I = Laser(x(cell, i), y(cell, i), z(cell, i), t, V_c);
					// else
					// 	I = 0;
					// S = pi*r(1, i)*r(1, i);
					// Ee(cell, i, q + 1) = Ee(cell, i, q) + (Q + K_ab*S*I)*dt/(rho_s*(4/3)*pi*r(1, i)^3);
					// T(cell, i, q + 1) = Ee(cell, i, q + 1)/C_s;
				}
			}
		}
	}

	return TP;
}

