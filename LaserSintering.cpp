#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SLS.h"
#include <cstring>

using namespace std;

TempProfile LaserSintering(PowderBed PB, PowderBed PB_BIG)
{
	// Find the closest particle to each large element
	int closest_particle[PB_BIG.cell_count][PB_BIG.particle_count];
	float mid_dist;
	for (int cell = 0; cell < PB_BIG.cell_count; ++cell)
	{
		for (int i = 0; i < PB_BIG.particle_count; ++i)
		{
			mid_dist = 1;
			for (int cell2 = 0; cell2 < PB.cell_count; ++cell2)
			{
				for (int j = 0; j < PB.particle_count; ++j)
				{
					if (pow(pow((PB_BIG.x_particles[cell][i] - PB.x_particles[cell2][j]), 2) + pow((PB_BIG.y_particles[cell][i] - PB.y_particles[cell2][j]), 2) + pow((PB_BIG.z_particles[cell][i] - PB.z_particles[cell2][j]), 2), 0.5) < mid_dist)
					{
						mid_dist = pow(pow((PB_BIG.x_particles[cell][i] - PB.x_particles[cell2][j]), 2) + pow((PB_BIG.y_particles[cell][i] - PB.y_particles[cell2][j]), 2) + pow((PB_BIG.z_particles[cell][i] - PB.z_particles[cell2][j]), 2), 0.5);
						closest_particle[cell][i] = cell2*1000 + j;
					}
				}
			}
		}
	}

	
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
	float delta_t = 0.001;
	LaserPath LP;
	// Finding laser path and total number of required time steps and laser location at each time by reading the gcode
	LP = GcodeReader(delta_t);

	// Solving the actual laser sintering equations
	float K, Q, I, S;	// Heat transfer coefficient between two particles in contact
	float K_ab = 0.3;
	float rho = 7800;
	float C_s = 477;
	int cel, par; // Intermediate variables for finding the closest particles to the large elements
	int neigh, cell2;	// Middle parameter for setting the neighboring particle
	int Adaptive_decision_maker;	// Variable that finds out if we want to perform the small scale analysis on this cell or not
	for (int t = 0; t < LP.time_steps; ++t)
	{
		cout << t << endl;
		// Solving for the small (particle size) packing
		for (int cell = 0; cell < PB.cell_count; ++cell)
	    {
	    	// Finding out if we want to do the small scale solution on this cell or not
	    	Adaptive_decision_maker = ADM(cell, PB.x_particles[cell], PB.y_particles[cell], TP.T[cell], LP.x_laser[t], LP.y_laser[t], t);
	    	if (Adaptive_decision_maker == 1)
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
						// Laser power getting into the bed
						I = LaserBeam(PB.x_particles[cell][i], PB.y_particles[cell][i], PB.z_particles[cell][i], PB.r_particles[i], LP.laser_speed, LP.x_laser[t], LP.y_laser[t]);
						S = 4.0*atan(1)*PB.r_particles[i]*PB.r_particles[i];		// Particle surface absorbing the laser powder
						TP.E[cell][i] = TP.E[cell][i] + (Q + K_ab*S*I)*delta_t/(rho*(4.0/3.0)*4.0*atan(1)*pow(PB.r_particles[i], 3)); //particle energy increase by laser
						TP.T_temp[cell][i] = TP.E[cell][i]/C_s;	// Particle temperature change
					}
				}
			}
		}
		// Solving for the large (element size) packing
		for (int cell = 0; cell < PB_BIG.cell_count; ++cell)
	    {
		    for (int i = 0; i < PB_BIG.particle_count; ++i)
		    {   
		    	cel = (closest_particle[cell][i]/1000);
		    	par = (closest_particle[cell][i] % 1000);
				if (TP.T[cel][par] > 500)
				{
					TP.T_BIG[cell][i] = TP.T[cel][par];
					TP.E_BIG[cell][i] = TP.T_BIG[cell][i]*C_s*15;
				}
		    	Q = 0;
				for (int j = 0; j < 15; ++j)
				{
					if (PB_BIG.neighbors[cell][i][j] != 0)
					{
						if  (PB_BIG.neighbors[cell][i][j] > 1000)
						{
							neigh = PB_BIG.neighbors[cell][i][j];
							K = CondCoeff(PB_BIG.x_particles[cell][i], PB_BIG.y_particles[cell][i], PB_BIG.z_particles[cell][i], PB_BIG.r_particles[i], PB_BIG.x_particles[cell][neigh], PB_BIG.y_particles[cell][neigh], PB_BIG.z_particles[cell][neigh], PB_BIG.r_particles[neigh], TP.T[cell][i], TP.T[cell][neigh]);
							Q = Q + K*(TP.T_BIG[cell][i] - TP.T_BIG[cell][neigh]);
						}
						else
						{
							cell2 = (PB_BIG.neighbors[cell][i][j]/1000) - 1;
							if (cell2 != 0)
								neigh = (PB_BIG.neighbors[cell][i][j] % (cell2*1000));
							K = CondCoeff(PB_BIG.x_particles[cell][i], PB_BIG.y_particles[cell][i], PB_BIG.z_particles[cell][i], PB_BIG.r_particles[i], PB_BIG.x_particles[cell2][neigh], PB_BIG.y_particles[cell2][neigh], PB_BIG.z_particles[cell2][neigh], PB_BIG.r_particles[neigh], TP.T[cell][i], TP.T[cell2][neigh]);
							Q = Q + K*(TP.T_BIG[cell][i] - TP.T_BIG[cell2][neigh]);
						}
					}
					// Laser power getting into the bed
					I = LaserBeam(PB_BIG.x_particles[cell][i], PB_BIG.y_particles[cell][i], PB_BIG.z_particles[cell][i], PB_BIG.r_particles[i], LP.laser_speed, LP.x_laser[t], LP.y_laser[t]);
					S = 4.0*atan(1)*PB_BIG.r_particles[i]*PB_BIG.r_particles[i];		// Particle surface absorbing the laser powder
					TP.E_BIG[cell][i] = TP.E_BIG[cell][i] + (Q + K_ab*S*I)*delta_t/(rho*(4.0/3.0)*4.0*atan(1)*pow(PB_BIG.r_particles[i], 3)); //particle energy increase by laser
					TP.T_temp_BIG[cell][i] = TP.E_BIG[cell][i]/C_s;	// Particle temperature change
				}
			}
		}
		// Replacing the actual temps with temporary temps
		for (int cell = 0; cell < PB.cell_count; ++cell)
	    {
	    	for (int i = 0; i < PB.particle_count; ++i)
	    	{
	    		TP.T[cell][i] = TP.T_temp[cell][i];
	    	}
	    }
		for (int cell = 0; cell < PB_BIG.cell_count; ++cell)
	    {
	    	for (int i = 0; i < PB_BIG.particle_count; ++i)
	    	{
	    		TP.T_BIG[cell][i] = TP.T_temp_BIG[cell][i];
	    	}
	    }
	}

	return TP;
}

