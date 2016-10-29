#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SLS.h"
#include <cstring>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
/* Function for determining if we want to perform small scale calculations on the cell or not */
int ADM(int cell, PowderBed PB, TempProfile TP, LaserPath LP, int t)
{
	int final_decision;
	final_decision = 0; // Don't perform small scale on this grid

	// Finding the average location of the grid cell
	float avg_x, avg_y; // Average location of the grid
	float sum_x, sum_y;
	float q;
	q = 0;
	sum_x = 0;
	sum_y = 0;
	for (int i = 0; i < 150; ++i)
	{
		if ((PB.x_particles[cell][i] != 0) && (PB.y_particles[cell][i] != 0))
		{
			sum_x += PB.x_particles[cell][i];
			sum_y += PB.y_particles[cell][i];
			q = q + 1;
		}
	}
	avg_x = ((float)sum_x)/q;
	avg_y = ((float)sum_y)/q;

	// Finding maximum temperature of the grid cell
	float max_temp;
	max_temp = 0;

	for (int i = 0; i < 150; ++i)
	{
		if (TP.T[cell][i] > max_temp)
		{
			max_temp = TP.T[cell][i];
			q = i;
		}
	}

	float dist = pow(pow(avg_x - LP.x_laser[t], 2) + pow(avg_y - LP.y_laser[t], 2), 0.5);
	if ((max_temp > 1000) || (dist < 0.00015))
		final_decision = 1;

	return final_decision;
}
