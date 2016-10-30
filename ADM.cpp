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
int ADM(int cell, float x_p[], float y_p[], float T[], float x_l, float y_l, int t)
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
		if ((x_p[i] != 0) && (y_p[i] != 0))
		{
			sum_x += x_p[i];
			sum_y += y_p[i];
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
		if (T[i] > max_temp)
		{
			max_temp = T[i];
			q = i;
		}
	}

	float dist = pow(pow(avg_x - x_l, 2) + pow(avg_y - y_l, 2), 0.5);
	if ((max_temp > 1000) || (dist < 0.00015))
		final_decision = 1;

	return final_decision;
}
