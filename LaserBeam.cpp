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
/* Function for calculating the laser beam power for the location of a particle at a certain point in time */
float LaserBeam(float x_p, float y_p, float z_p, float r_p, float v_l, float x_l, float y_l)
{
	float W = 10;
	float r_laser = 0.00005;
	float I = floor(z_p/0.00003)*((2*W)/(4.0*atan(1)*r_laser*r_laser))*pow(2.7183, (-(2*(pow((x_p - x_l), 2) + pow((y_p - y_l), 2)))/(pow(r_laser, 2))));

	return I;
}