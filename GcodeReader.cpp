#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SLS.h"
#include <cstring>
#include <fstream>
#include <string>
#include </usr/include/python2.7/Python.h>
#include <sstream>

using namespace std;

LaserPath GcodeReader(float delta_t)
{
	// Define laserpath coordinates
	LaserPath LP;

	// initialize required parameters


	// First read the gcode provided by user
	// ifstream gcode("gcode.txt");
	// string str; 
	// int q;
	// q = 0;
	// while(getline(gcode, str))
	// {
	// 	if (q == 0)
	// 	{
	// 		LP.laser_speed = atof(str.c_str());;
	// 	}
	// 	q = q + 1;
	// 	// cout << str << endl;
	// }
	// gcode.close();

	float distance = sqrt(pow(0.0002 - 0.000, 2) + pow(0.0005 - 0.0005, 2));
	// LP.time_steps = int((distance/(LP.laser_speed))/(delta_t));
	LP.laser_speed = 0.001;
	LP.time_steps = 1000;

	for (int i = 0; i < LP.time_steps; ++i)
	{
		if (i == 0)
		{
			LP.x_laser[0] = 0;
			LP.y_laser[0] = 0.0001;
		}
		else
		{
			LP.x_laser[i] = LP.x_laser[i - 1] + (0.0002 - 0.000)/(float(LP.time_steps));
			LP.y_laser[i] = LP.y_laser[i - 1] + (0.0005 - 0.0005)/(float(LP.time_steps));
		}
	}
	return LP;

}

