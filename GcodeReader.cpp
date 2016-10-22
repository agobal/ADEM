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
	ifstream gcode("gcode.txt");
	string str; 
	string s;
	while (getline(gcode, str))
	{
		istringstream iss(str);
		while (iss >> s)
		{
			string mid = s[0];
			if (mid == "V")
			{
				// LP.laser_speed = atof(s.erase(0, 1));
				cout << s.erase(0, 1) << endl;
			}
		}
	}
	gcode.close();
	return LP;

}

