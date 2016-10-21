#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SLS.h"
#include <cstring>
#include <fstream>
#include <string>

using namespace std;

LaserPath GcodeReader(float delta_t)
{
	// Define laserpath coordinates
	LaserPath LP;

	// First read the gcode provided by user
	ifstream gcode("Gcode.txt");
	string line;

	while( getline(gcode, line) )
	{
		cout << line << endl;
	}
	return LP;
}

