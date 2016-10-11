#include <fstream>
#include <iostream>
#include <stdio.h>
#include "SLS.h"
#include <ctime>

using namespace std;

int main()
{
/*	int start_s = clock();
	//Defining bed and particle dimensions (all dimensions in SI)
	float bed_x, bed_y, bed_z;

	//Bed geometry
	bed_x = 0.2;
	bed_y = 0.2;
	bed_z = 0.2;

	//Defining powder characteristics
	ParticleChar PC1;

	PC1.avgrd = 15.0/1000000.0;
	PC1.stddev = 5.0/1000000.0;
	PC1.packfrac = 0.575;
	PC1.nmin = 700;

	//Making up the powder bed
	PowderBed Bed;
	Bed = PackingGenerator(PC1);

	//Making up the smaller packing within the large one
	SPowderBed SBed;
	SBed = SmallPackingGenerator(Bed);

	cout<<"Number of particles:"<<SBed.cnt<<endl;

	//Generating the shape function matrices and element stiffness and capacity matrices
	ElementCoefficients EC;
	EC = ShFcnGen(Bed, SBed);

	//Finding elements on sides



	//Writing powder bed data into .txt files
	OutputWriter(Bed);

	int stop_s = clock();
	cout<<"Stiffness matrices created"<<endl;
	cout << "time: " << (stop_s - start_s)/double(CLOCKS_PER_SEC) << endl;
*/
	return 0;
}