#include <fstream>
#include <iostream>
#include <stdio.h>
#include "SLS.h"
#include <ctime>

using namespace std;

/* 
This program tries to perform a 3D full scale simulation of the selective laser sintering process.
A multiscale adaptive discrete element method is used for performing thermal simulations and the laser 
beam is modeled as a Gaussian heat input on top of the powder bed.
*/

int main()
{	
	/* Starting the chronometer */
	int start_s = clock();

	/* Powder bed dimensions and layer thickness (all in milimeters) */
	float bed_x, bed_y, bed_z, layer_thickness;

	bed_x = 1.0;
	bed_y = 1.0;
	bed_z = 2.0;
	layer_thickness = 0.05;

	/* Powder packing properties */
	ParticleChar PC1;	// Loading the powder packing data structre from SLS.h

	PC1.avgrd = 7.5/1000000.0;		// Average diameter of particles inside the packing
	PC1.stddev = 2.5/1000000.0;		// Standard deviation of the particles inside packing
	PC1.packfrac = 0.63;	// Packing fraction of the bed
	PC1.nmin = 700;		// Minimum number of particles inside the bed

	/* Making up the powder bed */
	PowderBed Bed;		// Loading the data structure for storing particles inside the powder bed
	Bed = PackingGenerator(PC1);	// Generating the desired packing 

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

	return 0;
}