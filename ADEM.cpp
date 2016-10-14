#include <fstream>
#include <iostream>
#include <stdio.h>
#include "SLS.h"
#include <ctime>
#include <math.h>
#include <typeinfo>
#include <vector>

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

	/* Grid size and number of grids (all in milimeters) */
	float grid_x, grid_y, grid_z, grid_volume;
	int num_grid_x, num_grid_y, num_grid_z;

	grid_x = 0.1;	// Size of grid-x
	grid_y = 0.1;	// Size of grid-y
	grid_z = 0.05;	// Size of grid-z

	grid_volume = grid_x*grid_y*grid_z;		// Volume of each grid

	num_grid_x = int(bed_x/grid_x);		// Number of grids in x-direction
	num_grid_y = int(bed_y/grid_y);		// Number of grids in y-direction
	num_grid_z = int(bed_z/grid_z);		// Number of grids in z-direction

	static int num_grid = num_grid_x*num_grid_y*num_grid_z;	// Total number of grids in the bed

	float bed_volume;
	bed_volume = bed_x*bed_y*bed_z;

	/* Powder packing properties */
	ParticleChar PC1;	// Loading the powder packing data structre from SLS.h

	PC1.avgrd = 15.0/1000.0;		// Average diameter of particles inside the packing
	PC1.stddev = 5.0/1000.0;		// Standard deviation of the particles inside packing
	PC1.packfrac = 0.63;	// Packing fraction of the bed
	PC1.nmin = grid_volume/((4.0/3.0)*4.0*atan(1.0)*pow(((PC1.avgrd)/2.0), 3.0));		// Minimum number of particles inside each grid

	/* Finding out the number of powders inside each grid */
	float* p_d;
	p_d = DiameterFinder(PC1, grid_volume);

	static int num_particle_grid = int(p_d[0]);

	// Creating the new array for particle diameters
	float particle_diameter[num_particle_grid];
	for (int i = 0; i < num_particle_grid; ++i)
	{
		particle_diameter[i] = p_d[i + 1];
	}

	// Initialize the arrays of particle locations and their neighbors
	struct PowderBed{
		PowderBed(): x_p(num_grid, num_particle_grid), y_p(num_grid, num_particle_grid), z_p(num_grid, num_particle_grid)  { }
		vector<float> x_p;
		vector<float> y_p;
		vector<float> z_p;
	};
	PowderBed PB;
	cout << PB.x_p[num_grid*num_particle_grid] << endl;
	/* Making up the powder bed */
	// Bed = PackingGenerator(PC1);	// Generating the desired packing 

/*	//Making up the smaller packing within the large one
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
	cout<<"Stiffness matrices created"<<endl;*/
	// cout << "time: " << (stop_s - start_s)/double(CLOCKS_PER_SEC) << endl;

	return 0;
}