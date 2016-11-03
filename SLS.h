struct PowderBed{
	float x_particles[100][150];
	float y_particles[100][150];
	float z_particles[100][150];
	float r_particles[150];
	int neighbors[100][150][20];
	int particle_count;
	int cell_count;
};

struct TempProfile{
	float T[100][150];
	float T_BIG[100][150];
	float T_time[100][150];
	float T_BIG_time[100][150];
	float T_temp[100][150];
	float T_temp_BIG[100][150];
	float E[100][150];
	float E_BIG[100][150];
	float E_temp[100][150];
	float E_temp_BIG[100][150];
};

struct LaserPath{
	int time_steps;
	float laser_speed;
	float x_laser[1000];
	float y_laser[1000];
};

struct ParticleChar{
	float avgrd;
	float stddev;
	float packfrac;
	int nmin;
};

struct BedGeometry{
	float bed_x;
	float bed_y;
	float bed_z;
	float layer_thickness;
	float layer_thickness_big;
	float bed_volume;
	float grid_x;
	float grid_y;
	float grid_z;
	float grid_z_big;
	float grid_volume;
	float grid_volume_big;
	int num_grid_x;
	int num_grid_y;
	int num_grid_z;
	int num_grid_z_big;
	int num_grid;
	int num_grid_big;
};




// PowderBed PackingGenerator(struct ParticleChar PC, struct BedGeometry BG1, struct PowderBed PB);
PowderBed PackingGenerator(ParticleChar, BedGeometry, PowderBed);

TempProfile LaserSintering(PowderBed PB, PowderBed PB_BIG, int output_timestep);

LaserPath GcodeReader(float delta_t);

/*SPowderBed SmallPackingGenerator(struct PowderBed Struct);

ElementCoefficients ShFcnGen(struct PowderBed PB, struct SPowderBed SPB);

ElementEdge EdgeFinder(struct PowderBed PB);*/

// int* CornerLocator(struct PowderBed PB, struct SPowderBed SPB);

float* RadiusFinder(struct ParticleChar PC, float grid_volume);

// float** PsiGenerator(struct PowderBed PB, struct ElementEdge);

// float** Inverse(float** PsiC);

// float** ElementStiffness(float* Ge, float** Ke);

// void OutputWriter(PowderBed);

// Function for calculating the heat transfer coefficient between two particles based on their position and temperatures
float CondCoeff(float x1, float y1, float z1, float r1, float x2, float y2, float z2, float r2, float T1, float T2);

// Function for calculating the laser beam powder for a certain particle
float LaserBeam(float x_p, float y_p, float z_p, float r_p, float v_l, float x_l, float y_l);

// Function for finding out if we want to do a small scale analysis on the cell or not
int ADM(int cell, float x_p[], float y_p[], float T[], float x_l, float y_l, int t);