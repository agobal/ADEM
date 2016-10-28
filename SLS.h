struct PowderBed{
	float x_particles[200][130];
	float y_particles[200][130];
	float z_particles[200][130];
	float r_particles[130];
	int neighbors[200][130][20];
	int particle_count;
	int cell_count;
};

struct TempProfile{
	float T[200][130];
	float T_temp[200][130];
	float E[200][130];
	float E_temp[200][130];
};

struct LaserPath{
	int time_steps;
	float laser_speed;
	float x_laser[100000];
	float y_laser[100000];
};

struct SPowderBed{
	float x_particles[800];
	float y_particles[800];
	float z_particles[800];
	float r_particles[800];
	int neighbors[800][10];
	int numbers[800];
	float x_el;
	float y_el;
	float z_el;
	int cnt;
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
	float bed_volume;
	float grid_x;
	float grid_y;
	float grid_z;
	float grid_volume;
	int num_grid_x;
	int num_grid_y;
	int num_grid_z;
	int num_grid;
};

struct ElementCoefficients{
	float K[8][8];
	float C[8][8];
};

struct ElementEdge{
	int x0y0[15];
	int x0yf[15];
	int x0z0[15];
	int x0zf[15];
	int xfy0[15];
	int xfyf[15];
	int xfz0[15];
	int xfzf[15];
	int y0z0[15];
	int y0zf[15];
	int yfz0[15];
	int yfzf[15];
};



// PowderBed PackingGenerator(struct ParticleChar PC, struct BedGeometry BG1, struct PowderBed PB);
PowderBed PackingGenerator(ParticleChar, BedGeometry, PowderBed);

TempProfile LaserSintering(struct PowderBed PB);

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