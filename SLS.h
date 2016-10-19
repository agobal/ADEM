struct PowderBed{
	float x_particles[8][130];
	float y_particles[8][130];
	float z_particles[8][130];
	float r_particles[130];
	int neighbors[8][130][15];
	int particle_count;
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

	extern float solid_conductivity;
	extern float air_conductivity;
	extern float Youngs_modulus;
	extern float Poissons_ratio;
	extern float E_prime;
	extern float outside_pressure;
	extern float average_particle_radius;
	extern float F;
	extern float M;
	extern float contact_radius;
	extern float contact_diameter;
	extern float center_distance;
	extern float contact_area;
	extern float solid_contact_area;
	extern float air_contact_area;
	extern float solid_heat_capacity;
	extern float liquid_heat_capacity;
	extern float solid_density;
	extern float liquid_density;
	extern float melting_specific_heat;
	extern float evaporation_specific_heat;
	extern float melting_temperature;
	extern float evaporation_temperature;
	extern float mass;

PowderBed PackingGenerator(struct ParticleChar PC, struct BedGeometry BG1, struct PowderBed PB);

/*SPowderBed SmallPackingGenerator(struct PowderBed Struct);

ElementCoefficients ShFcnGen(struct PowderBed PB, struct SPowderBed SPB);

ElementEdge EdgeFinder(struct PowderBed PB);*/

// int* CornerLocator(struct PowderBed PB, struct SPowderBed SPB);

float* RadiusFinder(struct ParticleChar PC, float grid_volume);

// float** PsiGenerator(struct PowderBed PB, struct ElementEdge);

// float** Inverse(float** PsiC);

// float** ElementStiffness(float* Ge, float** Ke);

// void OutputWriter(PowderBed);
