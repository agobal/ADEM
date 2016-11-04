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
/* Calculation of the heat conduction coefficient during the sintering process is based on the following works :
	A theoretical model for effective thermal conductivity (ETC) of particulate beds under compression, Galit Weidenfeld et al, 2004
	Simulation of powder sintering using a discrete element model, Jerzy rojek et al, 2013
*/
ParticleDiffusion CondCoeff(float x1, float y1, float z1, float r1, float x2, float y2, float z2, float r2, int sinter_flag, float T1, float T2, float delta_t)
{
	ParticleDiffusion PD;
	// This function calculates the heat conduction coefficient between two powder particles based on their location and temperatures
	// Thermal properties
	float k1 = 40; // Thermal conductivity of particle 1
	float k2 = 40; // Thermal conductivity of particle 2
	float k_air = 0.02; // Thermal conductivity of air
	float alpha1 = 1;
	float alpha2 = 1;
	float k = (2*k1*k2)/(k1 + k2); // Average thermal conductivity of the two surfaces;
	float F = 1.0; // The contact force between the two particles
	float E = 207000000000.0; // Young's modulus of the material
	float nu = 0.33; // Poisson's ratio of the material
	float E_prime = E/(1 - pow(nu, 2));
	float r_eq = (2*r1*r2)/(r1 + r2);
	float m = 0.1; // Average slope of surface roughness
	float sigma = 0.00002; // RMS roughness height
	float r_c, P, A_a, A_r, A_v, cc_dist;
	int s_flag;

	// Mechanical properties
	if (sinter_flag == 0)
	{
		r_c = pow(((3.0/16.0)*r_eq*F/E_prime) , 0.33); // Contact area radius between two circles
		P = F/(4.0*atan(1)*pow(r_c, 2)); // Compression pressure
		A_a = 4*atan(1)*pow(r_c, 2); // Total contact area
		A_r = A_a*(P*1.41)/(E_prime*m); // Contact point area
		A_v = A_a - A_r; // Void contact area

		cc_dist = r1*cos(asin(r_c/r1)) + r2*cos(asin(r_c/r2));	// Distance of center to center particles
		// Move particle 2 to distance
		PD.x_new = x2 + (r1 + r2 - cc_dist)/(r1 + r2)*(x1 - x2);
		PD.y_new = y2 + (r1 + r2 - cc_dist)/(r1 + r2)*(y1 - y2);
		PD.z_new = z2 + (r1 + r2 - cc_dist)/(r1 + r2)*(z1 - z2);
		PD.sintering_flag = 1;
	}
	if (sinter_flag == 1)
	{

		cc_dist = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
		r_c = r2*sin(acos(cc_dist/(2*r2)));
		P = F/(4.0*atan(1)*pow(r_c, 2)); // Compression pressure
		A_a = 4*atan(1)*pow(r_c, 2); // Total contact area
		A_r = A_a; // Contact point area
		A_v = A_a - A_r; // Void contact area
		PD.sintering_flag = 1;
	}

	// Material properties or whatever (used for sintering calculations)
	float rho = 7800; // Material density
	float k_n = 700000; // Contact stiffness (mechanical)
	float epsilon = 0.7; // Damping coefficient
	float miu = 0.05; // Friction coefficient
	float D_gDelta_g = 1.8*pow(10, -19); // Diffusion coefficient
	float omega = 1.21*pow(10, -29); // Atomic volume
	float gamma_s = 1.58; // Surface energy
	float psi = (5.0/6.0)*4.0*atan(1); // Dihedral angle
	float k_b = 1.38064852*pow(10, -23); // Boltzmann constant

	float D_b = (D_gDelta_g*omega)/(k_b*T1); // Effective diffusion coefficient
	float v_n = (8*D_b*gamma_s)/(pow(r_c, 4))*(4*r_eq*(1 - cos(psi/2)) + r_c*sin(psi/2));
	// float F_n = (4.0*atan(1)*pow(r_c, 4))/(8*D_b)*v_n - 4.0*atan(1)*gamma_s*(4*r_eq*(1 - cos(psi/2)) + r_c*sin(psi/2));
	float r_c_dot = -r_eq*v_n/r_c; // temporal change in r_c

	r_c = r_c + r_c_dot*delta_t;
	A_a = 4*atan(1)*pow(r_c, 2); // Total contact area
	A_r = A_a; // Contact point area
	A_v = A_a - A_r; // Void contact area
	cc_dist = r1*cos(asin(r_c/r1)) + r2*cos(asin(r_c/r2));
	PD.x_new = x2 + (r1 + r2 - cc_dist)/(r1 + r2)*(x1 - x2);
	PD.y_new = y2 + (r1 + r2 - cc_dist)/(r1 + r2)*(y1 - y2);
	PD.z_new = z2 + (r1 + r2 - cc_dist)/(r1 + r2)*(z1 - z2);

	float h_e = (1.55*m*k/sigma)*pow((P*1.41)/(E_prime*m) , 0.94);
	float h_c = h_e + (A_v/(A_a*sigma))*k_air;
	float G = h_c*A_a; // Overall heat transfer (in room temperature)
	PD.K = G;
	return PD;
}
