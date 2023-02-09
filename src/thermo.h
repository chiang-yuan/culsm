/* -*- C++ -*------------------------------------------------------------
 *
 *	thermo.h
 *
 * ---------------------------------------------------------------------- */

#ifndef THERMO_H
#define THERMO_H

#include "culsm.h"
#include "error.h"
#include "system.h"


// host vectors

// extern float* h_k;			// bond stiffness
// extern double* h_parsum_pe;

// device vectors

extern float* d_m;
extern double* d_x;
extern double* d_v;
extern float* d_ke;
extern float* d_r0;
extern float* d_ue;
extern float* d_fa;
extern int* d_atom_i;
extern int* d_atom_j;
extern float* d_stress;

// extern double* d_parsum_pe;

// extern __global__ void reduce_pe(
// 	double* x, float* k, float* r0, int* atom_i, int* atom_j,
// 	double* parsum_pe, // partial sum of potential energy
// 	int nbonds
// );

class Thermo
{
public:
	Thermo(Error *error_);
	~Thermo();
	
	int nthermo;

	int set(int nthermo_);
	double pe(System & sys);	// potential energy
	double ke(System & sys);	// kinetic energy
	double stress(System & sys, int comp);	// stress
	virtual int write_thermo(int timestep, System &sys);
protected:
	char* filepattern;
	FILE* file;
	Error *error;

private:
	
};


#endif // !THERMO_H

