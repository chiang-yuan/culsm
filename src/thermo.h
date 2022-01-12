/* -*- C++ -*------------------------------------------------------------
 *
 *	thermo.h
 *
 * ---------------------------------------------------------------------- */

#ifndef THERMO_H
#define THERMO_H

#include "error.h"
#include "system.h"
#include "culsm.h"

extern float* h_k;			// bond stiffness

class Thermo
{
public:
	Thermo(Error *error_);
	~Thermo();
	
	int nthermo;

	int set(int nthermo_);
	double pe(System & sys);	// potential energy
	double ke(System & sys);	// kinetic energy
	virtual int write_thermo(int timestep, System &sys);
protected:
	char* filepattern;
	FILE* file;
	Error *error;

private:
	
};


#endif // !THERMO_H

