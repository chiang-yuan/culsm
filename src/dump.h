/* -*- C++ -*------------------------------------------------------------
 *
 *	dump.h
 *
 * ---------------------------------------------------------------------- */

#ifndef DUMP_H
#define DUMP_H

#include "error.h"
#include "system.h"
#include "culsm.h"

extern float* h_stress;	// atom stresses

class Dump
{
public:
	Dump(Error *error_);
	~Dump();
	
	int ndump;

	int set(char* filepattern_, int ndump_);
	virtual int write_lmpdump(int timestep, System &sys);
protected:
	char* filepattern;
	FILE* file;
	Error *error;

private:
	
};


#endif // !DUMP_H

