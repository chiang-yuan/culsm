/* -*- C++ -*------------------------------------------------------------
 *
 *	run.h
 *
 * ---------------------------------------------------------------------- */

#ifndef RUN_H
#define RUN_H

#include "culsm.h"
#include "error.h"
#include "system.h"
#include "dump.h"
#include "thermo.h"
#include "fix.h"

class Run {
public:
	Run(Error *error_);
	~Run();
	
	virtual int verlet(
		float dt, int timesteps, 
		System &sys, Dump &dump, Thermo &thermo, std::vector<Fix> &fixes);
protected:
	Error *error;
	
private:
	
};

#endif // !RUN_H

