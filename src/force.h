/* -*- C++ -*------------------------------------------------------------
 *
 *	force.h
 *
 * ---------------------------------------------------------------------- */

#ifndef FORCE_H
#define FORCE_H

#include "error.h"
#include "system.h"
#include "culsm.h"

class Force {
public:
	Force(Error *error_);
	~Force();

protected:
	Error *error;
private:
	
};

#endif // !FORCE_H