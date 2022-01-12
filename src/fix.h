/* -*- C++ -*------------------------------------------------------------
 *
 *	fix.h
 *
 * ---------------------------------------------------------------------- */

#ifndef FIX_H
#define FIX_H

#include "error.h"
#include "system.h"
#include "culsm.h"

class Fix
{
public:
	Fix(Error *error_);
	Fix(Error *error_, int type_, double dispx_, double dispy_, double dispz_);
	~Fix();
	
	int type;
	double dispx, dispy, dispz;
protected:
	Error *error;

private:
	
};


#endif // !FIX_H

