/* -*- C++ -*------------------------------------------------------------
 *
 *	create_bonds.h
 *
 * ---------------------------------------------------------------------- */

#ifndef CREATE_BONDS_H
#define CREATE_BONDS_H

#include "error.h"
#include "system.h"
#include "culsm.h"

class BondCreator
{
public:
	BondCreator(Error *error_);
	~BondCreator();
	
	// int range_many(int bondtype, float rmin, float rmax, System &sys);
	int range_type(int btype, int ati, int atj, float rmin, float rmax, System &sys);
protected:
	Error *error;

private:
	
};


#endif // !CREATE_BONDS_H

