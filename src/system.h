/* -*- C++ -*------------------------------------------------------------
 *
 *	system.h
 *
 * ---------------------------------------------------------------------- */

#ifndef SYSTEM_H
#define SYSTEM_H

#define	MAX_ATOMS	100000
#define MAX_BONDS	100000
#define MAX_ANGLES	100000

#include "atom.h"
#include "bond.h"
#include "angle.h"
#include "error.h"
#include "culsm.h"

struct AtomType {
	char element[2];			// elemnet symbol
	double mass;
	char potential[MAX_NAME];	// potential name
	double coeff[2];			// potential coefficient
};

struct BondType {
	int IJType[2];				// types of i and j atom
	double coeff[3];			// bond coefficient (k, r0, rc)
	char style[MAX_NAME];		// bond style
	char name[MAX_NAME];		// bond name
};

struct AngleType {
	int types[3];				
	double coeff[2];			// angle coefficient
	char style[MAX_NAME];		// angle style
	char name[MAX_NAME];		// angle name
};

class System {
public:
	int triclinic_flag;

	double crystal[6];
	double box[3][3];								// xlo xhi yz
													// ylo yhi xz 
													// zlo zhi xy
	enum pbc { p, n, s, m }; pbc periodic[3];		// p: periodic, 
													// n: non-periodic and fixed
													// s: non-periodic and shrink-wrapped
													// m: non-periodic and shrink-wrapped with a minimal value

	bigint natoms;
	bigint nbonds;
	bigint nangles;
	bigint ndihedrals;
	bigint nimpropers;

	int no_atom_types;
	int no_bond_types;
	int no_angle_types;
	int no_dihedral_types;
	int no_improper_types;

	int no_molecular;

	double *x;	// particle coordinates
	double *v;	// particle velocities
	double *a;	// particle acceleration
	int *type;	// particle types

	int *bond_type;	//
	int *atom_i;	// 
	int *atom_j;	//

	std::vector<Atom> atoms;
	std::vector<Bond> bonds;
	std::vector<Angle> angles;

	std::vector<AtomType> atomTypes;
	std::vector<BondType> bondTypes;
	std::vector<AngleType> angleTypes;
	
	friend std::ostream& operator<<(std::ostream& out, const System &sys);
	
	System();
	System(Error *error_);
	~System();
protected:
	Error *error;
private:

};

std::ostream & operator<<(std::ostream& out, const System & sys);

#endif // !SYSTEM_H

