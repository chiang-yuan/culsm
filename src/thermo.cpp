#include "thermo.h"

Thermo::Thermo(Error *error_)
{
	error = error_;
}

Thermo::~Thermo()
{
}

int Thermo::set(int nthermo_)
{
    nthermo = nthermo_;
	return 0;
}

double Thermo::pe(System & sys)
{
	// TODO: parallel reduction
	double energy = 0;

	// TODO: other potentials (non-bonded, angles, etc)

	// bonds

	bigint b;
	int btype;
	bigint ai;
	bigint aj;

	double k;
	double r0;
	double rij;	

	for (b = 0; b < sys.nbonds; b++) {
		btype = sys.bond_type[b];
		// k = sys.bondTypes[itype].coeff[0];
		k = h_k[b];
		r0 = sys.bondTypes[btype - 1].coeff[1];
		
		ai = sys.atom_i[b];
		aj = sys.atom_j[b];
		
		rij = sqrt(pow(sys.x[ai*3] - sys.x[aj*3],2) + 
			pow(sys.x[ai*3 + 1] - sys.x[aj*3 + 1],2) +
			pow(sys.x[ai*3 + 2] - sys.x[aj*3 + 2],2));

		energy += 1.0/2.0*k*pow(rij - r0, 2);
	}

	return energy;
}

double Thermo::ke(System & sys)
{
	// TODO: parallel reduction
	double energy = 0;

	bigint a;
	int type;
	double mass;
	double v2;

	for (a = 0; a < sys.natoms; a++) {
		type = sys.type[a];
		mass = sys.atomTypes[type - 1].mass;
		
		v2 = pow(sys.v[a*3], 2) + 
			pow(sys.v[a*3 + 1], 2) +
			pow(sys.v[a*3 + 2], 2);

		energy += 1.0/2.0*mass*v2;
	}

	return energy;
}

int Thermo::write_thermo(int timestep, System & sys)
{
	printf("%10d \t %16.9e \t %16.9e\n", 
		timestep, pe(sys), ke(sys));

	return 0;
}


