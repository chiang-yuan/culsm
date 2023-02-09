#include "create_bonds.h"

// __global__ void add_range_type(
// 	int* type, int* x, 
// 	int ati, int atj, float rmin, float rmax,
// 	int natoms
// )
// {
// 	extern __shared__ int nbonds;
// 	bigint ai = blockIdx.x * blockDim.x + threadIdx.x;

// 	if (ai < natoms) {

// 		for (bigint aj = ai+1; aj < natoms; aj++) {

// 			if (!((type[ai] == ati && type[aj] == atj) ||
// 				(type[ai] == atj && type[aj] == ati))) 
// 				continue;

// 			int ai3 = ai*3;
// 			int aj3 = aj*3;

// 			double rsq = pow(x[ai3] - x[aj3], 2) +
// 				pow(x[ai3 + 1] - x[aj3 + 1], 2) +
// 				pow(x[ai3 + 2] - x[aj3 + 2], 2);

// 			if (rsq < pow(rmin, 2) || rsq > pow(rmax, 2)) 
// 				continue;

// 			if 

// 		}
// 	}

// }

BondCreator::BondCreator(Error *error_)
{
	error = error_;
}

BondCreator::~BondCreator()
{
}

int BondCreator::range_type(
	int btype, // types of added bonds
	int ati, int atj, // types of i, j atoms 
	float rmin, float rmax, // max and min of added range
	System& sys
	) {

	int* new_bond_type = NULL;
	int* new_atom_i = NULL;
	int* new_atom_j = NULL;
	
	for (bigint ai = 0; ai < sys.natoms; ai++) {
		for (bigint aj = ai+1; aj < sys.natoms; aj++) {

			if (ai == aj) continue;

			if (!((sys.type[ai] == ati && sys.type[aj] == atj) ||
				(sys.type[ai] == atj && sys.type[aj] == ati))) 
				continue;

			int ai3 = ai*3;
			int aj3 = aj*3;

			double rsq = pow(sys.x[ai3] - sys.x[aj3], 2) +
				pow(sys.x[ai3 + 1] - sys.x[aj3 + 1], 2) +
				pow(sys.x[ai3 + 2] - sys.x[aj3 + 2], 2);

			if (rsq < pow(rmin, 2) || rsq > pow(rmax, 2)) 
				continue;

			if (sys.nbonds == 0) {
				sys.nbonds++; 
				sys.bond_type = (int*) malloc(sys.nbonds*sizeof(int));
				sys.atom_i = (int*) malloc(sys.nbonds*sizeof(int));
				sys.atom_j = (int*) malloc(sys.nbonds*sizeof(int));

				sys.bond_type[sys.nbonds - 1] = btype;
				sys.atom_i[sys.nbonds - 1] = ai;
				sys.atom_j[sys.nbonds - 1] = aj;
			}
			else {
				sys.nbonds++;
				new_bond_type = (int*) realloc(sys.bond_type, sys.nbonds*sizeof(int));
				new_atom_i = (int*) realloc(sys.atom_i, sys.nbonds*sizeof(int));
				new_atom_j = (int*) realloc(sys.atom_j, sys.nbonds*sizeof(int));

				if (new_bond_type == NULL || new_atom_i == NULL || new_atom_j == NULL) {
					free(sys.bond_type);
					free(sys.atom_i);
					free(sys.atom_j);

					return error->message("Failed to reallocate memory", 1);
				} 
				else {
					sys.bond_type = new_bond_type;
					sys.atom_i = new_atom_i;
					sys.atom_j = new_atom_j;

					sys.bond_type[sys.nbonds - 1] = btype;
					sys.atom_i[sys.nbonds - 1] = ai;
					sys.atom_j[sys.nbonds - 1] = aj;
				}
				
			}

		}
	}

	return 0;
}

