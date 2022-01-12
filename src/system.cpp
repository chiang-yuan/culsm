#include "system.h"

System::System()
{
	natoms = 0;
	nbonds = 0;
	nangles = 0;
	ndihedrals = 0;
	nimpropers = 0;

	no_atom_types = 0;
	no_bond_types = 0;
	no_angle_types = 0;
	no_dihedral_types = 0;
	no_improper_types = 0;
	
	no_molecular = 0;
}

System::System(Error *error_)
{
	error = error_;

	natoms = 0;
	nbonds = 0;
	nangles = 0;
	ndihedrals = 0;
	nimpropers = 0;

	no_atom_types = 0;
	no_bond_types = 0;
	no_angle_types = 0;
	no_dihedral_types = 0;
	no_improper_types = 0;
	
	no_molecular = 0;
}

System::~System()
{
}

std::ostream & operator<<(std::ostream& out, const System & sys)
{
	char* output = new char[MAX_LINE*MAX_STRING];
	char* buffer = new char[MAX_STRING];

	sprintf(buffer, "System information\n");
	strcpy(output, buffer);
	// atom types
	sprintf(buffer, "%ld atom types\n", sys.atomTypes.size());
	strcat(output, buffer);
	for (int i = 0; i < sys.atomTypes.size(); i++) {
		int num = 0;
		for (int a = 0; a < sys.natoms; a++) {
			if (sys.type[a] == i + 1) num++;
		}
		// for (auto a : sys.atoms) {
		// 	if (a.type == i + 1) num++;
		// }
		buffer = new char[MAX_STRING];
		sprintf(buffer, "\ttype %2d %5s : %5d\n",
			i + 1, 
			sys.atomTypes[i].element,
			num);
		strcat(output, buffer);
	}

	// bond types
	sprintf(buffer, "%ld bond types\n", sys.bondTypes.size());
	strcat(output, buffer);
	for (int i = 0; i < sys.bondTypes.size(); i++) {
		int num = 0;
		for (int b = 0; b < sys.nbonds; b++) {
			if (sys.bond_type[b] == i + 1) num++;
		}
		// for (auto b : sys.bonds) {
		// 	if (b.type == i + 1) num++;
		// }
		buffer = new char[MAX_STRING];
		sprintf(buffer, "\ttype %2d %5s : %5d\n",
			i + 1,
			sys.bondTypes[i].name,
			num);
		strcat(output, buffer);
	}

	// // angle types
	// sprintf(buffer, "           %10s           \n--------------------------------\n", "Angle Types");
	// strcat(output, buffer);
	// for (int i = 0; i < sys.angleTypes.size(); i++) {
	// 	int num = 0;
	// 	for (auto a : sys.angles) {
	// 		if (a.type == i + 1) num++;
	// 	}
	// 	buffer = new char[MAX_STRING];
	// 	sprintf(buffer, "%5d %5s : %5d\n",
	// 		i + 1,
	// 		sys.angleTypes[i].name,
	// 		num);
	// 	strcat(output, buffer);
	// }

	std::string str = output;
	out << str;
	return out;
}
