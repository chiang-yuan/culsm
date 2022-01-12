#include "dump.h"

Dump::Dump(Error *error_)
{
	error = error_;
}

Dump::~Dump()
{
}

int Dump::set(char* filepattern_, int ndump_)
{
	filepattern = new char[MAX_STRING];
	strcpy(filepattern, filepattern_);
	
    ndump = ndump_;

	return 0;
}

int Dump::write_lmpdump(int timestep, System & sys)
{
	char* filename = new char[MAX_STRING];
	sprintf(filename, filepattern, timestep);

	file = fopen(filename, "a");
	if (file == NULL) return error->message("Cannot open file \"%s\"", -1, filename);

	// printf("Write into lammps dump file: %s ...\n", filename);

	// Section: System

	fprintf(file, "ITEM: TIMESTEP\n%d\n", timestep);
	fprintf(file, "ITEM: NUMBER OF ATOMS\n%lld atoms\n", sys.natoms);

	for (bigint a = 0; a < sys.natoms; a++) {
		bigint i3 = a * 3;
		sys.box[0][0] = std::min(sys.box[0][0], sys.x[i3]);
		sys.box[0][1] = std::max(sys.box[0][1], sys.x[i3]);
		sys.box[1][0] = std::min(sys.box[1][0], sys.x[i3 + 1]);
		sys.box[1][1] = std::max(sys.box[1][1], sys.x[i3 + 1]);
		sys.box[2][0] = std::min(sys.box[2][0], sys.x[i3 + 2]);
		sys.box[2][1] = std::max(sys.box[2][1], sys.x[i3 + 2]);
	}
	fprintf(file, "ITEM: BOX BOUNDS ss ss pp\n%23.16e %23.16e\n%23.16e %23.16e\n%23.16e %23.16e\n",
		sys.box[0][0], sys.box[0][1],
		sys.box[1][0], sys.box[1][1],
		sys.box[2][0], sys.box[2][1]
	);

	// Section : Atoms

	fprintf(file, "ITEM: ATOMS id type x y z sxx syy szz sxy syz szx\n");

	for (bigint a = 0; a < sys.natoms; a++) {
		bigint i3 = a * 3;
		bigint i6 = a * 6;
		fprintf(file, "%6lld %3d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
			a + 1, sys.type[a],
			sys.x[i3], sys.x[i3 + 1], sys.x[i3 + 2],
			h_stress[i6], h_stress[i6 + 1], h_stress[i6 + 2], 
			h_stress[i6 + 3], h_stress[i6 + 4], h_stress[i6 + 5]
		);
		if ( a % MAX_CHUNK == 0 ) fflush(file);
	}

	fclose(file);

	return 0;
}


