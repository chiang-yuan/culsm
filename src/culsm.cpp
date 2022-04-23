/* -*- C++ -*------------------------------------------------------------
 *
 *	culsm.cpp
 *
 * ---------------------------------------------------------------------- */

#include "culsm.h"
#include "system.h"
#include "error.h"
#include "read_data.h"
#include "run.h"
#include "dump.h"
#include "thermo.h"
#include "fix.h"
// #include "write_data.h"

using namespace std;

int parse_command(const char* line) {
    string clone = line;

    // printf("clone: %s\n", clone.c_str());

    // trim off comment after character #
    string trimmed = clone.substr(0, clone.find("#", 0));

    // printf("trimmed: %s\n", trimmed.c_str());

    char* token;
    char *trimmed_cstr = new char[strlen(trimmed.c_str()) + 1];
	strcpy(trimmed_cstr, trimmed.c_str());

    int n = 0; 

    if ( (token = strtok(trimmed_cstr, WHITESPACE)) == NULL) {
		return n;
	}

    while (token != NULL) {
        // printf("(%i) %s\t", ++n, token);
        token = strtok(NULL, WHITESPACE);
        n++;
    }

    return n;
}

int main(int argc, char *argv[]) {
	time_t tstr, tend;
	time(&tstr);
	printf("=====================================================================\n");
	printf("CULSM (%s) starts at %s", CULSM_VERSION, ctime(&tstr));
	printf("=====================================================================\n");
	

	Error error;
	System sys(&error);
	ReadData reader(&error);
	Run run(&error);
	Dump dump(&error);
	Thermo thermo(&error);
	vector<Fix> fixes;
	
	// WriteData writer(&error);

    char* command = new char[MAX_STRING];

    int i = 0;
    while ( fgets(command, MAX_STRING, stdin) != NULL) {

        argc = parse_command(command);

        if (argc == 0) continue;

        char* buffer = new char[MAX_STRING];
        sprintf(buffer, "[command %2d]: ", ++i);
        puts(strcat(buffer, command));

		strcpy(buffer, command);
		argv = new char*[argc];
		argv[0] = strtok(buffer, WHITESPACE);

        if (argv[0] == NULL) {
            // TODO: error 
			return error.message("Incorrect input commnad line", -1);
		}

		for (int j = 1; j < argc; j++) {
			argv[j] = strtok(NULL, WHITESPACE);
			if (argv[j] == NULL) {
                // TODO: error
				return error.message("Missing argument in command line", -1);
			}
		}

        // TODO: commend list
        if (strncmp(argv[0], "read_data", 9) == 0) {
            int nargs = 1;
            char** args = new char*[nargs];

            if (argc != 1 + nargs) error.message("Incorrect argument for read_data", -1);
            
            for (int j = 0; j < nargs; j++) args[j] = argv[1+j];

			reader.load_from_lmpdata(args[0], sys);
        }

		if (strncmp(argv[0], "run", 3) == 0) {
			float dt;
			int timesteps;
			
			if (sscanf(command, "%*s %f %d", &dt, &timesteps) == EOF) error.message("Incorrect argument for run", -1);
			
			run.verlet(dt, timesteps, sys, dump, thermo, fixes);
		}

		if (strncmp(argv[0], "fix", 3) == 0) {
			int type = atoi(argv[1]);
			double dispx = atof(argv[2]);
			double dispy = atof(argv[3]);
			double dispz = atof(argv[4]);
		
			// if (sscanf(command, "%*s %d %lf %lf %lf", &type, &dispx, &dispy, &dispz) == EOF) error.message("Incorrect argument for dump", -1);

			fixes.push_back(Fix(&error, type, dispx, dispy, dispz));
		}

		if (strncmp(argv[0], "dump", 4) == 0) {
			char* filepattern = new char[MAX_STRING];
			int ndump;
		
			if (sscanf(command, "%*s %s %d", filepattern, &ndump) == EOF) error.message("Incorrect argument for dump", -1);

			dump.set(filepattern, ndump);
		}

		if (strncmp(argv[0], "thermo", 4) == 0) {
			int nthermo;
		
			if (sscanf(command, "%*s %d", &nthermo) == EOF) error.message("Incorrect argument for thermo", -1);

			thermo.set(nthermo);
		}

		if (strncmp(argv[0], "mass", 4) == 0) {

			if (argc < 3) error.message("Incorrect argument for mass", -1);

			int atype = atoi(argv[1]);  // atom type (1 based)
			double m = atof(argv[2]);

			sys.atomTypes[atype - 1].mass = m;	
		}

		if (strncmp(argv[0], "bond", 4) == 0) {

			if (argc < 2) error.message("Incorrect argument for bond", -1);

			int btype = atoi(argv[1]);  // bond type (1 based)

			if (strncmp(argv[2], "leb", 3) == 0) {

				if (argc != 1 + 5) error.message("Incorrect argument for linear elastic brittle bond", -1);
				double ke = atof(argv[3]);
				double r0 = atof(argv[4]);
				double rc = atof(argv[5]);

				sys.bondTypes[btype - 1].no_bond_coeffs = 3;
				
				sys.bondTypes[btype - 1].coeff[0] = ke;	// linear elastic stiffness
				sys.bondTypes[btype - 1].coeff[1] = r0;	// equilibrium length
				sys.bondTypes[btype - 1].coeff[2] = rc;	// critical length

				strcpy(sys.bondTypes[btype - 1].name, argv[2]);
			}
			else if (strncmp(argv[2], "lecs", 4) == 0) {

				if (argc != 1 + 4) error.message("Incorrect argument for linear elastic critical strain bond", -1);
				double ke = atof(argv[3]);
				double cs = atof(argv[4]);

				sys.bondTypes[btype - 1].no_bond_coeffs = 2;
				
				sys.bondTypes[btype - 1].coeff[0] = ke;	// linear elastic stiffness
				sys.bondTypes[btype - 1].coeff[1] = cs;	// critical strain

				strcpy(sys.bondTypes[btype - 1].name, argv[2]);
			}
			else if (strncmp(argv[2], "bep", 3) == 0) {

				if (argc != 1 + 7) error.message("Incorrect argument for linear elastic critical strain bond", -1);
				double ke = atof(argv[3]);
				double kp = atof(argv[4]);
				double r0 = atof(argv[5]);
				double rc = atof(argv[6]);
				double fy = atof(argv[7]);

				sys.bondTypes[btype - 1].no_bond_coeffs = 5;
				
				sys.bondTypes[btype - 1].coeff[0] = ke;
				sys.bondTypes[btype - 1].coeff[1] = kp;
				sys.bondTypes[btype - 1].coeff[2] = r0;
				sys.bondTypes[btype - 1].coeff[3] = rc;
				sys.bondTypes[btype - 1].coeff[4] = fy;

				strcpy(sys.bondTypes[btype - 1].name, argv[2]);
			}
		}
    }
	
	time(&tend);
	int seconds = (int)difftime(tend, tstr);
	printf("Total time: %d:%02d:%02d\n", seconds/3600, (seconds%3600)/60, seconds%60);
	printf("=====================================================================\n");
	printf("CULSM (%s) ends at %s", CULSM_VERSION, ctime(&tend));
	printf("=====================================================================\n");
    return 0;
}


/*
int count_words(const char* line) {
	int n = strlen(line) + 1;
	char *copy = new char[n];
	strcpy(copy, line);

	char *ptr;
	if ((ptr = strchr(copy, '#'))) *ptr = ' ';

	if (strtok(copy, WHITESPACE) == NULL) {
		delete[] copy;
		return 0;
	}

	n = 1;
	while (strtok(NULL, WHITESPACE)) n++;
	delete[] copy;

bool check_arg(char **arg, const char *flag, int num, int argc) {

	if (num >= argc) {
		printf("Missing argument for \"%s\" flag\n", flag);
		return false;
	}

	if (arg[num][0] == '-' || arg[num][0] == '<' || arg[num][0] == '>') {
		printf("Incorrect argument to \"%s\" flag: %s", flag, arg[num]);
		return false;
	}
	return true;
}

int main(int argc, char *argv[]) {

	// system

	Error error;
	System sys(&error);
	// ReadData reader(&error);
	// WriteData writer(&error);

	if (argc < 2) {
		//printf("usage command: < (input file name) [-w delete type #] [-w add type #] [-s number (type#) #] [-h valence n/# (typeO) (typeH) (typeO-H)] [> (output file name) hint y/n]\n");
		char* buffer = new char[MAX_STRING];
		fgets(buffer, MAX_STRING, stdin);

		argc = count_words(buffer);
		argv = new char*[argc];

		argv[0] = strtok(buffer, WHITESPACE);
		if (argv[0] == NULL) {
			return error.message("Incorrect input commnad line", 1);
		}

		for (int i = 1; i < argc; i++) {
			argv[i] = strtok(NULL, WHITESPACE);
			if (argv[i] == NULL) {
				return error.message("Missing argument in command line", 2);
			}
		}

		int n = 0;
		while (n < argc)
		{
			if (strncmp(argv[n], "<", 1) == 0) {
				int ncomm = 1;
				char** commd = new char*[ncomm];

				for (int i = 0; i < ncomm; i++) {
					n++;
					if (!check_arg(argv, "read_data", n, argc)) return error.message("", 3);
					commd[i] = argv[n];
				}
				// printf("ReadData::command(): %d\n", reader.command(ncomm, commd, sys));

			}

			if (strncmp(argv[n], ">", 1) == 0) {
				int ncomm = 3;
				char** commd = new char*[ncomm];

				for (int i = 0; i < ncomm; i++) {
					n++;
					if (!check_arg(argv, "write_data", n, argc)) return error.message("", 6);
					commd[i] = argv[n];
				}

				// printf("WriteData::commd(): %d\n", writer.command(ncomm, commd, sys));
			}

			if (argv[n] == NULL) {
				return error.message("Missing argument in command line", 5);
			}
			n++;
		}
		system("pause");
	}
	else {
		time_t localtime;
		time(&localtime);
		printf("=====================================================================\n");
		printf("CULSM (%s) starts at %s", CULSM_VERSION, ctime(&localtime));
		printf("=====================================================================\n");

		int n = 1;
		while (n < argc)
		{
			if (strncmp(argv[n], "<", 1) == 0) {
				int ncomm = 1;
				char** commd = new char*[ncomm];

				for (int i = 0; i < ncomm; i++) {
					n++;
					if (!check_arg(argv, "read_data", n, argc)) return error.message("", 3);
					commd[i] = argv[n];
				}
				// printf("End ReadData: %d\n", reader.command(ncomm, commd, sys));

			}


			if (strncmp(argv[n], ">", 1) == 0) {
				int ncomm = 3;
				char** commd = new char*[ncomm];

				for (int i = 0; i < ncomm; i++) {
					n++;
					if (!check_arg(argv, "write_data", n, argc)) return error.message("", 6);
					commd[i] = argv[n];
				}

				// printf("End WriteData: %d:\n", writer.command(ncomm, commd, sys));
			}

			if (argv[n] == NULL) {
				return error.message("Missing argument in command line", 5);
			}
			n++;
		}
		time(&localtime);
		printf("=====================================================================\n");
		printf("CULSM (%s) ends at %s", CULSM_VERSION, ctime(&localtime));
		printf("=====================================================================\n");
	}

	return 0;
}
*/