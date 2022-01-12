/* -*- C++ -*------------------------------------------------------------
 *
 *	read_data.h
 *
 * ---------------------------------------------------------------------- */

#ifndef READ_DATA_H
#define READ_DATA_H

#include "error.h"
#include "system.h"
#include "culsm.h"

class ReadData
{
public:
	ReadData(Error *error_);
	~ReadData();
	
	// virtual int command(int argc, char * argv[], System &sys);
	virtual int load_from_lmpdata(char* filename_, System &sys);
protected:
	FILE* file;
	char buffer[MAX_STRING];	// buffer char array to store read line
	char header[MAX_STRING];	// header
	Error *error;

	// bool check_arg(char **arg, const char *flag, int num, int argc);
	int read_block(FILE* file_, int nrows_, int ncolumns_, char *buffer_);
	int count_words(const char* line);	
private:
	
};


#endif // !READ_DATA_H

