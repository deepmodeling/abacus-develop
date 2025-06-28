#include <stdio.h>
#include "binstream.h"

/**
 * @brief Construct a new Binstream:: Binstream object
 * 
 * @param filename 
 * @param op "r": read
 *           "a": add
 *           "w": write 
 */
Binstream::Binstream(const std::string filename,const char *op)
{
	fileptr=fopen(filename.c_str(),op);
}

Binstream::~Binstream()
{
	if(fileptr != NULL)	fclose(fileptr);
}

// close file
void Binstream:: close()
{
	fclose(fileptr);
	fileptr = NULL;
	return;
}

// open a file
void Binstream::open(const std::string filename,const char *op)
{
	fileptr=fopen(filename.c_str(),op);
}

// ! operator
// we can use if(!Binstream) ...
bool Binstream::operator!() const
{
	if (fileptr==NULL)
		return true;
	else
		return false;
}

// bool operator
// we can use if(Binstream) ...
Binstream::operator bool() const
{
	if (fileptr==NULL)
		return false;
	else
		return true;
}
