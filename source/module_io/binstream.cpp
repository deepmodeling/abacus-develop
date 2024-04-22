#include "binstream.h"

#include <cstdio>

/**
 * @brief Construct a new Binstream:: Binstream object
 * 
 * @param filename 
 * @param op "r": read
 *           "a": add
 *           "w": write 
 */
Binstream::Binstream(const std::string& filename,const char *op)
{
	fileptr=fopen(filename.c_str(),op);
}

Binstream::~Binstream()
{
	if(fileptr != nullptr)	fclose(fileptr);
}

// close file
void Binstream:: close()
{
	fclose(fileptr);
	fileptr = nullptr;
}

// open a file
void Binstream::open(const std::string& filename,const char *op)
{
	fileptr=fopen(filename.c_str(),op);
}

// ! operator
// we can use if(!Binstream) ...
bool Binstream::operator!() const
{
	return fileptr==nullptr;
}

// bool operator
// we can use if(Binstream) ...
Binstream::operator bool() const
{
	return fileptr != NULL;
}
