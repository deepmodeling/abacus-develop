#include <stdio.h>
#include <string>
#include "binstream.h"

namespace
{
// Binstream is always a *binary* stream. On Windows, fopen mode "r"/"w"/"a"
// opens in text mode, which translates CRLF and treats 0x1A as EOF, corrupting
// binary data (e.g. wavefunction / charge files) -> "Some data didn't be read".
// Append "b" if the caller didn't, so binary mode is always used. On POSIX the
// "b" flag is a harmless no-op, so the Linux behaviour is unchanged.
std::string ensure_binary_mode(const char* op)
{
	std::string mode(op ? op : "");
	if (mode.find('b') == std::string::npos)
	{
		mode += 'b';
	}
	return mode;
}
} // namespace

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
	fileptr=fopen(filename.c_str(),ensure_binary_mode(op).c_str());
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
	fileptr=fopen(filename.c_str(),ensure_binary_mode(op).c_str());
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
