#include <stdio.h>
#include <string>
#include "source_base/module_out/binstream.h"
#include "source_base/tool_quit.h"

namespace
{
// Binstream is always a *binary* stream. On Windows, fopen mode "r"/"w"/"a"
// opens in text mode, which translates CRLF and treats 0x1A as EOF, corrupting
// binary data (e.g. wavefunction / charge files) -> "Some data couldn't be read".
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
Binstream::Binstream(const std::string filename, const char *op)
{
    fileptr = fopen(filename.c_str(), ensure_binary_mode(op).c_str());
}

Binstream::~Binstream()
{
    if(fileptr != NULL)
    {
        // A destructor must not terminate the program, so a delayed write
        // failure surfacing here can only be reported, not quit on. Callers
        // that need the guarantee should use close(), which checks the result.
        if(fclose(fileptr) != 0)
        {
            ModuleBase::WARNING("Binstream::~Binstream",
                "fclose failed: some buffered data may not have been written.");
        }
        fileptr = NULL;
    }
}

// close file
void Binstream::close()
{
    if(fileptr == NULL)
    {
        return;
    }
    // fclose() flushes the buffer; its failure is the last chance to detect
    // delayed write errors (e.g. RLIMIT_FSIZE, full disk).
    if(fclose(fileptr) != 0)
    {
        fileptr = NULL;
        ModuleBase::WARNING_QUIT("Binstream::close",
            "fclose failed: some buffered data may not have been written.");
    }
    fileptr = NULL;
    return;
}

// open a file
void Binstream::open(const std::string filename, const char *op)
{
    // Close any previously opened file first; overwriting fileptr directly
    // would leak the old handle together with its buffered data.
    close();
    fileptr = fopen(filename.c_str(), ensure_binary_mode(op).c_str());
}

// ! operator
// we can use if(!Binstream) ...
bool Binstream::operator!() const
{
    if (fileptr == NULL)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// bool operator
// we can use if(Binstream) ...
Binstream::operator bool() const
{
    if (fileptr == NULL)
    {
        return false;
    }
    else
    {
        return true;
    }
}
