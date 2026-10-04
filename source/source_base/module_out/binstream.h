#ifndef RWSTREAM_H
#define RWSTREAM_H

#include <stdio.h>
#include <cstdlib>
#include <complex>
#include <iostream>

#include "source_base/tool_quit.h"

/**
 * @brief A stream to read or write binary data.
 * @author qianrui 2020-1-6
 */
class Binstream
{
    public:
        Binstream() = default;
        Binstream(const std::string, const char*);
        ~Binstream();
        // The destructor closes the owned FILE*, so copying would double-close.
        Binstream(const Binstream&) = delete;
        Binstream& operator=(const Binstream&) = delete;
        FILE* fileptr = nullptr;
        void close();
        void open(const std::string, const char*);
        bool operator!() const;
        operator bool() const;

        template<class T>
        Binstream& operator>>(T& data);

        template<class T>
        Binstream& operator<<(const T& data);

        template<class T>
        Binstream& read(T* data, const int n);

        template<class T>
        Binstream& write(const T* data, const int n);
};

// read a data from file
template<class T>
Binstream& Binstream::operator>>(T& data)
{
    const int size = sizeof(T);
    if(this->fileptr == NULL)
    {
        ModuleBase::WARNING_QUIT("Binstream::operator>>",
            "cannot read from an unopened file.");
    }
    size_t ch = fread(&data, size, 1, this->fileptr);
    if(ch < 1)
    {
        ModuleBase::WARNING_QUIT("Binstream::operator>>",
            "Some data couldn't be read. Please make sure you are using op: \"r\".");
    }
    return *this;
}

// write a data into file
template<class T>
Binstream& Binstream::operator<<(const T& data)
{
    const int size = sizeof(T);
    if(this->fileptr == NULL)
    {
        ModuleBase::WARNING_QUIT("Binstream::operator<<",
            "cannot write to an unopened file.");
    }
    size_t ch = fwrite(&data, size, 1, this->fileptr);
    // fwrite() reports success while data is still buffered; fflush() forces
    // the buffer out so delayed write errors (e.g. RLIMIT_FSIZE, full disk)
    // are detected here instead of being silently dropped by fclose().
    if(ch < 1 || fflush(this->fileptr) != 0)
    {
        ModuleBase::WARNING_QUIT("Binstream::operator<<",
            "Some data couldn't be written.");
    }
    return *this;
}

// read an array of data
template<class T>
Binstream& Binstream::read(T* data, const int n)
{
    const int size = sizeof(T);
    if(this->fileptr == NULL)
    {
        ModuleBase::WARNING_QUIT("Binstream::read",
            "cannot read from an unopened file.");
    }
    size_t ch = fread(data, size, n, this->fileptr);
    if(ch < static_cast<size_t>(n))
    {
        ModuleBase::WARNING_QUIT("Binstream::read",
            "Some array elements couldn't be read. Please make sure you are using op: \"r\".");
    }
    return *this;
}

// write an array of data
template<class T>
Binstream& Binstream::write(const T* data, const int n)
{
    const int size = sizeof(T);
    if(this->fileptr == NULL)
    {
        ModuleBase::WARNING_QUIT("Binstream::write",
            "cannot write to an unopened file.");
    }
    size_t ch = fwrite(data, size, n, this->fileptr);
    // fwrite() reports success while data is still buffered; fflush() forces
    // the buffer out so delayed write errors (e.g. RLIMIT_FSIZE, full disk)
    // are detected here instead of being silently dropped by fclose().
    if(ch < static_cast<size_t>(n) || fflush(this->fileptr) != 0)
    {
        ModuleBase::WARNING_QUIT("Binstream::write",
            "Some array elements couldn't be written.");
    }
    return *this;
}

#endif
