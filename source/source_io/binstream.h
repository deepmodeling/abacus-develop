#ifndef RWSTREAM_H
#define RWSTREAM_H

#include <stdio.h>
#include <cstdlib>
#include <complex>
#include <iostream>


/**
 * @brief A stream to read or write binary data.
 * @author qianrui 2020-1-6
 */
class Binstream
{
	public:
		Binstream(){
			fileptr=NULL; //we should use NULL (not nullptr) here because FILE use NULL.
		};
		Binstream(const std::string,const char*);
		~Binstream();
		FILE* fileptr;
		void close();
		void open(const std::string,const char*);
		bool operator!() const;
		operator bool() const;

		template<class T>
		Binstream& operator>>( T& data);

		template<class T>
		Binstream& operator<<(const T& data);

		template<class T>
		Binstream& read(T* data,const int n);

		template<class T>
		Binstream& write(const T* data,const int n);

};

// read a data from file
template<class T>
Binstream& Binstream:: operator>>(T& data)
{
	const int size=sizeof(T);
	size_t ch = fread(&data,size,1,this->fileptr);
	if(ch<1)
	{
		std::cout<<"Error in Binstream: Some data didn't be read."<<std::endl;
		std::cout<<"Please make you are using op: \"r\""<<std::endl;
    	exit(0);
	}
	return *this;
}

// write a data into file
template<class T>
Binstream& Binstream:: operator<<(const T& data)
{
    const int size=sizeof(T);
    fwrite(&data,size,1,this->fileptr);
    return *this;
}

//read an array of data
template<class T>
Binstream& Binstream::read(T* data, const int n)
{
	const int size=sizeof(T);
	size_t ch = fread(data,size,n,this->fileptr);
    if(ch<n)
    {
        std::cout<<"Error in Binstream: Some dynamic memory didn't be read."<<std::endl;
		std::cout<<"Please make you are using op: \"r\""<<std::endl;
        exit(0);
    }
    return *this;
}

//write an array of data
template<class T>
Binstream& Binstream::write(const T* data, const int n)
{
	const int size=sizeof(T);
    fwrite(data,size,n,this->fileptr);
    return *this;
}



/*//for dynamic memory
//malloc_usable_size has problem!
template<class T>
Binstream& operator<<(Binstream& wstream,T* &data)
{
    int size=sizeof(T);
	int n=malloc_usable_size(data)/sizeof(T);
    fwrite(data,size,n,wstream.fileptr);
    return wstream;
}


//for dynamic memory
template<class T>
Binstream& operator>>(Binstream& rstream,T* &data)
{
	int size=sizeof(T);
	int n=malloc_usable_size(data)/sizeof(T);
	std::cout<<malloc_usable_size(data)<<' '<<sizeof(T)<<' '<<n<<std::endl;
	size_t ch;
    ch=fread(data,size,n,rstream.fileptr);
	if(ch<n)
	{
		std::cout<<"Error in Binstream: Some dynamic memory didn't be read."<<std::endl;
		exit(0);
    }
	return rstream;
}
*/

#endif

