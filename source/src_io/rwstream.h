#ifndef RWSTREAM_H
#define RWSTREAM_H

#include <stdio.h>
#include <cstdlib>
#include <complex>
#include <iostream>
using namespace std;


//A class to read or write binary data.
//By qianrui 2020-1-6
class Rwstream
{
	public:
		Rwstream(){
			fileptr=NULL;
		};
		Rwstream(const string,const char*);
		Rwstream(FILE* &ptr);
		~Rwstream(){};
		void setptr(FILE* &ptr);
		FILE* fileptr;
		void close();
		void open(const string,const char*);
		bool operator!() const;
		operator bool() const;
};


//For template operator, you had better write defination in .h file!!
//for variation or array
template<class T>
Rwstream& operator<<(Rwstream& wstream,const T data)
{
    int size=sizeof(T);
	int n=sizeof(data)/sizeof(T);
    fwrite(&data,size,n,wstream.fileptr);
    return wstream;
}
/*//for dynamic memory
//malloc_usable_size has problem!
template<class T>
Rwstream& operator<<(Rwstream& wstream,T* &data)
{
    int size=sizeof(T);
	int n=malloc_usable_size(data)/sizeof(T);
    fwrite(data,size,n,wstream.fileptr);
    return wstream;
}*/
template<class T>
Rwstream& operator>>(Rwstream& rstream,T& data)
{
	int size=sizeof(T);
	int n=sizeof(data)/sizeof(T);
	size_t ch;
    ch=fread(&data,size,n,rstream.fileptr);
	if(ch<n)
	{
		cout<<"Error in Rwstream: Some data didn't be read."<<endl;
    	exit(0);
	}
	return rstream;
}


/*//for dynamic memory
template<class T>
Rwstream& operator>>(Rwstream& rstream,T* &data)
{
	int size=sizeof(T);
	int n=malloc_usable_size(data)/sizeof(T);
	cout<<malloc_usable_size(data)<<' '<<sizeof(T)<<' '<<n<<endl;
	size_t ch;
    ch=fread(data,size,n,rstream.fileptr);
	if(ch<n)
	{
		cout<<"Error in Rwstream: Some dynamic memory didn't be read."<<endl;
		exit(0);
    }
	return rstream;
}*/
template<class T>
void rwread(Rwstream& rstream,T* &data,int n)
{
	int size=sizeof(T);
	size_t ch;
    ch=fread(data,size,n,rstream.fileptr);
    if(ch<n)
    {
        cout<<"Error in Rwstream: Some dynamic memory didn't be read."<<endl;
        exit(0);
    }
    return;
}

#endif

