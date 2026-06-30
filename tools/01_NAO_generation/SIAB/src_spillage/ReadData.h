#ifndef ReadData_H
#define ReadData_H

#include "common.h"

class ReadData
{
	public:

	ComplexArray Qin;
	ComplexArray *Sq1q2;
	ComplexArray Sinv;

	//===================================
	// read in data Q=<J|psi> , S=<J|J>
	//===================================
	void OverlapQandS( const string &name);

	//========================================================
	// inverse matrix of S directly, 
	// used for checking the correctness
	// of spillage calculation
	//========================================================
	void OverlapSinv( const string &name);

	double *weight;

	private:
	void OverlapQ( ifstream &ifs);
	void OverlapSq1q2( ifstream &ifs );

	int test;

	public:
	ReadData();
	~ReadData();

};

#endif
