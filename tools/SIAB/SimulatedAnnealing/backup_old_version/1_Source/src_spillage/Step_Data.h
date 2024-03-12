#ifndef STEP_DATA_H
#define STEP_DATA_H

#include "common.h"

class Step_Data
{
	public:
	// be given value in SpillageStep::init_QS_matrix()
	ComplexMatrix *Qoverlap;
	ComplexMatrix *Soverlap;
	ComplexMatrix *Qtrial;
	ComplexMatrix *Strial;
	ComplexMatrix *Sinv;
	Inverse_Matrix inverse;

    public:
    ComplexArray Qin;

	double *Mk;//[nks]; mohan add 2010-05-02
	matrix Mkb;//[nks][nbands]; mohan add 2010-05-02
	matrix Mkb_used;
	double *weight; //[nks]
    int nks;
    int nbands;
    int nwfc;
    int ne;

	int nwfc2;
	int nwfc_all;
	int test;
	double spillage0;

    Step_Data();
    ~Step_Data();
	void init( const int &nks_in, const double *weight_in, 
			const int &nbands_in, const int &ne_in, 
			const int &nwfc_in, const int &nwfc2_in,
			const int &nwfc_all);

	void initQS();

	private:
    void allocate();
};


#endif
