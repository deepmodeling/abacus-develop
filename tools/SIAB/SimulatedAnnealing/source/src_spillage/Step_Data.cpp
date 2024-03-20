#include "Step_Data.h"
#include "../src_tools/complexmatrix_inline.h"

Step_Data::Step_Data()
{
//	TITLE("Step_Data","Step_Data");
	Qoverlap = new ComplexMatrix[1];
	Soverlap = new ComplexMatrix[1];
	Qtrial = new ComplexMatrix[1];
	Strial = new ComplexMatrix[1];
	weight = new double[1];
	Sinv = new ComplexMatrix[1];
	inverse_S = new Inverse_Matrix_S[1];
	Mk = new double[1];
	test = 0;
}

Step_Data::~Step_Data()
{
//	TITLE("Step_Data","~Step_Data");
	delete[] Qoverlap;
	delete[] Soverlap;
	delete[] Qtrial;
	delete[] Strial;
	delete[] weight;
	delete[] Sinv;
	delete[] inverse_S;
	delete[] Mk;
}

void Step_Data::init(
	const int &nks_in, 
	const double* weight_in,
	const int &nbands_in, 
	const int &ne_in, 
	const int &nwfc_in,
	const int &nwfc2_in,
	const int &nwfc_all_in)
{
	if(test==1)TITLE("Step_Data","init");
	this->nks = nks_in;
	assert(nks > 0);
	
	delete[] Mk; //mohan add 2010-05-02
	this->Mk = new double[nks];
	ZEROS( Mk, nks );
	
	delete[] weight;
	this->weight = new double[nks];
	for(int ik=0; ik<nks; ik++)
	{
		weight[ik] = weight_in[ik];
	}
	this->nbands = nbands_in;
	
	assert(nbands > 0);
	this->Mkb.create(nks, nbands);//mohan add 2010-05-02
	this->Mkb_used.create(nks, nbands);
	this->Mkb.zero_out();
	this->Mkb_used.zero_out();
	
	this->ne = ne_in;
	this->nwfc = nwfc_in;
	this->nwfc2 = nwfc2_in;
	this->nwfc_all = nwfc_all_in;
	
	
	if(test==1)
	{
	cout << "\n" << setw(10) << "nwfc"
		<< setw(10) << "nwfc2" 
		<< setw(10) << "nwfc_all";
	
	cout << "\n" << setw(10) << nwfc
		<< setw(10) << nwfc2
		<< setw(10) << nwfc_all;
	}

	this->allocate();

	return;
}

void Step_Data::allocate()
{
    this->Qin.create(nks, nbands, nwfc_all, ne);
    return;
}

// be called in Multi_Zeta.cpp
void Step_Data::initQS(void)
{
	if(test==1) TITLE("Step_Data","initQS");
	assert(nks > 0);
	assert(nbands > 0);
	assert(nwfc > 0); // nwfc is for <J|Bloch>, <J|J>
	assert(ne > 0);

	delete[] Qtrial;
	delete[] Strial;
    delete[] Qoverlap;
    delete[] Soverlap;
    delete[] Sinv;
    delete[] inverse_S;

    this->Qoverlap = new ComplexMatrix[nks];
    this->Soverlap = new ComplexMatrix[nks];
    this->Qtrial = new ComplexMatrix[nks];
    this->Strial = new ComplexMatrix[nks];
    this->Sinv = new ComplexMatrix[nks];
    this->inverse_S = new Inverse_Matrix_S[nks];

    for(int ik=0; ik<nks; ik++)
    {
        Qoverlap[ik].create(nwfc2, nbands);
        Soverlap[ik].create(nwfc2, nwfc2);
        Qtrial[ik].create(nwfc2, nbands);
        Strial[ik].create(nwfc2, nwfc2);
        Sinv[ik].create(nwfc2, nwfc2);
		inverse_S[ik].init(nwfc2);					// Peize Lin update 2015-12-05
    }

	return;
}



