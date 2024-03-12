#include "Type_Information.h"

Type_Information::Type_Information()
{
	n = new int[1];
	nbase = new int[1];
	fa = new string[1];
	state = "new";	// mohan 2009-08-27
}

Type_Information::~Type_Information()
{
	delete[] n;
	delete[] nbase;
	delete[] fa;
}

void Type_Information::init() 
{ 
	assert(lmax>=-1); 
	delete[] n;
	delete[] nbase;
	delete[] fa;
	n = new int[lmax+1];
	nbase = new int[lmax+1];
	fa = new string[lmax+1];
	for(int i=0; i<lmax+1; i++)
	{
		n[i] = nbase[i] = 0;
	}
	return;
}

void Type_Information::cal_nmax()
{
	this->nmax = 0;
	for(int l=0; l<lmax+1; l++)
	{
		this->nmax = std::max( nmax, n[l]+nbase[l] );
	}	
//	cout << "\n nmax = " << nmax;
	return;
}

void Type_Information::cal_nw(void)
{
	TITLE("Type_Information","cal_nw");
	// where is this lmax from ?
	// in Multi_Zeta,
	// read in from file directly
	// for each level and each type.

	// be called in SpillageStep. 

	// nw2: number of total local orbitals.
	// (ia, l, n, m)
	// Including multi-zeta orbitals,
	// which are made up from same Jl(q).
    this->nw2 = 0;
    for(int l=0; l<lmax+1; l++)
    {
		if( fa[l] == "f" )
		{
        	nw2 += (2*l+1)*n[l]*na;
		}
		else if( fa[l] == "a" )
		{
			nw2 += 1*na; // 1 means average orbital used,
			// for example, d orbitals 5 ==> 1
			// f orbitals, 7 ==> 1
		}
        //cout << "\n nw = " << nw;
    }

	// nw: number of overlap we must save.
	// doesn't consider n[l],
	// the only condition is n[l] > 0 in this atom species.
	this->nw = 0;
	for(int l=0; l<lmax+1; l++)
	{
		// warning : in different step, n[l] value may be different.
		if(n[l]>0)
		{
			nw += (2*l+1)*na;
		}
	}

	cout << " cal_nw::nw = " << nw;
	cout << " cal_nw::nw2 = " << nw2 << endl;
    return;
}


