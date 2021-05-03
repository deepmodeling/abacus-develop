#include "xc_type.h"
#include "src_global/global_function.h"
#include "src_pw/global.h"
#include "src_pw/exx_global.h"

xcfunc::xcfunc()
{
}

xcfunc::~xcfunc()
{
}

// mohan update 2009-12-15
const string exc[8] = { "NOX", "SLA", "SL1", "RXC", "OEP", "HF", "PB0X", "B3LP"};
const string corr[11] = { "NOC", "PZ", "VWN", "LYP", "PW", "WIG", "HL", "OBZ",
                          "OBW", "GL", "B3LP" };
const string gradx[10] = { "NOGX", "B88", "GGX", "PBX",  "RPB", "HCTH", "OPTX", "META", "PB0X", "B3LP"};
const string gradc[8] = { "NOGC", "P86", "GGC", "BLYP", "PBC", "HCTH", "META", "B3LP"};

// from function.f90
//-----------------------------------------------------------------------
void xcfunc::which_dft(const string *dft)
{
	//-----------------------------------------------------------------------
	// translates a string containing the exchange-correlation name
	// into internal indices iexch, icorr, igcx, igcc

	const int nxc = 8; // number of exchange functional
	const int ncc = 11; // number of correlation functional
	const int ngcx = 10; // number of gradient correction for exchange functional
	const int ngcc = 8; // number of gradient correction for correlation functional

	//int l=0;
	int i=0;
	int notset = -1;

	// (1) exchange
	this->iexch = notset;

	for (i = 0;i < nxc;i++) 
	{
		if (exc[i] == dft[0])
		{
			set_dft_value(iexch, i);
		}
	} 

	// (2) correlation
	this->icorr = notset;

	for (i = 0;i < ncc;i++)
	{
		if (corr[i] == dft[1])
		{	
			set_dft_value(icorr, i);
		}
	} 

	// (3) gradient correction, exchange
	this->igcx = notset;
	for (i = 0;i < ngcx;i++)
	{
		if (gradx[i] == dft[2])
		{
			set_dft_value(igcx, i);
		}
	}

	// (4) gradient correction, correlation
	this->igcc = notset;
	for (i = 0;i < ngcc;i++)
	{
		if (gradc[i] == dft[3])
		{
			set_dft_value(igcc, i);
		}
	}

	//======================= Second Part ===============================
	// special case : BLYP => B88 for gradient correction on exchange

	static int itype = 0;
	++itype;
	stringstream ss;
	ss << " ELEMENT " << itype << " FUNCTIONAL : "; 
//	cout << ss.str() << dft[0] << " " << dft[1] << " " << dft[2] << " " << dft[3] << endl;
	
	if( match_one( dft, "PBE0"))
	{
		set_dft_value(iexch,6);
		set_dft_value(icorr,4);
		set_dft_value(igcx,8);
		set_dft_value(igcc,4);
	}	
	if( match_one( dft, "LDA"))
	{
		set_dft_value(iexch,1);
		set_dft_value(icorr,1);
		set_dft_value(igcx,0);
		set_dft_value(igcc,0);
	}	
	else if ( match_one(dft, "PBE") )
	{
		// special case : PBE
		set_dft_value(iexch, 1);
		set_dft_value(icorr, 4);
		set_dft_value(igcx, 3);
		set_dft_value(igcc, 4);
	} 
	else if( match_one( dft, "revPBE" ) )
	{
		// special case : revPBE
		set_dft_value(iexch,1);
		set_dft_value(icorr,4);
		set_dft_value(igcx, 4);
		set_dft_value(igcc, 4);
	}
	else if ( match_one(dft, "PBEsol") )
	{
		set_dft_value(iexch, 1);
		set_dft_value(icorr, 4);
		set_dft_value(igcx, 10);
		set_dft_value(igcc, 8);
	}
	else if ( match_one(dft, "WC") )
	{
		set_dft_value(iexch, 1);
		set_dft_value(icorr, 4);
		set_dft_value(igcx, 11);
		set_dft_value(igcc, 4);
	}	
	else if ( match_one( dft, "BLYP") )
	{
		set_dft_value(iexch, 1);
		set_dft_value(icorr, 3);
		set_dft_value(igcx, 1);
		set_dft_value(igcc, 3);
	}
	else if ( match_one(dft, "BP") )
	{
		// special case : BP = B88 + P86
		set_dft_value(igcx, 1);
		set_dft_value(igcc, 1);
	} 
	else if ( match_one(dft, "PW91") )
	{
		// special case : PW91 = GGX + GGC
		set_dft_value(iexch, 1);
		set_dft_value(icorr, 4);
		set_dft_value(igcx, 2);
		set_dft_value(igcc, 2);
	} 
	else if ( match_one(dft, "HCTH") )
	{
		// special case : HCTH already contains LDA exchange and correlation
		set_dft_value(iexch, 0);
		set_dft_value(icorr, 0);
	}
	else if ( match_one(dft, "OPTX") )
	{
		// special case : OPTX already contains LDA exchange
		set_dft_value(iexch, 0);
	}
	else if (match_one(dft, "OLYP") )
	{
		// special case : OLYP = OPTX + LYP
		set_dft_value(iexch, 0);
		set_dft_value(icorr, 3);
		set_dft_value(igcx, 6);
		set_dft_value(igcc, 3);
	}

	if (igcx == 6)
	{
		cerr << "\n which_dft,OPTX untested please test,";	//-igcx);
	}

	if (iexch == notset)
	{
		// Default value: Slater exchange
		set_dft_value(iexch, 1);
	}
	if (icorr == notset)
	{
		// Default value: Perdew-Zunger correlation
		set_dft_value(icorr, 1);
	}
	if (igcx == notset)
	{
		// Default value: no gradient correction on exchange
		set_dft_value(igcx, 0);
	}
	if (igcc == notset)
	{
		// Default value: no gradient correction on correlation
		set_dft_value(igcc, 0);
	}

//	cout << "\n iexch = " << iexch; 
//	cout << "\n icorr = " << icorr; 
//	cout << "\n igcx = " << igcx; 
//	cout << "\n igcc = " << igcc << endl; 

	//cout << "\n corr = " << exc[icorr];
	//cout << "  corr = " << corr[icorr];
	//cout << "  gradx = " << gradx[igcx];
	//cout << "  gradc = " << gradc[igcc] << endl;
	//'-'//corr (icorr) //'-'//gradx (igcx) //'-'//gradc (igcc)
	//      WRITE( stdout,'(a)') dftout
	
	copy_to_now();

	hybrid_first();
	
	return;
} // end subroutine which_dft

//-----------------------------------------------------------------------
void xcfunc::set_dft_value(int &m,const int i)
{
	int notset = - 1;

	if (m != notset && m != i)
	{
		cerr << "\n set_dft_value, two conflicting matching values,";	// 1);
	}
	m = i;
	return;
} //  end subroutine set_dft_value

void xcfunc::printdft(ofstream &ofs)
{
	ofs << "\n iexch = " << iexch
	<< "  -> " << exc [iexch];
	ofs	<< "\n  icorr = " << icorr
	<< " -> " << corr[icorr];
	ofs	<< "\n  igcx = " << igcx
	<< "  -> " << gradx[igcx];
	ofs	<< "\n  igcc = " << igcc
	<< "  -> " << gradc[igcc];
}


void xcfunc::ostreamdft(ostream &ofs) // zws add 20150108
{
	if ( iexch == 1 && icorr == 1 && igcx == 0 && igcc == 0 )
	{ 
		ofs << "PZ-LDA";
	}
	else if (iexch == 1 && icorr == 4 && igcx == 3 && igcc == 4 )
	{
		ofs << "PBE";
	}
	else
	{
		ofs <<  exc [iexch] ;
		ofs	<<  " " << corr[icorr] ;
		ofs	<<  " " << gradx[igcx] ;
		ofs	<<  " " << gradc[igcc] ;
	}	
}


bool xcfunc::match_one(const string* dft, const string &name)const
{
	for(int i=0; i<4; i++)
	{
		if(dft[i]==name) return 1; // match one of the four string.
	}
	return 0; // no one match
}

// Peize Lin add 2016-12-03
void xcfunc::copy_to_now()
{
	iexch_now = iexch;
	icorr_now = icorr;
	igcx_now  = igcx ;
	igcc_now  = igcc ;
}

// Peize Lin add 2016-12-03
void xcfunc::hybrid_first()
{
	// may do something
	WARNING("functional","file "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__)+" may error ");
	if(Exx_Global::Hybrid_Type::HF==exx_global.info.hybrid_type)
	{
		iexch_now = 1;
		igcx_now = 3;
	}
	else if(Exx_Global::Hybrid_Type::PBE0==exx_global.info.hybrid_type)
	{
		iexch_now = 1;
		igcx_now = 3;
	}
	else if(Exx_Global::Hybrid_Type::HSE==exx_global.info.hybrid_type)
	{
		iexch_now = 1;
		igcx_now = 3;
	}
}
