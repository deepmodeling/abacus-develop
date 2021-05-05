#include "../src_io/output.h"
#include "pseudo_us.h"

pseudo_us::pseudo_us()
{
//	cout << "\n === pseudo_us === ";
	rinner = new double[1];
}

pseudo_us::~pseudo_us()
{
//	cout << "\n exit pseudo_us ";
	delete[] rinner;
}

//---------------------------------------------------------------------
void pseudo_us::set_pseudo_us(const Pseudopot_upf &upf)
{
//	cout << "\n === set_pseudo_us() === ";
	//int i;

	this->set_pseudo_nc(upf);

/*	nqlc = upf.nqlc;
	nqf  = upf.nqf;

	delete[] rinner;
	rinner = new double[nqlc];

	for (i = 0;i < nqlc;i++)
	{
		this->rinner[i] = upf.rinner[i];
	}

	if (nbeta > 0)
	{
		this->qqq.create(nbeta, nbeta);
		this->qqq = upf.qqq;
	}

	if (mesh > 0 && nbeta > 0)
	{
		this->qfunc.create(nbeta, nbeta, mesh);
		this->qfunc = upf.qfunc;
	}

	this->qfcoef.create(nbeta, nbeta, upf.qfcoef.getBound3(), upf.qfcoef.getBound4() );
	this->qfcoef = upf.qfcoef;
*/

//	if (lspinorb && !upf.has_so)	// see USE spin_orb
//		cout << "\n upf_to_internal,At least one non s.o. pseudo,";	//-1);
//
//	lspinorb=lspinorb && upf.has_so;
//	cout << "\n End set_pseudo_us() " << endl;
} // end subroutine set_pseudo_upf

void pseudo_us::print_pseudo_us(ofstream &ofs)
{
	ofs << "\n\n internal format for pseudopotential ";
	output out;
	print_pseudo_nc(ofs);
	ofs << "\n pseudo_us : ";
	ofs << "\n nqf	" << nqf
	<< "\n nqlc	" << nqlc;
	out.printr1_d(ofs, " rinner : ", rinner, nbeta);
	out.printrm(ofs, " qqq : ", qqq);
	out.printr3_d(ofs, " qfunc : ", qfunc);
	out.printr4_d(ofs, " qfcoef : ", qfcoef);
	ofs << "\n ----------------------\n";
}


