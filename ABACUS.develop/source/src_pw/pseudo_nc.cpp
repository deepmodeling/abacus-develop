/* pseudo_nc.cpp */
#include "pseudo_nc.h"
#include "global.h"

pseudo_nc::pseudo_nc()
{
//	cout << "\n === pseudo_nc === ";
	lll = new int[1];
}

pseudo_nc::~pseudo_nc()
{
//	cout << "\n exit pseudo_nc ";
	delete[] lll;
}

//---------------------------------------------------------------------
void pseudo_nc::set_pseudo_nc(const Pseudopot_upf &upf)
{
	//-----------------------------------------------------------------

	int i, nb;

	this->set_pseudo_vl(upf);

	delete[] lll;
	lll = new int[nbeta];
	assert(lll != 0);

	for (i = 0;i < nbeta;i++)
	{
		lll[i] = upf.lll[i];
	}

	kkbeta = 0;

	for (nb = 0;nb < nbeta;nb++)
	{
		kkbeta = (upf.kkbeta[nb] > kkbeta) ? upf.kkbeta[nb] : kkbeta;
	}

	betar.create(upf.beta.nr, upf.beta.nc);
	
//	OUT("betar.nr",upf.beta.nr); // nbeta
//	OUT("betar.nc",upf.beta.nc); // mesh

	dion.create(nbeta, nbeta);

	betar = upf.beta;

	dion = upf.dion;

	nh = 0;

	for (nb = 0; nb < nbeta;nb++)
	{
		nh += 2 * lll [nb] + 1;
	} 

	return;
} // end subroutine set_pseudo_upf


void pseudo_nc::print_pseudo_nc(ofstream &ofs)
{
	print_pseudo_vl(ofs);
	ofs << "\n pseudo_nc : ";
	ofs << "\n kkbeta	" << kkbeta;
	ofs << "\n nh  " << nh;
	out.printr1_d(ofs, " lll : ", lll, nbeta);
	out.printrm(ofs, " betar : ", betar);
	out.printrm(ofs, " dion : ", dion);
	ofs << "\n ----------------------";
}

