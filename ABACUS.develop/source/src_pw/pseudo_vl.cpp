/* pseudo_vl.cpp */

#include "pseudo_vl.h"
#include "global.h"

pseudo_vl::pseudo_vl()
{
	vloc_at =  new double[1];
}

pseudo_vl::~pseudo_vl()
{
	delete[] vloc_at;
}

//---------------------------------------------------------------------
void pseudo_vl::set_pseudo_vl(const Pseudopot_upf &upf)
{
//	cout << "\n === set_pseudo_vl() === ";
	this->set_pseudo_at(upf);

	delete[] vloc_at;
	vloc_at = new double[mesh];
	assert(vloc_at != 0);

	for (int i = 0;i < mesh;i++)
	{
		vloc_at[i] = upf.vloc[i];
	}
//	cout << "\n End set_pseudo_vl() " << endl;
} 

void pseudo_vl::print_pseudo_vl(ofstream &ofs)
{
	ofs << "\n pseudo_vl:";
	print_pseudo_at(ofs);
	out.printr1_d(ofs, "vloc_at : ", vloc_at, mesh);
	ofs << "\n ----------------------------------- ";
}
