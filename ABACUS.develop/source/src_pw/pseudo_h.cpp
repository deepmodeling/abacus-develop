//==========================================================
// AUTHOR : Lixin He,mohan
// DATE : 2008-11-11
//==========================================================
#include "xc_type.h"
#include "global.h"
#include "pseudo_h.h"
#include "tools.h"

pseudo_h::pseudo_h()
{
	els = new string[1];
	lchi = new int[1];
	oc = new double[1];
	jjj = new double[1];
	jchi = new double[1];
	nn = new int[1];
	has_so = false;
	zv = 0;
}

pseudo_h::~pseudo_h()
{
	delete[] els;
	delete[] lchi;
	delete[] oc;
	delete[] jjj;
	delete[] jchi;
	delete[] nn;
}

//---------------------------------------------------------------------
void pseudo_h::set_pseudo_h(const Pseudopot_upf &upf)
{
	//	TITLE("pseudo_h::set_pseudo_h");
	//   set "is"-th pseudopotential using the Unified Pseudopotential Format

	int i=0;
	this->nv = upf.nv;// UPF file version number
	this->psd = upf.psd;
	this->pp_type = upf.pp_type;

	this->tvanp = upf.tvanp;// if USPP
	this->nlcc = upf.nlcc;// Non linear core corrections( bool ?)
	
	for(int i=0; i<4; i++)
	{
		this->dft[i] = upf.dft[i];
	}

	this->zv = upf.zp;
	this->etotps = upf.etotps;
	this->ecutwfc = upf.ecutwfc;
	this->ecutrho = upf.ecutrho;

	this->lmax = upf.lmax;
	this->mesh = upf.mesh;
	
	// mohan update 2021-02-22
	//  max number of points in the atomic radial mesh
	int ndmx = 2000; 
	if (this->mesh > ndmx)
	{
		cout << "\n set_pseudo_h, too many grid points,";
	}

	this->nchi = upf.nwfc;
	this->nbeta = upf.nbeta;

	delete[] els;
	this->els = new string[nchi];
	assert(this->els != 0);

	delete[] lchi;
	this->lchi = new int[this->nchi];
	assert(this->lchi != 0);

	delete[] oc;
	this->oc = new double[nchi];
	assert(this->oc != 0);

	for (i = 0;i < nchi;i++)
	{
		this->els[i] = upf.els[i];
		this->lchi[i] = upf.lchi[i];
		this->oc[i] = upf.oc[i];
	}

	delete[] jjj;
	this->jjj = new double[nbeta];
	assert(this->jjj != 0);

	delete[] nn;
	this->nn = new int[nchi];
	assert(this->nn != 0);

	delete[] jchi;
	this->jchi = new double[nchi];
	assert(this->jchi != 0);

	this->has_so = upf.has_so;//added by zhengdy-soc
	if (this->has_so)
	{ 
		for (i = 0;i < nchi;i++){
			this->nn[i] = upf.nn[i];
			this->jchi[i] = upf.jchi[i];
		}
		for (i = 0;i < upf.nbeta;i++)
		{
			this->jjj[i]  = upf.jjj [i];
		}
	}
	else
	{
		for (i = 0;i < nchi;i++){
			this->nn[i] = 0;
			this->jchi[i] = 0;
		}
		for (i = 0;i < upf.nbeta;i++)
		{
			this->jjj[i]  = 0;
		}
	} 

	/*
		jchi = new double[nwfc];
		if (upf.has_so){
			for(i=0;i<upf.nwfc;i++){
				jchi[i] = upf.jchi[i];
			}
		}else{
			for(i=0;i<upf.nwfc;i++){
				jchi[i] = 0;
			}
		} // endif
	*/

//	cout << "\n End set_pseudo_h() " << endl;
	return;
} // end subroutine set_pseudo_upf

void pseudo_h::print_pseudo_h(ofstream &ofs)
{
	ofs << "\n pseudo_info :";
	ofs	<< "\n nv		" << nv;
	ofs	<< "\n psd	" << psd;
	ofs	<< "\n pp_type	" << pp_type;
	ofs	<< "\n tvanp	" << tvanp;
	ofs	<< "\n nlcc	" << nlcc;
	ofs	<< "\n dft	" << dft;
	ofs	<< "\n zv		" << zv;
	ofs	<< "\n etotps	" << etotps;
	ofs	<< "\n ecutwfc " << ecutwfc;
	ofs	<< "\n ecutrho " << ecutrho;
	ofs	<< "\n lmax	" << lmax;
	ofs	<< "\n mesh	" << mesh;
	ofs	<< "\n nchi	" << nchi;
	ofs	<< "\n nbeta	" << nbeta;
//	out.printr1_d(ofs," els: ", els, nchi);
	out.printr1_d(ofs, " lchi: ", lchi, nchi);
	out.printr1_d(ofs, " oc: ", oc, nchi);
	ofs << "\n ----------------------";
}

