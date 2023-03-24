#ifndef DOS_TEST_H
#define DOS_TEST_H

#include<iostream>
#include<fstream>
#include"module_base/constants.h"

class DosPrepare
{
public:
	DosPrepare(int is_in, std::string fa_in, std::string fa1_in,
		double de_ev_in, double emax_ev_in,
		double emin_ev_in, double bcoeff_in,
		int nks_in, int nkstot_in, int nbands_in):
		is(is_in),fa(fa_in),fa1(fa1_in),
		de_ev(de_ev_in),emax_ev(emax_ev_in),
		emin_ev(emin_ev_in),bcoeff(bcoeff_in),
		nks(nks_in),nkstot(nkstot_in),nbands(nbands_in){}
	int is;
	std::string fa;
	std::string fa1;
	double de_ev;
	double emax_ev;
	double emin_ev;
	double bcoeff;
	int nks;
	int nkstot;
	int nbands;
	std::vector<int> isk;
	std::vector<double> wk;
	ModuleBase::matrix ekb;
	ModuleBase::matrix wg;
	void set_isk()
	{
		this->isk.reserve(nks);
		for(int i=0;i<nks;i++)
		{
			isk[i] = 0; //spin-unpolarized case, only 1 spin
		}
	}
	void read_wk()
	{
		this->wk.reserve(nks);
		std::ifstream ifs;
		std::string tmpstring;
		int dummy;
		double kx,ky,kz;
		ifs.open("./support/kpoints");
		while(ifs.good())
		{
			getline(ifs,tmpstring);
			getline(ifs,tmpstring);
			for(int ik=0; ik<nks; ++ik)
			{
				ifs>>dummy >>kx >>ky >>kz >>this->wk[ik];
				ifs.ignore(150,'\n');
				ifs.rdstate();
			}
		}
		ifs.close();
		for(int ik=0; ik<nks; ++ik)
		{
			wk[ik] *= 2.0;
		}
	}
	void read_istate_info()
	{
		this->ekb.create(nks,nbands);
		this->wg.create(nks,nbands);
		std::ifstream ifs;
		std::string tmpstring;
		ifs.open("./support/istate.info");
		int dummy;
		while(ifs.good())
		{
			for(int ik=0; ik<nks; ++ik)
			{
				if(ik==0)
				{
					getline(ifs,tmpstring);
				}
				else
				{
					getline(ifs,tmpstring);
					getline(ifs,tmpstring);
					getline(ifs,tmpstring);
				}
				for(int ib=0; ib<nbands; ++ib)
				{
					ifs>> dummy >> this->ekb(ik,ib) >> this->wg(ik,ib);
					ifs.ignore(150,'\n');
				}
				ifs.rdstate();
			}
		}
		ifs.close();
		this->ekb *= 1.0/ModuleBase::Ry_to_eV; 
	}
};

#endif
