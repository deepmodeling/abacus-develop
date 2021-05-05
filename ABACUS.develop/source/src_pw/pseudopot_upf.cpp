//==========================================================
// Author: Lixin He,mohan
// DATE : 2009-02-26
//==========================================================
#include "../src_io/output.h"
#include "pseudopot_upf.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring>		// Peize Lin fix bug about strcpy 2016-08-02

int Number[2];

using namespace std;

Pseudopot_upf::Pseudopot_upf()
{
	this->els = new string[1];
	this->lchi = new int[1];
	this->oc = new double[1];

	this->r = new double[1];
	this->rab = new double[1];
	this->vloc = new double[1];

	this->kkbeta = new int[1];
	this->lll = new int[1];

	this->rho_at = new double[1];
	this->rho_atc = new double[1];

	this->nn = new int[1];//zhengdy-soc
	this->jchi = new double[1];
	this->jjj = new double[1];

        functional_error = 0;//xiaohui add 2015-03-24
}

Pseudopot_upf::~Pseudopot_upf()
{
	delete [] els;  //header_15
	delete [] lchi; //header_16
	delete [] oc;   //header_17

	delete [] r;    //mesh_1
	delete [] rab;  //mesh_2
	delete [] vloc;  //local_1

	delete [] kkbeta; // nl_1
	delete [] lll; // nl_2

	delete [] rho_at;// psrhoatom_1
	delete [] rho_atc;

	delete [] nn;
	delete [] jjj;
	delete [] jchi;
}

int Pseudopot_upf::init_pseudo_reader(const string &fn)
{
    if(test_pp) TITLE("Pseudopot_upf","init");
    // First check if this pseudo-potential has spin-orbit information
    ifstream ifs(fn.c_str(), ios::in);

    if (!ifs)
    {
        // can't find the file.
        return 1;
    }

    //cout << "global_pseudo_type =" << global_pseudo_type << endl;
    if(global_pseudo_type=="auto") //{zws
	{
		set_pseudo_type(fn);
	} //}

	// read in the .UPF type of pseudopotentials
	if(global_pseudo_type=="upf")
	{
		int info = read_pseudo_upf(ifs);
		//this->print_pseudo_upf(ofs_running);
		return info;
	}
	// read in the .vwr type of pseudopotentials
	else if(global_pseudo_type=="vwr")
	{
		int info = read_pseudo_vwr(ifs);
	//	this->print_pseudo_upf(ofs_running);
		return info;
	}
	else if(global_pseudo_type=="upf201")
	{
		int info = read_pseudo_upf201(ifs);
		//this->print_pseudo_upf(ofs_running);
		return info;
	}




	return 0;
}


//----------------------------------------------------------
// setting the type of the pseudopotential file
//----------------------------------------------------------
int Pseudopot_upf::set_pseudo_type(const string &fn) //{zws add
{
    ifstream pptype_ifs(fn.c_str(), ios::in);
    string dummy, strversion;
	//int ierr = 0;

	if (pptype_ifs.good())
	{
		getline(pptype_ifs,dummy);
		//cout << "dummy : " << dummy << endl;
		//if ( dummy == "<PP_INFO>" )
		//{
		//	global_pseudo_type = "upf";
		//}

		stringstream wdsstream(dummy);
		getline(wdsstream,strversion,'"');
		getline(wdsstream,strversion,'"');
		//cout << "strversi:" << strversion << endl;

		if ( trim(strversion) == "2.0.1" )
		{
			global_pseudo_type = "upf201";
		}
		else
		{
			global_pseudo_type = "upf";
		}
	}
	return 0;
}

string& Pseudopot_upf::trim(string &in_str)
{
    static const string deltri = " \t" ; // delete tab or space
    string::size_type position = in_str.find_first_of(deltri, position);
    if (position == string::npos)
        return in_str;
    return trim(in_str.erase(position, 1) );
}

string Pseudopot_upf::trimend(string &in_str)
{
    const string &deltri =" \t" ;
    string::size_type position = in_str.find_last_not_of(deltri)+1;
    //cout << "--" << position << "--" << endl;
    string tmpstr=in_str.erase(position);
    //cout << "--" << tmpstr << "--" << endl;
    return tmpstr.erase(0,tmpstr.find_first_not_of(deltri));
} //}zws



//----------------------------------------------------------
// This code is used to read in vwr pseudopotential format,
// which is used in PEtot mostly. Now I only use LDA, so if 
// PBE or other functionals are used, one needs to change 
// the following code. The vwr format needs we to generate 
// NL projectors by ourself. One way to check this is correct
// is to use opium to generate .ncpp pseudopotential first, 
// which contains the same informaiton in vwr, then we had
// both UPF format from upftools in Quantum Espresso and
// we can write a short code to transform ncpp to vwr.
// Then compare the two results.
// mohan 2013-05-25
//-----------------------------------------------------------
int Pseudopot_upf::read_pseudo_vwr(ifstream &ifs)
{
	ofs_running << " -------------------------------------------------" << endl;
	cout << " READ IN VWR TYPE PSEUDOPOTENTIALS." << endl;
	ofs_running << " Read in vwr type pseudopotentials " << endl;


	// --------------------------------------
	// (1) read in data 
	// --------------------------------------
	this->dft[0]="SLA";
	this->dft[1]="PZ";
	this->dft[2]="NOGX";
	this->dft[3]="NOGC";
	this->pp_type="NC";
	this->tvanp=false;
	ofs_running << " Always use PZ-LDA by now." << endl;
	

	// (1) read in mesh
	string value;
	int length=0;
	ifs >> value; length = value.find(","); value.erase(length,1);
	mesh = std::atoi( value.c_str() );
	//the mesh should be odd, which is forced in Simpson integration
	if(mesh%2==0) 
	{
		mesh=mesh-1;
		ofs_running << " Mesh number - 1, because we need odd number, \n this may affect some polar atomic orbitals." << endl;
	}
	ofs_running << setw(15) << "MESH" << setw(15) << mesh << endl;
	// (2) read in nlcc: nonlinear core correction
	ifs >> value; length = value.find(","); value.erase(length,1);
	nlcc = std::atoi( value.c_str() );
	ofs_running << setw(15) << "NLCC" << setw(15) << nlcc << endl;
	// (3) iatom : index for atom 
	ifs >> value; length = value.find(","); value.erase(length,1);
	psd = value;
	ofs_running << setw(15) << "ATOM" << setw(15) << psd << endl;
	// (4) valence electron number
	ifs >> value; length = value.find(","); value.erase(length,1);
	zp = std::atoi( value.c_str() );
	ofs_running << setw(15) << "Z(VALENCE)" << setw(15) << zp << endl;
	// (5) spd_loc, which local pseudopotential should I choose
	ifs >> value; length = value.find(","); value.erase(length,1);
	spd_loc = std::atoi( value.c_str() );
	ofs_running << setw(15) << "LOC(spd)" << setw(15) << spd_loc << endl;
	// (6) read in the occupations
	double* tmp_oc = new double[3];
	ifs >> value; length = value.find(","); value.erase(length,1);
	tmp_oc[0]= std::atoi( value.c_str() );
	ifs >> value; length = value.find(","); value.erase(length,1);
	tmp_oc[1]= std::atoi( value.c_str() );
	ifs >> value; length = value.find(","); value.erase(length,1);
	tmp_oc[2]= std::atoi( value.c_str() );
	ofs_running << setw(15) << "OCCUPATION" << setw(15) << tmp_oc[0] 
	<< setw(15) << tmp_oc[1] << setw(15) << tmp_oc[2] << endl;
	// (7) spin orbital
	ifs >> has_so;


	// label to count the projector or atomic wave functions	
	getline(ifs,value);
	int iref_s, iref_p, iref_d;
	ifs >> iref_s >> iref_p >> iref_d;
	ofs_running << setw(15) << "Vnl_USED" << setw(15) << iref_s 
	<< setw(15) << iref_p << setw(15) << iref_d << endl;
	if(spd_loc==1) iref_s=0;
	else if(spd_loc==2) iref_p=0;
	else if(spd_loc==3) iref_d=0;
	ifs >> iTB_s >> iTB_p >> iTB_d;
	ofs_running << setw(15) << "Orb_USED" << setw(15) << iTB_s 
	<< setw(15) << iTB_p << setw(15) << iTB_d << endl;
	
	
	// calculate the number of wave functions
	this->nwfc = 0;
	if(iTB_s) ++nwfc;
	if(iTB_p) ++nwfc;
	if(iTB_d) ++nwfc;
	ofs_running << setw(15) << "NWFC" << setw(15) << nwfc << endl;
	// allocate occupation number array for wave functions
	delete[] this->oc;
	delete[] els;
	this->oc = new double[nwfc];
	els = new string[nwfc];
	// set the value of occupations
	delete[] lchi;
	lchi = new int[nwfc];
	int iwfc=0;
	if(iTB_s){oc[iwfc]=tmp_oc[0];lchi[iwfc]=0;els[iwfc]="S";++iwfc;}
	if(iTB_p){oc[iwfc]=tmp_oc[1];lchi[iwfc]=1;els[iwfc]="P";++iwfc;}
	if(iTB_d){oc[iwfc]=tmp_oc[2];lchi[iwfc]=2;els[iwfc]="D";++iwfc;}
	delete[] tmp_oc;
	getline(ifs,value);


	// global variables that will be used
	// in other classes.
	delete[] r;
	delete[] rab;
	delete[] vloc;
	this->r = new double[mesh];
	this->rab = new double[mesh];
	this->vloc = new double[mesh];
	ZEROS(r,mesh);
	ZEROS(rab,mesh);
	ZEROS(vloc,mesh);
	delete[] rho_at;
	delete[] rho_atc;
	rho_at = new double[mesh];
	rho_atc = new double[mesh];
	ZEROS(rho_at, mesh);
	ZEROS(rho_atc, mesh);
	// local variables in this function
	this->vs = new double[mesh];
	this->vp = new double[mesh];
	this->vd = new double[mesh];
	this->ws = new double[mesh];
	this->wp = new double[mesh];
	this->wd = new double[mesh];
	ZEROS(vs,mesh);
	ZEROS(vp,mesh);
	ZEROS(vd,mesh);
	ZEROS(ws,mesh);
	ZEROS(wp,mesh);
	ZEROS(wd,mesh);
	string line;
	if(spd_loc>0 && nlcc==0)
	{
		for(int ir=0; ir<mesh; ++ir)
		{
			// it's an interesting question whether
			// ws[ir] has 1/sqrt(4pi)
			ifs >> r[ir] >> vs[ir] >> vp[ir] >> vd[ir] 
				>> ws[ir] >> wp[ir] >> wd[ir];
			getline(ifs, line);
		}
	}
	else if(spd_loc==0 && nlcc==0)
	{
		for(int ir=0; ir<mesh; ++ir)
		{
			ifs >> r[ir] >> vs[ir] >> vp[ir] >> vd[ir] 
				>> ws[ir] >> wp[ir] >> wd[ir] >> vloc[ir];
			getline(ifs, line);
		}
	}
	else if(spd_loc>0 && nlcc==1)
	{
		for(int ir=0; ir<mesh; ++ir)
		{
			ifs >> r[ir] >> vs[ir] >> vp[ir] >> vd[ir] 
				>> ws[ir] >> wp[ir] >> wd[ir] >> rho_atc[ir];
			getline(ifs, line);
		}
	}
	else if(spd_loc==0 && nlcc==1)
	{
		for(int ir=0; ir<mesh; ++ir)
		{
			ifs >> r[ir] >> vs[ir] >> vp[ir] >> vd[ir] 
				>> ws[ir] >> wp[ir] >> wd[ir] >> vloc[ir] >> rho_atc[ir];
			getline(ifs, line);
		}
	}
	// Hartree to Rydberg
	for(int ir=0; ir<mesh; ++ir)
	{
		vs[ir] *= 2.0;
		vp[ir] *= 2.0;
		vd[ir] *= 2.0;
		vloc[ir] *= 2.0;
	}


	// --------------------------------------
	// (2) check unit
	// --------------------------------------
	// calculate rab;
	// rab may not be accurate enough
	rab[0] = r[0];
	for(int ir=1; ir<mesh-1; ++ir)
	{
		rab[ir]=(r[ir+1]-r[ir-1])/2.0;
	}
	// check unit of vs, vp, vd
	double units = 0.0;
	double unitp = 0.0;
	double unitd = 0.0;
	for(int ir=1; ir<mesh-1; ++ir)
	{
		double dr = (r[ir+1]-r[ir-1])/2.0;
		units += ws[ir] * ws[ir] * r[ir] * r[ir] * dr;
		unitp += wp[ir] * wp[ir] * r[ir] * r[ir] * dr;
		unitd += wd[ir] * wd[ir] * r[ir] * r[ir] * dr;
	}
	ofs_running << setw(15) << "WFC_UNIT" << setw(15) << units 
	<< setw(15) << unitp << setw(15) << unitd << endl; 


	// because only the rank=0 procesor read the pseudopotential
	// information, in order to make all the processors to stop
	// the job, we need to return the error information first.
	// we need to choose a threshold for the deviation of the
	// norm of pseudo atomic orbitals, I set 0.2
	// mohan 2013-06-28
	if( abs(units-1.0) > 0.2 && (iTB_s==1 || iref_s==1)) {return 3;}
	if( abs(unitp-1.0) > 0.2 && (iTB_p==1 || iref_p==1)) {return 3;}
	if( abs(unitd-1.0) > 0.2 && (iTB_d==1 || iref_d==1)) {return 3;}


	// calculate the phi*r*sqrt(4pi)
	this->chi.create(nwfc,mesh);
	for(int ir=0; ir<mesh; ++ir)
	{
		int iwfc=0;
		if(iTB_s==1){chi(iwfc,ir) = ws[ir]*r[ir];++iwfc;}
		if(iTB_p==1){chi(iwfc,ir) = wp[ir]*r[ir];++iwfc;}
		if(iTB_d==1){chi(iwfc,ir) = wd[ir]*r[ir];++iwfc;}
	}
	// rho atom 
	for(int ir=0; ir<mesh; ++ir)
	{
		for(int iwfc=0; iwfc<nwfc; ++iwfc)
		{
			rho_at[ir] += oc[iwfc]*chi(iwfc,ir)*chi(iwfc,ir);
		}
	}



	// --------------------------------------
	// (4) local pseudopotential 
	// --------------------------------------
	if(spd_loc==0) for(int ir=0; ir<mesh; ++ir)
	{
			// do nothing	
	}
	else if(spd_loc==1) for(int ir=0; ir<mesh; ++ir) vloc[ir] = vs[ir]; 
	else if(spd_loc==2) for(int ir=0; ir<mesh; ++ir) vloc[ir] = vp[ir];
	else if(spd_loc==3) for(int ir=0; ir<mesh; ++ir) vloc[ir] = vd[ir];
	else if(spd_loc==12) for(int ir=0; ir<mesh; ++ir) vloc[ir] = (vs[ir]+vp[ir])/2.0;
	else if(spd_loc==13) for(int ir=0; ir<mesh; ++ir) vloc[ir] = (vs[ir]+vd[ir])/2.0;
	else if(spd_loc==23) for(int ir=0; ir<mesh; ++ir) vloc[ir] = (vp[ir]+vd[ir])/2.0;


	// --------------------------------------
	// (5) setup nonlocal pseudopotentials
	// --------------------------------------
	// for non-local pseudopotentials.
	if(iref_d==1) lmax=2;
	else if(iref_p==1) lmax=1;
	else if(iref_s==1) lmax=0;
	else
	{
		cout << "\n !!! READ THIS FIRST !!!" << endl;
		cout << " Could not decide which is the max angular momentum from .vwr pseudopotential file." << endl;
		cout << " No reference states in .vwr pseudopotential file." << endl;
		cout << " That's incorrect, please check the refenrece states in .vwr file.";
		cout << "\n !!! READ THIS FIRST !!!" << endl;
		return 3;	
	}
	// no projectors now
	this->nbeta = 0;
	if(iref_s==1) ++nbeta; // add one s projector
	if(iref_p==1) ++nbeta; // add one p projector
	if(iref_d==1) ++nbeta; // add one p projector
	ofs_running << setw(15) << "NPROJ" << setw(15) << nbeta << endl;
	this->nd = this->nbeta;
	ofs_running << setw(15) << "N-Dij" << setw(15) << nd << endl;
	// calculate the angular momentum for each beta
	delete[] lll;
	lll = new int[nbeta]; 
	int icount=0;
	if(iref_s==1) {lll[icount]=0; ++icount;}// s projector
	if(iref_p==1) {lll[icount]=1; ++icount;}// p projector
	if(iref_d==1) {lll[icount]=2; ++icount;}// p projector
	for(int i=0; i<nbeta; ++i)
	{
		ofs_running << " lll[" << i << "]=" << lll[i] << endl;
	}
	// kkbeta(nbeta): number of mesh points for projector i (must be .le.mesh )
	delete[] kkbeta;
	kkbeta = new int[nbeta];
	for(int ib=0; ib<nbeta; ++ib)
	{
		kkbeta[ib] = mesh;
	}
	// nonlocal projector
	beta.create(nbeta,mesh);
	// coefficients
	dion.create(nbeta,nbeta);


	// --------------------------------------
	// (6) generate nonlocal pseudopotentials
	// --------------------------------------
	// tmp function to evaluate < beta | delta_v | beta>
	double* func = new double[mesh];
	// tmp value (vs, vp or vd)
	double* vl = new double[mesh];
	// tmp wave function (ws, wp or wd with r)
	double* wlr = new double[mesh];
	ZEROS(func, mesh);
	ZEROS(vl, mesh);
	ZEROS(wlr, mesh);
	double rcut = 5.0/1.03;
	ofs_running << setw(15) << "RCUT_NL" << setw(15) << rcut << endl;
	for(int ib=0; ib<nbeta; ++ib)
	{
		double coef = 0.0;
		const int lnow = lll[ib];
		if(lnow==0) for(int ir=0; ir<mesh; ++ir){vl[ir]=vs[ir]; wlr[ir]=ws[ir]*r[ir];}
		else if(lnow==1) for(int ir=0; ir<mesh; ++ir){vl[ir]=vp[ir]; wlr[ir]=wp[ir]*r[ir];}
		else if(lnow==2) for(int ir=0; ir<mesh; ++ir){vl[ir]=vd[ir]; wlr[ir]=wd[ir]*r[ir];}
		// for non-local projectors
		// note that < phi | dV | phi > integration must have 4pi,
		// this 4pi is also needed in < phi | phi > = 1 integration.
		// However, this phi has sqrt(sphi) already because I
		// found < phi | phi > = 1 directly.
		ofs_running << " Projector index = " << ib+1 << ", L = " << lnow << endl;
		for(int ir=2; ir<mesh-1; ++ir)
		{
			// p nl
			beta(ib,ir)=(vl[ir]-vloc[ir])*wlr[ir];
			if(r[ir]<rcut)
			{
				coef=coef+(vl[ir]-vloc[ir])*wlr[ir]*wlr[ir]*(r[ir+1]-r[ir-1])/2.0;
			}
	//		ofs_running << setw(15) << r[ir] << setw(25) << beta(ib,ir) << endl;
		}


// In quantum espresso they did this:
//		dion(ib,ib)=1.0/coef;

        if(coef<0.0) dion(ib,ib) = -1.0;
        if(coef>=0.0) dion(ib,ib) = 1.0;
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// suppose wave function have sqrt(4pi) already
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        coef=1.0/sqrt(abs(coef));
		ofs_running << setw(25) << "1/sqrt(<phi|deltaV|phi>)" << setw(15) << coef << endl;
		for(int ir=0; ir<mesh; ++ir)
		{
			beta(ib,ir) *= coef;
			// --------- FOR TEST ---------
			if(ib>2)
			{
//				beta(ib,ir) *= 0.0; // for test, disable Non-local
			}
			// --------- FOR TEST ---------
		}
		
	}

	// print out the projector.
	/*
	ofs_running << " Nonlocal projector : " << endl;
	for(int ir=0; ir<mesh; ++ir)
	{
		if(r[ir]<rcut)
		{
			ofs_running << " " << setw(15) << r[ir]; 
			for(int ib=0; ib<nbeta; ++ib)
			{
				ofs_running << setw(25) << beta(ib,ir);
			}
			ofs_running << endl;
		}
	}
	*/
	ofs_running << " -------------------------------------------------" << endl;


	// --------------------------------------
	// (5) clean up 
	// --------------------------------------
	delete[] func;
	delete[] vl;
	delete[] wlr;

	delete[] vs;
	delete[] vp;
	delete[] vd;

	delete[] ws;
	delete[] wp;
	delete[] wd;

	return 0;
}

//  read pseudopot_upf potential "upf" in the Unified
//  Pseudopot_upfpotential Format
int Pseudopot_upf::read_pseudo_upf(ifstream &ifs)
{
	string dummy;
	int ierr = 0;
	this->has_so = false;

	//addinfo_loop
	ifs.rdstate();

	while (ifs.good())
	{
		ifs >> dummy;
		if(dummy=="<PP_ADDINFO>")
		{
			ierr = 1;
			break;
		}
		ifs.rdstate();
	}

	if (ierr == 1) 
	{
		has_so = true;
	}
	
	// Search for Header
	// This version doesn't use the new routine SCAN_BEGIN
	// because this search must set extra flags for
	// compatibility with other pp format reading

	ierr = 0;

	ifs.clear();
	ifs.seekg(0);
	ifs.rdstate();

	// header_loop:
	while (ifs.good())
	{
		ifs >> dummy;
		if(dummy=="<PP_HEADER>")
		{
			ierr = 1;
			read_pseudo_header(ifs);
			SCAN_END(ifs, "</PP_HEADER>");
			break;
		}
	}

	if (ierr == 0) 
	{
		// 2: something in pseudopotential file not match.
		return 2;
	}

	// Search for mesh information
	if ( SCAN_BEGIN(ifs, "<PP_MESH>") )
	{
		read_pseudo_mesh(ifs);
		SCAN_END(ifs, "</PP_MESH>");
	}

	// If  present, search for nlcc
	if (this->nlcc)
	{
		SCAN_BEGIN(ifs, "<PP_NLCC>"); 
		read_pseudo_nlcc(ifs);
		SCAN_END(ifs, "</PP_NLCC>");
	}

	// Search for Local potential
	SCAN_BEGIN(ifs, "<PP_LOCAL>");
	read_pseudo_local(ifs);
	SCAN_END(ifs, "</PP_LOCAL>");

	// Search for Nonlocal potential
	SCAN_BEGIN(ifs, "<PP_NONLOCAL>");
	read_pseudo_nl(ifs);
	SCAN_END(ifs, "</PP_NONLOCAL>");

	// Search for atomic wavefunctions
	SCAN_BEGIN(ifs, "<PP_PSWFC>");
	read_pseudo_pswfc(ifs);
	SCAN_END(ifs, "</PP_PSWFC>");

	// Search for atomic charge
	SCAN_BEGIN(ifs, "<PP_RHOATOM>");
	read_pseudo_rhoatom(ifs);
	SCAN_END(ifs, "</PP_RHOATOM>");

	// Search for add_info
	if (has_so)
	{
//		WARNING_QUIT("read_pseudo_upf","Can't read UPF containing spin-orbital term.");
		SCAN_BEGIN (ifs,"<PP_ADDINFO>");//added by zhengdy-soc
		read_pseudo_so (ifs);
		SCAN_END (ifs,"</PP_ADDINFO>");
	}

	ifs.clear();
	ifs.seekg(0);

	// return 0: read in sucessfully.
	return 0;
}// end subroutine read_pseudopot_upf

void Pseudopot_upf::read_pseudo_header(ifstream &ifs)
{
	READ_VALUE(ifs, this->nv);// Version number
	READ_VALUE(ifs, this->psd);// Element label

	// Type of pseudo : NC or US
	READ_VALUE(ifs, this->pp_type);
	if(pp_type=="US")
	{
		this->tvanp = true;
	}
	else if(pp_type=="NC")
	{
		this->tvanp = false;
	}
	else
	{
		// A bug here!!! can't quit together.
		cout << " pp_type=" << pp_type << endl;
		WARNING_QUIT("Pseudopot_upf::read_pseudo_header","unknown pseudo type");
	}

	// If use nlcc
	string nlc;
	READ_VALUE(ifs, nlc);

	if (nlc == "T")
	{
		this->nlcc = true;
	}
	else
	{
		this->nlcc = false;
	}

	// mohan modify 2009-12-15
	ifs >> dft[0] >> dft[1] >> dft[2];
	READ_VALUE(ifs, this->dft[3]);
	
	// dft functional enforced to modify
	// mohan add 2010-07-15
	if(DFT_FUNCTIONAL!="none")
	{
		/*xiaohui modify 2015-03-24
		dft[0] = DFT_FUNCTIONAL;
		dft[1] = DFT_FUNCTIONAL;
		dft[2] = DFT_FUNCTIONAL;
		dft[3] = DFT_FUNCTIONAL;
		xiaohui modify 2015-03-24*/

		//xiaohui add 2015-03-23
		string dft_functional;
		if(dft[1] == "PZ")
		{
			dft_functional = "lda";
		}
		else if(dft[1] == "PW")
		{
			dft_functional = "pbe";
		}

		if(dft_functional != DFT_FUNCTIONAL)
		{
			functional_error = 1;

			cout << " dft_functional readin is: " << DFT_FUNCTIONAL << endl;
			cout << " dft_functional in pseudopot file is: " << dft_functional << endl;
			ofs_warning << " dft_functional readin is: " << DFT_FUNCTIONAL << endl;
			ofs_warning << " dft_functional in pseudopot file is: " << dft_functional << endl;
			//WARNING_QUIT("Pseudopot_upf::read_pseudo_header","input xc functional does not match that in pseudopot file");
		}
	}
	
	READ_VALUE(ifs, this->zp);
	READ_VALUE(ifs, this->etotps);

	ifs >> this->ecutwfc >> this->ecutrho;
	ifs.ignore(75, '\n');

	READ_VALUE(ifs, this->lmax);
	READ_VALUE(ifs, this->mesh);

	ifs >> this->nwfc >> this->nbeta ;
	ifs.ignore(75, '\n');
	ifs.ignore(75, '\n');

	delete[] els;
	delete[] lchi;
	delete[] oc;
	this->els = new string[nwfc];
	this->lchi = new int[nwfc];
	this->oc = new double[nwfc];

	ZEROS(lchi, nwfc); // angular momentum of each orbital
	ZEROS(oc, nwfc);//occupation of each orbital

	for(int i=0;i<nwfc;i++)
	{
		ifs >> els[i] >> this->lchi[i] >> this->oc[i];
	}
	return;
}

void Pseudopot_upf::read_pseudo_mesh(ifstream &ifs)
{
	assert(mesh>0);

	delete[] r;
	delete[] rab;
	this->r = new double[mesh];
	this->rab = new double[mesh];
	ZEROS(r,mesh);
	ZEROS(rab,mesh);

	int ir = 0;

	if( SCAN_BEGIN(ifs, "<PP_R>", false) )
	{
		for (ir = 0;ir < mesh;ir++)
		{
			ifs >> this->r[ir];
		}
		SCAN_END(ifs, "</PP_R>");
	}

	if( SCAN_BEGIN(ifs, "<PP_RAB>", false) )
	{
		for (ir = 0;ir < mesh;ir++)
		{
			ifs >> this->rab[ir];
		}
		SCAN_END(ifs, "</PP_RAB>");
	}
	return;
}

void Pseudopot_upf::read_pseudo_nlcc(ifstream &ifs)
{
	assert(mesh>0);
	delete[] rho_atc;
	this->rho_atc = new double[mesh];
	ZEROS(rho_atc, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rho_atc[ir];
	}
	return;
}

void Pseudopot_upf::read_pseudo_local(ifstream &ifs)
{
	assert(mesh>0);
	delete[] vloc;
	this->vloc = new double[mesh];
	ZEROS(vloc, mesh);

	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->vloc[ir];
	}

	
	// for test
/*
	ofstream ofs("vloc_UPF.dat");
	for(int ir=0; ir<mesh; ++ir)
	{
		ofs << r[ir] << " " << vloc[ir] << endl;
	}
	ofs.close();
*/

	return;
}

void Pseudopot_upf::read_pseudo_nl(ifstream &ifs)
{
//	int nb, mb, n, ir, idum, ldum, lp, i, ikk;
	int nb, mb, ir, idum;

	if (nbeta == 0)
	{
//		this->nqf = 0;
//		this->nqlc = 0;
		delete[] kkbeta;
		delete[] lll;
		this->kkbeta = new int[1];
		this->lll = new int[1];
		this->beta.create(1, 1);
		this->dion.create(1, 1);
//		this->rinner = new double[lmax];
//		this->qqq.create(1, 1);
//		this->qfunc.create(1, 1, mesh);
//		this->qfcoef.create(1, 1, 1, lmax);
		return;
	}
	else
	{
		delete[] kkbeta;
		delete[] lll;
		this->kkbeta = new int[nbeta]; 
		this->lll = new int[nbeta]; 
		this->beta.create(nbeta , mesh);
		this->dion.create(nbeta , nbeta);

		for(int i=0;i<nbeta;i++)
		{
			SCAN_BEGIN(ifs, "<PP_BETA>", false);
			ifs >> idum;
			READ_VALUE(ifs, this->lll[i]);// nl_1
			READ_VALUE(ifs, this->kkbeta[i]);// nl_2
			// number of mesh points for projectors

			for (ir=0;ir<kkbeta[i];ir++)
			{
				ifs >> this->beta(i, ir);// nl_3

				// --------- FOR TEST ---------
				// mohan test bug
				if(i>2)
				{
					beta(i,ir) = 0.0;
				}
				// --------- FOR TEST ---------
			}
			SCAN_END(ifs, "</PP_BETA>");
		}

		// DIJ
		SCAN_BEGIN(ifs, "<PP_DIJ>", false);
		READ_VALUE(ifs, this->nd);//nl_4
		for (int i=0; i<this->nd; i++)
		{
			double swap;
			ifs >> nb >> mb >> swap;
			nb--;
			mb--;
			assert( nb >= 0); // mohan add 2011-03-10
			assert( mb >= 0);
			this->dion(mb, nb) = swap;// nl_5
			this->dion(nb, mb) = swap;
		}
		SCAN_END(ifs, "</PP_DIJ>");

		// QIJ
		if (tvanp)
		{
			WARNING_QUIT("read_pseudo_nl","Ultra Soft Pseudopotential not available yet.");
			/*
			SCAN_BEGIN(ifs, "<PP_QIJ>", false);
			ifs >> this->nqf;//nl_6
			int xx = (this->nqf > 1) ? this->nqf : 1;
			ifs.ignore(75, '\n');
			this->nqlc = 2 * lmax  + 1;
			this->rinner = new double[nqlc]();

			this->qqq.create(nbeta , nbeta);
			this->qqq.zero_out();

			this->qfunc.create(nbeta , nbeta , mesh);
			this->qfcoef.create(nbeta , nbeta , nqlc , xx);

			if (nqf != 0)
			{
				SCAN_BEGIN(ifs, "<PP_RINNER>", false);

				for (i = 0;i < nqlc;i++)
				{
					ifs >> idum >> rinner[i];
				}

				SCAN_END(ifs, "</PP_RINNER>");
			}

			for (int  nb = 0;nb < nbeta;nb++)
			{
				for (int  mb = nb;mb < nbeta;mb++)
				{
					ifs >> idum >> idum >> ldum;
					ifs.ignore(75, '\n');

					if (ldum != lll[mb])
					{
						cerr << " read_pseudo_nl inconsistent angular momentum for Q_ij \n ";
					}

					ifs >> x;

					this->qqq(nb, mb) = x;
					ifs.ignore(75, '\n');

					// "Q_int"
					// qqq.assign(mb , nb , x);
					this->qqq(mb, nb) = x;

					for (n = 0;n < mesh;n++)
					{
						ifs >> x;
						this->qfunc(mb , nb , n) = x ;
					}

					for (n = 0;n < mesh;n++)
					{
						x = this->qfunc(mb,nb,n); // bad idea!
						this->qfunc(nb,mb,n) =  x;
					}

					// QFCOEF
					if (nqf > 0)
					{
						SCAN_BEGIN(ifs, "<PP_QFCOEF>", false);
						for (lp = 0;lp < nqlc;lp++)
						{
							for (i = 1;i < nqf;i++)
							{
								ifs >> x;
								this->qfcoef(nb , mb , lp , i) = x;
							}
						}
						SCAN_END(ifs, "</PP_QFCOEF>");
					}
				}
			}
			SCAN_END(ifs, "</PP_QIJ>");
			*/
		}
		else // not tvanp
		{
//			this->nqf  = 1;//nl_6
//			this->nqlc = 2 * lmax  + 1;//nl_7
//			this->rinner = new double[nqlc]();//nl_8
//			this->qqq.create(nbeta , nbeta);//nl_9
//			this->qqq.zero_out();
//			this->qfunc.create(nbeta , nbeta , mesh);//nl_10
//			this->qfcoef.create(nbeta , nbeta , nqlc , nqf);//nl_11
		}
	}
	return;
}

void Pseudopot_upf::read_pseudo_pswfc(ifstream &ifs)
{
	this->chi.create(this->nwfc, this->mesh);
	for (int i=0;i<nwfc;i++)
	{
		string OrbitalName;
		int BelongToL = 0;
		double occupation = 0.0;
		string dummy;
		ifs >> OrbitalName >> BelongToL >> occupation >> dummy;
		for (int ir = 0;ir < mesh;ir++)
		{
			ifs >> this->chi(i, ir);
		}
	}
	return;
}

void Pseudopot_upf::read_pseudo_rhoatom(ifstream &ifs)
{
	delete[] rho_at;
	this->rho_at = new double[mesh];
	ZEROS(rho_at, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rho_at[ir];
	}
	return;
}

void Pseudopot_upf::read_pseudo_so(ifstream &ifs)
{
       //read soc info from upf, added by zhengdy-soc
       if(!this->has_so) return;
       delete[] nn;
       delete[] jchi;
       delete[] jjj;
       this->nn = new int[nwfc];
       this->jchi = new double[nwfc];
       this->jjj = new double[nbeta];
       ZEROS(nn,nwfc);
       ZEROS(jchi,nwfc);
       ZEROS(jjj,nbeta);
       //RELWFC
       for(int nw=0;nw< nwfc;nw++)
       {
             ifs >> this->els[nw] >>this->nn[nw] >> this->lchi[nw] >> this->jchi[nw] >> this->oc[nw];
             if(this->lchi[nw]-this->jchi[nw]-0.5>1e-7 && this->lchi[nw]-this->jchi[nw]-0.5<1e-7)
             {
                  cout<<"Ignore ADDINFO section"<<endl;
                  this->has_so = 0;
             }
       }
       //RELBETA
       for(int nb = 0;nb < nbeta;nb++)
       {
             ifs >> this->lll[nb] >> this->jjj[nb];
             if(this->lll[nb]-this->jjj[nb]-0.5>1e-7 && this->lll[nb]-this->jjj[nb]-0.5<1e-7)
             {
                  cout<<"Ignore ADDINFO section"<<endl;
                  this->has_so = 0;
             }
       }
       return;
}

/*
//
//	This routine reads from the new UPF file,
//  and the total angual momentum jjj of the beta and jchi of the
//  wave-functions.
void Pseudopot_upf::read_pseudo_addinfo(ifstream &ifs)
{
	for (int i = 0;i < 6;i++)
	{
		ifs.ignore(75, '\n');
	}

	int nb = 0;

	this->nn = new int[nwfc]();
	this->rcut = new double[nwfc]();
	this->rcutus = new double[nwfc]();
	this->epseu = new double[nwfc]();
	this->jchi = new double[nwfc]();
	this->jjj = new double[nbeta]();

	char s[10];

	for (nb = 0; nb < nwfc;nb++)
	{
		ifs >> s >> this->nn[nb] >> this->lchi[nb] >> this->jchi[nb] >> this->oc[nb];
		this->els[nb] = s;
	}


	for (nb = 0; nb < nbeta;nb++)
	{
		ifs >> lll[nb] >> jjj[nb];
	}

	ifs >> xmin >> rmax >> zmesh >> dx;

//	cerr<< " read_pseudo_addinfo Reading pseudo file \n";
}
*/

void Pseudopot_upf::print_pseudo_upf(ofstream &ofs)
{
	TITLE("Pseudopot_upf","print_pseudo_upf");
	ofs << " ==== read_pseudo_upf === " << endl;

	// print header
	ofs << " has_so: " << has_so << endl;
	ofs << " Version number : " << nv << endl;
	ofs << " Element label : " << psd << endl;
	ofs << " pp_type: " << pp_type << endl;
	ofs << " tvanp: " << tvanp << endl;
	ofs << " nlcc: " << nlcc << endl; 
	ofs << " dft: " << dft[0] << " " << dft[1] << " " << dft[2] << " " << dft[3] << endl;
	ofs << " zp: " << zp << endl;
	ofs << " etotps: " << etotps << endl;
	ofs << " ecutwfc: " << ecutwfc << endl;
	ofs << " ecutrho: " << ecutrho << endl;
	ofs << " lmax: " << lmax << endl;
	ofs << " mesh: " << mesh << endl;
	ofs << " nwfc: " << nwfc << endl;
	ofs << " nbeta: " << nbeta << endl;
	for(int i=0; i<nwfc; ++i)
	{
		ofs << " iw=" << i << " els=" << els[i] << " lchi=" << lchi[i] << " oc=" << oc[i] << endl;
	}

	// print PP_MESH
	/*
	ofs.setf(ios::scientific, ios::floatfield);

	for(int i=0; i<100; ++i)
	{
		ofs << setw(25) << r[i] << setw(25) << rab[i] << setw(25) << vloc[i] << setw(25) << rho_at[i] << endl;
	}

	ofs << "chi:" << endl;
	for(int ir=0; ir<100; ++ir)
	{
		for(int i=0; i<nwfc; ++i)
		{
			ofs << " " << chi(i,ir);
		}
		ofs << endl;
	}


	if (nlcc)
	{
		ofs << setw(25) << rho_atc[i] << endl;
	}


	// PRINT PP_NONLOCAL
	ofs << " nd: " << nd << endl;
	ofs << " nbeta: " << nbeta << endl;
	for(int i=0; i<nbeta; ++i)
	{
		ofs << lll[i] << " " << kkbeta[i] << endl;
	}



	ofs << " dion: " << endl;
	for(int i=0; i<nbeta; ++i)
	{
		for(int j=0; j<nbeta; ++j)
		{
			ofs << dion(i,j) << endl;
		}
	}


	for(int i=0; i<nbeta; ++i)
	{
		ofs << " beta: " << i+1 << endl;
		for(int ir=0; ir<kkbeta[i]; ++ir)
		{
			ofs << r[ir] << " " << beta(i,ir) << endl;
		}
	}


	ofs.unsetf(ios::scientific);
    */
	ofs << " End of pseudopot_upf." << endl;

	return;

}


int Pseudopot_upf::read_pseudo_upf201(ifstream &ifs)
{
    string dummy, word;
    //int i, j, idum, ir, ONCVPSP;
    int i, j, ir, ONCVPSP;
    
        
        while (ifs.good())
        {
                ifs >> dummy;
                // We start from PP_Header
                if(dummy=="<PP_HEADER")
                {
                        // Read header
                        READ_VALUE(ifs, word);   // generated
                        READ_VALUE(ifs, word);   // author
                        READ_VALUE(ifs, word);   // date
                        READ_VALUE(ifs, word);   // comment

                        //ifs >> word; //element //replaced with the code below
                        //{
                        //     if(word == "element=\"")
                        //     {
                        //        ifs >> word;
                        //        get_char(word);
                        //        this->psd = word.substr(0,Number[0]);
                        //     }
                        //     else
                        //     {
                        //        get_char(word);
                        //        this->psd = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
                        //     }
                        //     cout << "psd = " << this->psd << endl;
                        //}
                        READ_VALUE(ifs, word);   //{zws add
                        stringstream wdsstream(word);
                        getline(wdsstream,this->psd,'"');
                        getline(wdsstream,this->psd,'"'); //}zws
                        //cout << " this->psd = " << this->psd << endl;

                        ifs >> word;   // pseudo_type
                        {
                             if(word == "pseudo_type=\"")
                             {
                                ifs >> word;
                                get_char(word);
                                this->pp_type = word.substr(0,Number[0]);
                             }
                             else
                             {
                                get_char(word);
                                this->pp_type = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
                             }
                             //cout << "pp_type = " << this->pp_type << endl;
                        }

                             if(pp_type!="NC") 
                             {
                                //cout << " pp_type=" << pp_type << endl;
                                WARNING_QUIT("Pseudopot_upf::read_pseudo_header","unknown pseudo type");
                             }

                        READ_VALUE(ifs, word);   // relativistic
                        READ_VALUE(ifs, word);   // is_ultrasoft
                        if ( word.find("\"T\"") < word.length() ) // zws add 20160108
                        {
                        	cout << "\n WARNING: ULTRASOFT PSEUDOPOTENTIAL IS NOT SUPPORTED !!! \n" << endl;
                        }
                        READ_VALUE(ifs, word);   // is_paw
                        if ( word.find("\"T\"") < word.length() )
                        {
                        	cout << "\n WARNING: PAW PSEUDOPOTENTIAL IS NOT SUPPORTED !!! \n" << endl;
                        }

                        READ_VALUE(ifs, word);   // is_coulomb
                        ifs >> word;   // has_so
                             string so;
                             //cout << "word = " << word << endl;
                             {
                                  if(word == "has_so=\"")
                                  {
                                     ifs >> word;
                                     get_char(word);
                                     so = word.substr(0,Number[0]);
                                  }
                                  else
                                  {
                                     get_char(word);
                                     so = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
                                  }
                                  //cout << "so = " << so << endl;
                              } 

                             if (so == "T")
                             {
                                 this->has_so = true;
                             }
                             else
                             {
                                 this->has_so = false;
                             }

                        READ_VALUE(ifs, word);   // has_wfc
                        READ_VALUE(ifs, word);   // has_gipaw
                        string nlc;
                        //char p[13] = "paw_as_gipaw";
                        ifs >> word;             // paw_as_gipaw?
                        //cout << "word.substr(0,30) = " << word.substr(0,30) << "."<< endl;
                             if( word.substr(0,13) == "paw_as_gipaw" )
                             {
                                ONCVPSP = 0;
                                ifs >> word;     // core_correction
                                if(word == "core_correction=\"")
                                {
                                     ifs >> word;
                                     get_char(word);
                                     nlc = word.substr(0,Number[0]);
                                }
                                else
                                {
                                     get_char(word);
                                     nlc = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
                                }

                             }
                             else
                             {
                                ONCVPSP = 1; // Generated using ONCVPSP code by D. R. Hamann, SG15 DOJO
                                if(word == "core_correction=\"")
                                {
                                     ifs >> word;
                                     get_char(word);
                                     nlc = word.substr(0,Number[0]);
                                }
                                else
                                {
                                     get_char(word);
                                     nlc = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
                                }

                             }

                             //cout << "nlc = " << nlc << endl;

                             if (nlc == "T")
                             {
                                 this->nlcc = true;
                             }
                             else
                             {
                                 this->nlcc = false;
                             }

                        READ_VALUE(ifs, word);   // functional
                        //cout << "word = " << word << endl;
//                        this->dft[0]="SLA";
//                        this->dft[1]="PZ";
//                        this->dft[2]="NOGX";
//                        this->dft[3]="NOGC";
                        string funcstr;  //{zws 01-06-16
                        wdsstream.str("");
                        wdsstream.clear();
                        wdsstream << word;
                        for ( int idft = 0; idft < 2; idft++)
                        {
                        	getline(wdsstream,funcstr,'"');
                        }
                        wdsstream.str("");
                        wdsstream.clear();
					    wdsstream << funcstr;
					    //cout << " DFT                  : " ;
                        for( int idft = 0; idft < 4; idft++ )
                        {
                        	getline(wdsstream,dft[idft],'-');
                        	//cout << dft[idft] << " " ;
                        }
                        //cout << endl;  //}zws 01-06-16


                        //getline(ifs, word);
                        //cout << "word       = " << word << endl;
                        do  // zws begin from cutoff (rho/wfc) 
                        {
                            getline(ifs, word);
                            //cout << "word       = " << word << endl;
                            word.erase(0,word.find_first_not_of(" ") );
                            word.erase(word.find_last_not_of(" ")+1 );
                            //word = trim(word);
                            //cout << "trim(word) = " << word << endl;
                            get_char(word);
                            //cout << " Number = " << Number[0] << ", " << Number[1] << endl;
                            //cout << word.substr(0,Number[0])  << "__" << word.substr(Number[0]+1, Number[1]-Number[0]-1) << endl;
                            dummy = word.substr(0,Number[0]) ;
                            //cout << " dummy = " << dummy << endl;
                            if      ( dummy == "z_valence=" ) {
                                     this->zp = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                                     //cout << "zp = " << this->zp << endl;
                            }
                            else if ( dummy == "total_psenergy=" ) {
                                     this->etotps = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                                     //cout << "etotps = " << this->etotps << endl;
                            }
                            else if ( dummy == "rho_cutoff=" ){
                                     //cout << " rho_cutoff = " << word.substr(Number[0]+1, Number[1]-Number[0]-1) << endl;
                            }
                            else if ( dummy == "wfc_cutoff=" ) {
                                     //cout << " wfc_cutoff = " << word.substr(Number[0]+1, Number[1]-Number[0]-1) << endl;
                            }
                            else if ( dummy == "l_max=" ) {
                                     this->lmax = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                            }
                            else if ( dummy == "mesh_size=" ) {
                                     this->mesh = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                            }
                            else if ( dummy == "number_of_wfc=" ) {
                                     this->nwfc = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                            }
                            else if ( dummy == "number_of_proj=" ) {
                                     this->nbeta = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                            }
                            //break;

                        }while( word.substr(word.length()-1, 1) !=">" ); 
                        //cout << "word.substr(word.length()-1, 1)=" <<  word.substr(word.length()-1, 1)  << endl;
                        //exit(0);


                        //ifs >> word;   // zp
                        ////cout << "word = " << word << endl;
                        //{
                        //     if(word == "z_valence=\"")
                        //     {
                        //        ifs >> word;
                        //        get_char(word);
                        //        this->zp = atoi(word.substr(0,Number[0]).c_str());
                        //     }
                        //     else
                        //     {
                        //        get_char(word);
                        //        this->zp = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                        //     }
                        //     //cout << "zp = " << this->zp << endl;
                        //}

                        //ifs >> word;   // total_psenergy
                        //{
                        //     if(word == "total_psenergy=\"")
                        //     {
                        //        ifs >> word;
                        //        get_char(word);
                        //        this->etotps = atof(word.substr(0,Number[0]).c_str());
                        //     }
                        //     else
                        //     {
                        //        get_char(word);
                        //        this->etotps = atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                        //     }
                        //     //cout << "etotps = " << this->etotps << endl;
                        //}
                        ////cout << " word (total_psenergy) = " << word << endl;
                  

                        //if(ONCVPSP == 0)    //zws modify 20160108
                        //{
                        //	READ_VALUE(ifs, word);   // wfc_cutoff
                        //	//cout << "word = " << word << endl;
                        //}
                        //READ_VALUE(ifs, word); // rho_cutoff
                        //cout << "word (cutoff) = " << word << endl;


                        //ifs >> word;             // lmax
                        ////cout << "word (lmax) = " << word << endl;
                        //{
                        //        if(word == "l_max=\"")
                        //        {
                        //             ifs >> word;
                        //             get_char(word);
                        //             this->lmax = atoi(word.substr(0,Number[0]).c_str());
                        //        }
                        //        else
                        //        {
                        //             get_char(word);
                        //             this->lmax = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                        //        }

                        //}

                        ////cout << "lmax = " << this->lmax << endl;

                        //if(ONCVPSP == 0)
                        //{
                        //   READ_VALUE(ifs, word);   // l_max_rho
                        //}

                        //READ_VALUE(ifs, word);   // l_local

                        //ifs >> word;   // mesh_size
                        ////cout << "word (mesh) = " << word << endl;
                        //{
                        //     if(word == "mesh_size=\"")
                        //     {
                        //             ifs >> word;
                        //             get_char(word);
                        //             this->mesh = atoi(word.substr(0,Number[0]).c_str());
                        //     }
                        //     else
                        //     {
                        //             get_char(word);
                        //             this->mesh = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                        //     }
                        //     //cout << "mesh = " << this->mesh << endl;
                        //}
 


                        //ifs >> word;  // number_of_wfc
                        ////cout << "word = " << word << endl;
                        //{
                        //     if(word == "number_of_wfc=\"")
                        //     {
                        //             ifs >> word;
                        //             get_char(word);
                        //             this->nwfc = atoi(word.substr(0,Number[0]).c_str());

                        //     }
                        //     else
                        //     {
                        //             get_char(word);
                        //             this->nwfc = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                        //     }
                        //     //cout << "nwfc = " << this->nwfc << endl;
                        //}
                        //     
                        //ifs >> word;   // number_of_proj
                        ////cout << "word = " << word << endl;
                        //{
                        //     if(word == "number_of_proj=\"")
                        //     {
                        //             ifs >> word;
                        //             get_char(word);
                        //             this->nbeta = atoi(word.substr(0,Number[0]).c_str());

                        //     }
                        //     else
                        //     {
                        //             get_char(word);
                        //             this->nbeta = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                        //     }
                        //     //cout << "nbeta = " << this->nbeta << endl;
                        //}
                        

                        // READ Mesh
                        if(ONCVPSP == 0)
                        {
                           SCAN_BEGIN(ifs, "<PP_MESH");
                        }
                        else
                        {
                           SCAN_BEGIN(ifs, "<PP_MESH>");
                        }
    
                        assert(mesh>0);
                        if(ONCVPSP == 0)
                        {
                           ifs >> word;             // dx
                           ifs >> word;             // mesh
                           ifs >> word;             // xmin
                           ifs >> word;             // rmax
                           ifs >> word;             // zmesh
                        }

                        SCAN_BEGIN(ifs, "<PP_R"); 
                        READ_VALUE(ifs, word);    // type  size  columns

//                        double  rmesh0 = 1;    //{zws add160108 delete160328
//                        int 	nmeshdel = 0;
//                        ifs >> rmesh0;
//                        if ( abs(rmesh0) < 1.0e-15 )
//                        {
//                            mesh       -= 1;
//                            nmeshdel   += 1;
//                        }
//                        cout << " mesh =" << mesh << endl;
//                    	if (mesh%2 == 0)
//                    	{
//                    	    mesh     -= 1;
//                    	    nmeshdel += 1;
//                    	}    //}zws add 20160108
//                    	cout << " nmeshdel =" << nmeshdel << endl;


                        delete[] r;
                        delete[] rab;
                        this->r = new double[mesh];
                        this->rab = new double[mesh];
                        ZEROS(r,mesh);
                        ZEROS(rab,mesh);


//                        if (nmeshdel == 0)    //{zws add160108 delete160328
//                        {
//                            this->r[0] = rmesh0;
//                            for (ir = 1;ir < mesh;ir++)
//                            {
//                                ifs >> this->r[ir];
//                            }
//                        }
//                        else
//                        {
//                            for ( int idel=0; idel < nmeshdel-1; idel++)
//                        	{
//                            	cout << "skip " << nmeshdel << "grid point(s) in PP mesh" << endl;
//                        	    double	tmpdel;
//                        	    ifs >> tmpdel;
//                        	}
//                            for (ir = 0;ir < mesh;ir++)
//                            {
//                                 ifs >> this->r[ir];
//                            }
//                        }    //}zws 20160108
						for (ir = 0;ir < mesh;ir++)
						{
							 ifs >> this->r[ir];
						}
                        SCAN_END(ifs, "</PP_R>");

                        SCAN_BEGIN(ifs, "<PP_RAB");
                        READ_VALUE(ifs, word);    // type size columns

//                        for ( int idel=0; idel < nmeshdel; idel++)    //{zws add 20160108
//                    	{
//                    	    double	tmpdel;
//                    	    ifs >> tmpdel;
//                    	}    //}zws add 20160108
                        for (ir = 0;ir < mesh;ir++)
                        {
                             ifs >> this->rab[ir];
                        }
                        SCAN_END(ifs, "</PP_RAB>");
                        SCAN_END(ifs, "</PP_MESH>");

                        // READ NLCC
                        if (this->nlcc)
                        {
                            SCAN_BEGIN(ifs, "<PP_NLCC");
                            READ_VALUE(ifs, word);    // type size columns
                            
                            assert(mesh>0);
                            delete[] rho_atc;
                            this->rho_atc = new double[mesh];
                            ZEROS(rho_atc, mesh);
//                            for ( int idel=0; idel < nmeshdel; idel++)    //{zws add 20160108
//                        	{
//                        	    double	tmpdel;
//                        	    ifs >> tmpdel;
//                        	}    //}zws add 20160108
                            for (ir = 0;ir < mesh;ir++)
                            {
                                 ifs >> this->rho_atc[ir];
                            }
                            SCAN_END(ifs, "</PP_NLCC>");

                        }

                        // READ VLOCAL
                        SCAN_BEGIN(ifs, "<PP_LOCAL");
                        READ_VALUE(ifs, word);    // type size columns
                        
                        assert(mesh>0);
                        delete[] vloc;
                        this->vloc = new double[mesh];
                        ZEROS(vloc, mesh);

//                        for ( int idel=0; idel < nmeshdel; idel++)    //{zws add 20160108
//                    	{
//                    	    double	tmpdel;
//                    	    ifs >> tmpdel;
//                    	}    //}zws add 20160108
                        for (ir = 0;ir < mesh;ir++)
                        {
                             ifs >> this->vloc[ir];
                        }

                        SCAN_END(ifs, "</PP_LOCAL>");

                        // READ NONLOCAL
                        SCAN_BEGIN(ifs, "<PP_NONLOCAL>");

                        delete[] kkbeta;
                        delete[] lll;
                        this->kkbeta = new int[nbeta];
                        this->lll = new int[nbeta];
                        this->beta.create(nbeta , mesh);
                        this->dion.create(nbeta , nbeta);

                        for(i=0;i<nbeta;i++)
                        {
                            ifs >> word;  //number
                            ifs >> word;  //type
                            //cout << "word = " << word << endl;
                                 if(word == "type=\"")
                                 {
                                    ifs >> word;
                                 }
                            ifs >> word;  //size
                            //cout << "word = " << word << endl;
                                 if(word == "size=\"")
                                 {
                                    ifs >> word;
                                 }

                            ifs >> word;  //columns
                            //cout << "word = " << word << endl;
                                 if(word == "columns=\"")
                                 {
                                    ifs >> word;
                                 }

                            ifs >> word;  //index
                            //cout << "word = " << word << endl;
                            {
                                 if(word == "index=\"")
                                 {
                                     ifs >> word;
                                     get_char(word);
                                     //idum = atoi(word.substr(0,Number[0]).c_str());
                                 }
                                 else
                                 {
                                     get_char(word);
                                     //idum = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                                 }
                                 //cout << "idum = " << idum << endl;
                            }
  
                            if(ONCVPSP == 0)
                            {
                               ifs >> word;  //label
                            }

                            ifs >> word;  //angular_momentum
                            //cout << "word = " << word << endl;
                            {
                                 if(word == "angular_momentum=\"")
                                 {
                                    ifs >> word;
                                    get_char(word);
                                    this->lll[i] = atoi(word.substr(0,Number[0]).c_str());
            
                                 }
                                 else
                                 {
                                    get_char(word);
                                    this->lll[i] = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                                 }
                                 //cout << "lll[i] = " << this->lll[i] << endl;
                            }

                            ifs >> word;  //cutoff_radius_index
                            //cout << "word = " << word << endl;
                            {
                                 if(word == "cutoff_radius_index=\"")
                                 {
                                    ifs >> word;
                                    get_char(word);
                                    this->kkbeta[i] = atoi(word.substr(0,Number[0]).c_str());

                                 }
                                 else
                                 {
                                    get_char(word);
                                    this->kkbeta[i] = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                                 }
                                 //cout << "kkbeta[i] = " << this->kkbeta[i] << endl;
                            }
  
                            if(ONCVPSP ==0) 
                            {
                               ifs >> word;  //cutoff_radius
                               ifs >> word;  //ultrasoft_cutoff_radius
                            }
                            else
                            {
                                READ_VALUE(ifs, word); // cutoff_radius
                            }


//                            for ( int idel=0; idel < nmeshdel; idel++)    //{zws add 20160108
//                        	{
//                        	    double	tmpdel;
//                        	    ifs >> tmpdel;
//                        	}    //}zws add 20160108
                            for (ir=0;ir<mesh;ir++)
                            {
                                ifs >> this->beta(i, ir);

                            }
       
                            ifs >> word;  //number
                            //cout << "word = " << word << endl;

                        }
                        
                        // READ DIJ
                        SCAN_BEGIN(ifs, "<PP_DIJ");
                        READ_VALUE(ifs, word);  // type size columns
                        
                        this->nd = nbeta * nbeta;
                        for(i=0;i<nbeta;i++)
                            for(j=0;j<nbeta;j++)
                            {
                                ifs >> dion(i,j);
								if ( i != j  && dion(i,j) != 0.0 )
								{
									cout << " error: for i != j, Dij of Pseudopotential must be 0.0 " << endl;
									exit(1);
								}
                            }
                        SCAN_END(ifs, "</PP_DIJ>");

                        SCAN_END(ifs, "</PP_NONLOCAL>");

                        // READ PSWFC
                        SCAN_BEGIN(ifs, "<PP_PSWFC>");

                        delete[] els;
                        delete[] lchi;
                        delete[] oc;
                        this->els = new string[nwfc];
                        this->lchi = new int[nwfc];
                        this->oc = new double[nwfc];
                        ZEROS(lchi, nwfc); // angular momentum of each orbital
                        ZEROS(oc, nwfc);//occupation of each orbital

                        //cout << " nwfc = " << this->nwfc << ", mesh = " << this->mesh << endl;
                        this->chi.create(this->nwfc, this->mesh);
                        for (i=0;i<nwfc;i++)
                        {
                             ifs >> word;  // number
                             //cout << " word (number) = " << word << endl;
                             ifs >> word;  // type
                             ifs >> word;  // size 
                             //cout << " word (size) = " << word << endl;
                             {
                                 if(word == "size=\"")
                                 {
                                    ifs >> word;
                                    word = "\"" + word ;
                                 }
                             }
                             //cout << " word (size) = " << word << endl;
                             ifs >> word;  // columns
                             //cout << " word (columns) = " << word << endl;
                             ifs >> word;  // index

                             ifs >> word;  // occupation
                             {
                                 if(word == "occupation=\"")
                                 {
                                    ifs >> word;
                                    word = "\"" + word ;
                                 }
                                 get_char(word);
                                 oc[i] = atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                                 //cout << "oc[i] = " << oc[i] << endl;
                             }

                             ifs >> word;  // pseudo_energy
                             //cout << " word (pseudo) = " << word << endl;
                             {
                                 if(word == "pseudo_energy=\"")
                                 {
                                    ifs >> word;
                                    word = "\"" + word ;
                                 }
                                 get_char(word);
                                 //cout << "word pseudo_energy = " << word << endl;
                             }

                             ifs >> word;  // label
                             //cout << "word = " << word << endl;
                             {
                                 if(word == "label=\"")
                                 {
                                    ifs >> word;
                                    word = "\"" + word ;
                                 }
                                 get_char(word);
                                 els[i] = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
                                 //cout << "els[i] = " << els[i] << endl;
                             }


                             ifs >> word;  // l
                             //cout << " word (l) = " << word << endl;
                             {
                                 if(word == "l=\"")
                                 {
                                    ifs >> word;
                                    word = "\"" + word ;
                                 }
                                 get_char(word);
                                 lchi[i] = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
                                 //cout << " lchi[i] = " << lchi[i] << endl;
                             }

                             ifs >> word; // >
                             //cout << " word (>) = " << word << endl;
                             if ( word !=  ">" )
                             {
                            	 cout << " error: bad end while reading CHI" << i <<  " of PSWFC" << endl;
                            	 exit(1);
                             }

                             //ifs >> word; // n
                             //cout << "word n = " << word << endl;
                             //ifs >> word; // pseudo_energy
                             //cout << "word pseudo_energy = " << word << endl;
                             //ifs >> word; // cutoff_radius
                             //cout << "word cutoff_radius = " << word << endl;
                             //ifs >> word; // ultrasoft_cutoff_radius
                             //cout << "word ultrasoft_cutoff_radius = " << word << endl;


//                           for ( int idel=0; idel < nmeshdel; idel++)    //{zws add 20160108
//                         	 {
//                         	     double	tmpdel;
//                                 ifs >> tmpdel;
//                         	 }    //}zws add 20160108
                             for (ir = 0;ir < mesh;ir++)
                             {
                                  ifs >> this->chi(i, ir);
                             }
                             ifs >> word;  // number
                        }

                        SCAN_END(ifs, "</PP_PSWFC>");

                        // READ RHOATOM
                        SCAN_BEGIN(ifs, "<PP_RHOATOM");
                        READ_VALUE(ifs, word); // type size columns
    
                        delete[] rho_at;
                        this->rho_at = new double[mesh];
                        ZEROS(rho_at, mesh);

//                        for ( int idel=0; idel < nmeshdel; idel++)    //{zws add 20160108
//                        {
//                            double	tmpdel;
//                            ifs >> tmpdel;
//                        }    //}zws add 20160108
                        for (ir = 0;ir < mesh;ir++)
                        {
                             ifs >> this->rho_at[ir];
                        }
                        SCAN_END(ifs, "</PP_RHOATOM>");

			SCAN_BEGIN(ifs, "<PP_SPIN_ORB>");
//added by zhengdy-soc
			delete[] this->jchi;
			delete[] this->jjj;
			delete[] this->nn;
			this->jchi = new double [nwfc];
			this->jjj = new double [nbeta];
			this->nn = new int [nwfc];
			ZEROS(jchi,nwfc);
			ZEROS(jjj,nbeta);
			ZEROS(nn,nwfc);
			for(int round=0;round<2;round++){
				ifs>>word;
				if(word=="<PP_RELBETA.1")
				{
					for(int nb = 0;nb<nbeta;nb++)
					{
						if(nb!=0) ifs>>word; //RELBETA
						ifs>>word;           //index
						ifs>>word;           //lll
						get_char(word);
						this->lll[nb] = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
						ifs>>word;           //jjj
						get_char(word);
						this->jjj[nb] = atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
					}
				}
				else if(word=="<PP_RELWFC.1")
				{
					if(round==0){
						for(int nw = 0;nw<nwfc;nw++)
						{
							if(nw!=0) ifs>>word;     //RELWFC
							ifs>>word;               //index
							ifs>>word;               //els
							get_char(word);
							this->els[nw]= word.substr(Number[0]+1,(Number[1]-Number[0]-1));
							ifs>>word;               //nn
							get_char(word);
							this->nn[nw]= atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;               //lchi
							get_char(word);
							this->lchi[nw]= atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;               //jchi
							get_char(word);
							this->jchi[nw]= atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;               //oc
							get_char(word);
							this->oc[nw]= atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());

						}
					}
					else{
						for(int nw = 0;nw<nwfc;nw++)
						{
							if(nw!=0) ifs>>word;//RELWFC
							ifs>>word;          //index
							ifs>>word;          //lchi
							get_char(word);
							this->lchi[nw]= atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;          //jchi
							get_char(word);
							this->jchi[nw]= atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;
							get_char(word);
							this->nn[nw]= atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());

						}
					}
				}
				else if(round==0){
					this->has_so = 0;
//					cout<<"ignore SPIN_ORB part!"<<endl;
					break;
				}
			}
                        SCAN_END(ifs, "</PP_SPIN_ORB>");

                    	if (mesh%2 == 0) //{zws add 20160328
                    	{
                    	    mesh     -= 1;
                    	}
                    	//cout << " mesh =" << mesh << endl; //}


                        SCAN_END(ifs, "</UPF>");
                        break;
                }
        }
    return 0;
}


void Pseudopot_upf::get_char( string ss)
{
    int i, q;
    //char b[1]; //LiuXh 20171109
    char b='\"'; //LiuXh 20171109
    q =0;
    //strcpy(b,"\""); //LiuXh 20171109

    for(i=0;i<200;i++)
    {
        //if(ss[i]== b[0]) //LiuXh 20171109
        if(ss[i]== b) //LiuXh 20171109
        {
           Number[q] = i;
           q++;
        }

    }

    return;
}

/*string Pseudopot_upf::get_string( char ss[])
{
    int n[2], i, q;
    char b[1], ssss[200];
    string sss;
    q =0;
    strcpy(b,"\"");

    for(i=0;i<200;i++)
    {
        if(ss[i]== b[0])
        {
           n[q] = i;
           q++;
        }

    }

    for(i=0;i<(n[1]-n[0]-1);i++)
    {

        ssss[i] = ss[(n[0]+i+1)];
    }
    cout << "ssss = " << ssss << endl;
    sss.assign(ssss);
    cout << "sss = " << sss << endl;
    return sss;

}

int Pseudopot_upf::get_int( char ss[])
{
    int n[2], i, q, m;
    char b[1],sss[200];
    q =0;
    strcpy(b,"\"");

    for(i=0;i<200;i++)
    {
        if(ss[i]== b[0])
        {
           n[q] = i;
           q++;
        }

    }

    for(i=0;i<(n[1]-n[0]-1);i++)
    {

        sss[i] = ss[(n[0]+i+1)];
    }
    m = atoi(sss);

    return m;

}

double Pseudopot_upf::get_double( char ss[])
{
    int n[2], i, q;
    double l;
    char b[1],sss[200];
    q =0;
    strcpy(b,"\"");

    for(i=0;i<200;i++)
    {
        if(ss[i]== b[0])
        {
           n[q] = i;
           q++;
        }

    }

    for(i=0;i<(n[1]-n[0]-1);i++)
    {

        sss[i] = ss[(n[0]+i+1)];
    }
    l = atof(sss);

    return l;

}*/

int Pseudopot_upf::average_p()
{
	int error = 0;
	if(!this->has_so && LSPINORB) {error++; cout<<"warning_quit! no soc upf used for lspinorb calculation, error!"<<endl; return error;}
		//WARNING_QUIT("average_p", "no soc upf used for lspinorb calculation, error!");
	if(!this->has_so || LSPINORB) return error; 
	int new_nbeta = 0; //calculate the new nbeta
	for(int nb=0; nb< this->nbeta; nb++)
	{
		new_nbeta++;
		if(this->lll[nb] != 0 && abs(this->jjj[nb] - this->lll[nb] - 0.5) < 1e-6) //two J = l +- 0.5 average to one
			new_nbeta--;
	}
	this->nbeta = new_nbeta;
	matrix dion_new;
	dion_new.create(this->nbeta, this->nbeta);

	int old_nbeta=-1;
	for(int nb=0; nb<this->nbeta; nb++)
	{
		old_nbeta++;
		int l = this->lll[old_nbeta];
		int ind=0, ind1=0;
		if(l != 0)
		{
//			cout<<this->jjj[old_nbeta] <<" "<< this->lll[old_nbeta]<<" "<<abs(this->jjj[old_nbeta] - this->lll[old_nbeta] + 0.5)<<endl;
//			cout<<this->jjj[old_nbeta+1] <<" "<< this->lll[old_nbeta+1]<<" "<<abs(this->jjj[old_nbeta+1] - this->lll[old_nbeta+1] + 0.5)<<endl;
			if(abs(this->jjj[old_nbeta] - this->lll[old_nbeta] + 0.5) < 1e-6)
			{
				if(abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]-0.5)>1e-6) 
					WARNING_QUIT("average_p", "error beta function 1 !");
				ind = old_nbeta +1;
				ind1 = old_nbeta;
			}
			else
			{
				if(abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]+0.5)>1e-6)
					WARNING_QUIT("average_p", "error beta function 2 !");
				ind = old_nbeta;
				ind1 = old_nbeta +1;
			}
			double vion1 = ((l+1.0) * this->dion(ind,ind) + l * this->dion(ind1,ind1)) / (2.0*l+1.0);
			//average beta (betar)
			for(int ir = 0; ir<this->mesh;ir++)
			{
				this->beta(nb, ir) = 1.0 / (2.0 * l + 1.0) * 
						( (l + 1.0) * sqrt(this->dion(ind,ind) / vion1) *
						this->beta(ind, ir) + 
						l * sqrt(this->dion(ind1,ind1) / vion1) *
						this->beta(ind1, ir) ) ;
			}
			//average the dion matrix
			this->dion(nb, nb) = vion1;
			old_nbeta++;	
		}
		else{
			for(int ir = 0; ir<this->mesh;ir++)
				this->beta(nb, ir) = this->beta(old_nbeta, ir);
			this->dion(nb, nb) = this->dion(old_nbeta, old_nbeta);
		}
		this->lll[nb] = this->lll[old_nbeta]; //reset the lll index, ignore jjj index
	}
	//store the old dion and then recreate dion 
	for(int i=0;i<this->nbeta; i++)
		for(int j=0;j<this->nbeta;j++)
			dion_new(i,j) = this->dion(i,j);
	this->dion = dion_new;
//	this->dion.create(this->nbeta, this->nbeta);
//	for(int i=0;i<this->nbeta; i++)
//		for(int j=0;j<this->nbeta;j++)
//			this->dion(i,j) = dion_new(i,j);
	
	int new_nwfc = 0;
	for(int nb=0; nb<this->nwfc; nb++)
	{
		new_nwfc++;
		if(this->lchi[nb] != 0 && abs(this->jchi[nb] - this->lchi[nb] - 0.5)<1e-6)
			new_nwfc--;
	}
	this->nwfc = new_nwfc;
	int old_nwfc=0;
	for(int nb=0; nb<this->nwfc; nb++)
	{
		old_nwfc++;
		int l = this->lchi[old_nwfc];
		int ind=0, ind1=0;
		if(l!=0)
		{
			if(abs(this->jchi[old_nwfc] - this->lchi[old_nwfc] + 0.5) < 1e-6)
			{
				if(abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]-0.5)>1e-6) 
				{error++; cout<<"warning_quit! error chi function 1 !"<<endl; return error;}
//					WARNING_QUIT("average_p", "error chi function 1 !");
				ind = old_nwfc +1;
				ind1 = old_nwfc;
			}
			else
			{
				if(abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]+0.5)>1e-6)
				{error++; cout<<"warning_quit! error chi function 2 !"<<endl; return error;}
//					WARNING_QUIT("average_p", "error chi function 2 !");
				ind = old_nwfc;
				ind1 = old_nwfc +1;
			}
			//average chi
			for(int ir = 0; ir<this->mesh;ir++)
			{
				this->chi(nb, ir) = 1.0 / (2.0 * l + 1.0) *
					( (l+1.0)*this->chi(ind,ir) + (l*this->chi(ind1,ir)) );
				old_nwfc++;
			}
		}
		else{
			for(int ir = 0; ir<this->mesh;ir++)
				this->chi(nb, ir) = this->chi(old_nwfc, ir);
		}
		this->lchi[nb] = this->lchi[old_nwfc]; //reset lchi index
	}
	this->has_so = 0;	
	return error;
}
