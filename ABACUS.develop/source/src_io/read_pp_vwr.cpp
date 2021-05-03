#include "read_pp.h"

//----------------------------------------------------------
// This code is used to read in vwr pseudopotential format,
// Now we only use LDA, so if PBE or other functionals are used, 
// one needs to change the following code. The vwr format 
// needs we to generate NL projectors by ourself. 
// One way to check this is correct
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
		ofs_running << " Mesh number - 1, we need odd number, \n this may affect some polar atomic orbitals." << endl;
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
		}


// In pw they did this:
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
