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
int Pseudopot_upf::read_pseudo_vwr(std::ifstream &ifs, Atom_pseudo& pp)
{
	GlobalV::ofs_running << " -------------------------------------------------" << std::endl;
	std::cout << " READ IN VWR TYPE PSEUDOPOTENTIALS." << std::endl;
	GlobalV::ofs_running << " Read in vwr type pseudopotentials " << std::endl;


	// --------------------------------------
	// (1) read in data 
	// --------------------------------------
	pp.xc_func="PZ";
	pp.pp_type="NC";
    pp.tvanp = false;

    // (1) read in mesh
	std::string value;
	size_t length=0;
	ifs >> value; length = value.find(","); value.erase(length,1);
	pp.mesh = std::atoi( value.c_str() );
	//the mesh should be odd, which is forced in Simpson integration
	this->mesh_changed = false;
	if(pp.mesh%2==0) 
	{
		pp.mesh=pp.mesh-1;
		this->mesh_changed = true;
		GlobalV::ofs_running << " Mesh number - 1, we need odd number, \n this may affect some polar atomic orbitals." << std::endl;
	}
	GlobalV::ofs_running << std::setw(15) << "MESH" << std::setw(15) << pp.mesh << std::endl;
	// (2) read in nlcc: nonlinear core correction
	ifs >> value; length = value.find(","); value.erase(length,1);
	pp.nlcc = std::atoi( value.c_str() );
	GlobalV::ofs_running << std::setw(15) << "NLCC" << std::setw(15) << pp.nlcc << std::endl;
	// (3) iatom : index for atom 
	ifs >> value; length = value.find(","); value.erase(length,1);
	pp.psd = value;
	GlobalV::ofs_running << std::setw(15) << "ATOM" << std::setw(15) << pp.psd << std::endl;
	// (4) valence electron number
	ifs >> value; length = value.find(","); value.erase(length,1);
	pp.zv = std::stod( value );
	GlobalV::ofs_running << std::setw(15) << "Z(VALENCE)" << std::setw(15) << pp.zv << std::endl;
	// (5) spd_loc, which local pseudopotential should I choose
	ifs >> value; length = value.find(","); value.erase(length,1);
	spd_loc = std::atoi( value.c_str() );
	GlobalV::ofs_running << std::setw(15) << "LOC(spd)" << std::setw(15) << spd_loc << std::endl;
	// (6) read in the occupations
	std::vector<double> tmp_oc(3, 0.0);
	ifs >> value; length = value.find(","); value.erase(length,1);
	tmp_oc[0]= std::atoi( value.c_str() );
	ifs >> value; length = value.find(","); value.erase(length,1);
	tmp_oc[1]= std::atoi( value.c_str() );
	ifs >> value; length = value.find(","); value.erase(length,1);
	tmp_oc[2]= std::atoi( value.c_str() );
	GlobalV::ofs_running << std::setw(15) << "OCCUPATION" << std::setw(15) << tmp_oc[0] 
	<< std::setw(15) << tmp_oc[1] << std::setw(15) << tmp_oc[2] << std::endl;
	// (7) spin orbital
	ifs >> pp.has_so;


	// label to count the projector or atomic wave functions	
	getline(ifs,value);
	int iref_s, iref_p, iref_d;
	ifs >> iref_s >> iref_p >> iref_d;
	GlobalV::ofs_running << std::setw(15) << "Vnl_USED" << std::setw(15) << iref_s 
	<< std::setw(15) << iref_p << std::setw(15) << iref_d << std::endl;
	if(spd_loc==1) { iref_s=0;
	} else if(spd_loc==2) { iref_p=0;
	} else if(spd_loc==3) { iref_d=0;
}
	ifs >> iTB_s >> iTB_p >> iTB_d;
	GlobalV::ofs_running << std::setw(15) << "Orb_USED" << std::setw(15) << iTB_s 
	<< std::setw(15) << iTB_p << std::setw(15) << iTB_d << std::endl;
	
	
	// calculate the number of wave functions
	pp.nchi = 0;
	if(iTB_s) { ++pp.nchi;
}
	if(iTB_p) { ++pp.nchi;
}
	if(iTB_d) { ++pp.nchi;
}
	GlobalV::ofs_running << std::setw(15) << "NWFC" << std::setw(15) << pp.nchi << std::endl;
	// allocate occupation number array for wave functions
	pp.oc = std::vector<double>(pp.nchi, 0.0);
	pp.els = std::vector<std::string>(pp.nchi, "");
	// set the value of occupations
	pp.lchi = std::vector<int>(pp.nchi, 0);
	int iwfc=0;
	if(iTB_s){pp.oc[iwfc]=tmp_oc[0];pp.lchi[iwfc]=0;pp.els[iwfc]="S";++iwfc;}
	if(iTB_p){pp.oc[iwfc]=tmp_oc[1];pp.lchi[iwfc]=1;pp.els[iwfc]="P";++iwfc;}
	if(iTB_d){pp.oc[iwfc]=tmp_oc[2];pp.lchi[iwfc]=2;pp.els[iwfc]="D";++iwfc;}
	getline(ifs,value);


	// global variables that will be used
	// in other classes.
	pp.r = std::vector<double>(pp.mesh, 0.0);
	pp.rab = std::vector<double>(pp.mesh, 0.0);
	pp.vloc_at = std::vector<double>(pp.mesh, 0.0);
	pp.rho_at = std::vector<double>(pp.mesh, 0.0);
	pp.rho_atc = std::vector<double>(pp.mesh, 0.0);
	// local variables in this function
    std::vector<double> vs = std::vector<double>(pp.mesh, 0.0); // local pseudopotential for s, unit is Hartree
    std::vector<double> vp = std::vector<double>(pp.mesh, 0.0); // local pseudopotential for p
    std::vector<double> vd = std::vector<double>(pp.mesh, 0.0); // local pseudopotential for d
    std::vector<double> ws = std::vector<double>(pp.mesh, 0.0); // wave function for s
    std::vector<double> wp = std::vector<double>(pp.mesh, 0.0); // wave function for p
    std::vector<double> wd = std::vector<double>(pp.mesh, 0.0); // wave function for d
	std::string line;
	if(spd_loc>0 && pp.nlcc==0)
	{
		for(int ir=0; ir<pp.mesh; ++ir)
		{
			// it's an interesting question whether
			// ws[ir] has 1/sqrt(4pi)
			ifs >> pp.r[ir] >> vs[ir] >> vp[ir] >> vd[ir] 
				>> ws[ir] >> wp[ir] >> wd[ir];
			getline(ifs, line);
		}
	}
	else if(spd_loc==0 && pp.nlcc==0)
	{
		for(int ir=0; ir<pp.mesh; ++ir)
		{
			ifs >> pp.r[ir] >> vs[ir] >> vp[ir] >> vd[ir] 
				>> ws[ir] >> wp[ir] >> wd[ir] >> pp.vloc_at[ir];
			getline(ifs, line);
		}
	}
	else if(spd_loc>0 && pp.nlcc==1)
	{
		for(int ir=0; ir<pp.mesh; ++ir)
		{
			ifs >> pp.r[ir] >> vs[ir] >> vp[ir] >> vd[ir] 
				>> ws[ir] >> wp[ir] >> wd[ir] >> pp.rho_atc[ir];
			getline(ifs, line);
		}
	}
	else if(spd_loc==0 && pp.nlcc==1)
	{
		for(int ir=0; ir<pp.mesh; ++ir)
		{
			ifs >> pp.r[ir] >> vs[ir] >> vp[ir] >> vd[ir] 
				>> ws[ir] >> wp[ir] >> wd[ir] >> pp.vloc_at[ir] >> pp.rho_atc[ir];
			getline(ifs, line);
		}
	}
	// Hartree to Rydberg
	for(int ir=0; ir<pp.mesh; ++ir)
	{
		vs[ir] *= 2.0;
		vp[ir] *= 2.0;
		vd[ir] *= 2.0;
		pp.vloc_at[ir] *= 2.0;
	}


	// --------------------------------------
	// (2) check unit
	// --------------------------------------
	// calculate rab;
	// rab may not be accurate enough
	pp.rab[0] = pp.r[0];
	for(int ir=1; ir<pp.mesh-1; ++ir)
	{
		pp.rab[ir]=(pp.r[ir+1]-pp.r[ir-1])/2.0;
	}
	// check unit of vs, vp, vd
	double units = 0.0;
	double unitp = 0.0;
	double unitd = 0.0;
	for(int ir=1; ir<pp.mesh-1; ++ir)
	{
		double dr = (pp.r[ir+1]-pp.r[ir-1])/2.0;
		units += ws[ir] * ws[ir] * pp.r[ir] * pp.r[ir] * dr;
		unitp += wp[ir] * wp[ir] * pp.r[ir] * pp.r[ir] * dr;
		unitd += wd[ir] * wd[ir] * pp.r[ir] * pp.r[ir] * dr;
	}
	GlobalV::ofs_running << std::setw(15) << "WFC_UNIT" << std::setw(15) << units 
	<< std::setw(15) << unitp << std::setw(15) << unitd << std::endl; 


	// because only the rank=0 procesor read the pseudopotential
	// information, in order to make all the processors to stop
	// the job, we need to return the error information first.
	// we need to choose a threshold for the deviation of the
	// norm of pseudo atomic orbitals, I set 0.2
	// mohan 2013-06-28
	if( std::abs(units-1.0) > 0.2 && (iTB_s==1 || iref_s==1)) {return 3;}
	if( std::abs(unitp-1.0) > 0.2 && (iTB_p==1 || iref_p==1)) {return 3;}
	if( std::abs(unitd-1.0) > 0.2 && (iTB_d==1 || iref_d==1)) {return 3;}


	// calculate the phi*r*sqrt(4pi)
	pp.chi.create(pp.nchi,pp.mesh);
	for(int ir=0; ir<pp.mesh; ++ir)
	{
		int iwfc=0;
		if(iTB_s==1){pp.chi(iwfc,ir) = ws[ir]*pp.r[ir];++iwfc;}
		if(iTB_p==1){pp.chi(iwfc,ir) = wp[ir]*pp.r[ir];++iwfc;}
		if(iTB_d==1){pp.chi(iwfc,ir) = wd[ir]*pp.r[ir];++iwfc;}
	}
	// rho atom 
	for(int ir=0; ir<pp.mesh; ++ir)
	{
		for(int iwfc=0; iwfc<pp.nchi; ++iwfc)
		{
			pp.rho_at[ir] += pp.oc[iwfc]*pp.chi(iwfc,ir)*pp.chi(iwfc,ir);
		}
	}



	// --------------------------------------
	// (4) local pseudopotential 
	// --------------------------------------
	if(spd_loc==0) { for(int ir=0; ir<pp.mesh; ++ir)
	{
			// do nothing	
	}
	} else if(spd_loc==1) { for(int ir=0; ir<pp.mesh; ++ir) { pp.vloc_at[ir] = vs[ir]; 
}
	} else if(spd_loc==2) { for(int ir=0; ir<pp.mesh; ++ir) { pp.vloc_at[ir] = vp[ir];
}
	} else if(spd_loc==3) { for(int ir=0; ir<pp.mesh; ++ir) { pp.vloc_at[ir] = vd[ir];
}
	} else if(spd_loc==12) { for(int ir=0; ir<pp.mesh; ++ir) { pp.vloc_at[ir] = (vs[ir]+vp[ir])/2.0;
}
	} else if(spd_loc==13) { for(int ir=0; ir<pp.mesh; ++ir) { pp.vloc_at[ir] = (vs[ir]+vd[ir])/2.0;
}
	} else if(spd_loc==23) { for(int ir=0; ir<pp.mesh; ++ir) { pp.vloc_at[ir] = (vp[ir]+vd[ir])/2.0;
}
}


	// --------------------------------------
	// (5) setup nonlocal pseudopotentials
	// --------------------------------------
	// for non-local pseudopotentials.
	if(iref_d==1) { pp.lmax=2;
	} else if(iref_p==1) { pp.lmax=1;
	} else if(iref_s==1) { pp.lmax=0;
	} else
	{
		std::cout << "\n !!! READ THIS FIRST !!!" << std::endl;
		std::cout << " Could not decide which is the max angular momentum from .vwr pseudopotential file." << std::endl;
		std::cout << " No reference states in .vwr pseudopotential file." << std::endl;
		std::cout << " That's incorrect, please check the refenrece states in .vwr file.";
		std::cout << "\n !!! READ THIS FIRST !!!" << std::endl;
		return 3;	
	}
	// no projectors now
	pp.nbeta = 0;
	if(iref_s==1) { ++pp.nbeta; // add one s projector
}
	if(iref_p==1) { ++pp.nbeta; // add one p projector
}
	if(iref_d==1) { ++pp.nbeta; // add one p projector
}
	GlobalV::ofs_running << std::setw(15) << "NPROJ" << std::setw(15) << pp.nbeta << std::endl;
	this->nd = pp.nbeta;
	GlobalV::ofs_running << std::setw(15) << "N-Dij" << std::setw(15) << nd << std::endl;
	// calculate the angular momentum for each pp.betar
	pp.lll = std::vector<int>(pp.nbeta, 0); 
	int icount=0;
	if(iref_s==1) {pp.lll[icount]=0; ++icount;}// s projector
	if(iref_p==1) {pp.lll[icount]=1; ++icount;}// p projector
	if(iref_d==1) {pp.lll[icount]=2; ++icount;}// p projector
	for(int i=0; i<pp.nbeta; ++i)
	{
		GlobalV::ofs_running << " lll[" << i << "]=" << pp.lll[i] << std::endl;
	}
    // this->kbeta(pp.nbeta): number of mesh points for projector i (must be .le.mesh )
    this->kbeta = std::vector<int>(pp.nbeta, 0);
    pp.kkbeta = 0;
    for (int ib = 0; ib < pp.nbeta; ++ib)
    {
        this->kbeta[ib] = pp.mesh;
        pp.kkbeta = (this->kbeta[ib] > pp.kkbeta) ? this->kbeta[ib] : pp.kkbeta;
    }
    // nonlocal projector
    pp.betar.create(pp.nbeta,pp.mesh);
	// coefficients
	pp.dion.create(pp.nbeta,pp.nbeta);


	// --------------------------------------
	// (6) generate nonlocal pseudopotentials
	// --------------------------------------
	// tmp function to evaluate < pp.betar | delta_v | pp.betar>
	std::vector<double> func = std::vector<double>(pp.mesh, 0.0);
	// tmp value (vs, vp or vd)
	std::vector<double> vl = std::vector<double>(pp.mesh, 0.0);
	// tmp wave function (ws, wp or wd with r)
	std::vector<double> wlr = std::vector<double>(pp.mesh, 0.0);
	double rcut = 5.0/1.03;
	GlobalV::ofs_running << std::setw(15) << "RCUT_NL" << std::setw(15) << rcut << std::endl;
	for(int ib=0; ib<pp.nbeta; ++ib)
	{
		double coef = 0.0;
		const int lnow = pp.lll[ib];
		if(lnow==0) { for(int ir=0; ir<pp.mesh; ++ir){vl[ir]=vs[ir]; wlr[ir]=ws[ir]*pp.r[ir];}
		} else if(lnow==1) { for(int ir=0; ir<pp.mesh; ++ir){vl[ir]=vp[ir]; wlr[ir]=wp[ir]*pp.r[ir];}
		} else if(lnow==2) { for(int ir=0; ir<pp.mesh; ++ir){vl[ir]=vd[ir]; wlr[ir]=wd[ir]*pp.r[ir];}
}
		// for non-local projectors
		// note that < phi | dV | phi > integration must have 4pi,
		// this 4pi is also needed in < phi | phi > = 1 integration.
		// However, this phi has sqrt(sphi) already because I
		// found < phi | phi > = 1 directly.
		GlobalV::ofs_running << " Projector index = " << ib+1 << ", L = " << lnow << std::endl;
		for(int ir=2; ir<pp.mesh-1; ++ir)
		{
			// p nl
			pp.betar(ib,ir)=(vl[ir]-pp.vloc_at[ir])*wlr[ir];
			if(pp.r[ir]<rcut)
			{
				coef=coef+(vl[ir]-pp.vloc_at[ir])*wlr[ir]*wlr[ir]*(pp.r[ir+1]-pp.r[ir-1])/2.0;
			}
		}


// In pw they did this:
//		pp.dion(ib,ib)=1.0/coef;

        if(coef<0.0) { pp.dion(ib,ib) = -1.0;
}
        if(coef>=0.0) { pp.dion(ib,ib) = 1.0;
}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// suppose wave function have sqrt(4pi) already
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        coef=1.0/sqrt(std::abs(coef));
		GlobalV::ofs_running << std::setw(25) << "1/sqrt(<phi|deltaV|phi>)" << std::setw(15) << coef << std::endl;
		for(int ir=0; ir<pp.mesh; ++ir)
		{
			pp.betar(ib,ir) *= coef;
			// --------- FOR TEST ---------
			if(ib>2)
			{
//				pp.betar(ib,ir) *= 0.0; // for test, disable Non-local
			}
			// --------- FOR TEST ---------
		}
		
	}

	// print out the projector.
	/*
	GlobalV::ofs_running << " Nonlocal projector : " << std::endl;
	for(int ir=0; ir<mesh; ++ir)
	{
		if(r[ir]<rcut)
		{
			GlobalV::ofs_running << " " << std::setw(15) << r[ir]; 
			for(int ib=0; ib<pp.nbeta; ++ib)
			{
				GlobalV::ofs_running << std::setw(25) << pp.betar(ib,ir);
			}
			GlobalV::ofs_running << std::endl;
		}
	}
	*/
	GlobalV::ofs_running << " -------------------------------------------------" << std::endl;
	return 0;
}
