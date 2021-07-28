#include "read_pp.h"

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
			//---------------------
			// call member function
			//---------------------
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
		//---------------------
		// call member function
		//---------------------
		read_pseudo_mesh(ifs);
		SCAN_END(ifs, "</PP_MESH>");
	}

	// If  present, search for nlcc
	if (this->nlcc)
	{
		SCAN_BEGIN(ifs, "<PP_NLCC>"); 
		//---------------------
		// call member function
		//---------------------
		read_pseudo_nlcc(ifs);
		SCAN_END(ifs, "</PP_NLCC>");
	}

	// Search for Local potential
	SCAN_BEGIN(ifs, "<PP_LOCAL>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_local(ifs);
	SCAN_END(ifs, "</PP_LOCAL>");

	// Search for Nonlocal potential
	SCAN_BEGIN(ifs, "<PP_NONLOCAL>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_nl(ifs);
	SCAN_END(ifs, "</PP_NONLOCAL>");

	// Search for atomic wavefunctions
	SCAN_BEGIN(ifs, "<PP_PSWFC>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_pswfc(ifs);
	SCAN_END(ifs, "</PP_PSWFC>");

	// Search for atomic charge
	SCAN_BEGIN(ifs, "<PP_RHOATOM>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_rhoatom(ifs);
	SCAN_END(ifs, "</PP_RHOATOM>");

	// Search for add_info
	if (has_so)
	{
		SCAN_BEGIN (ifs,"<PP_ADDINFO>");//added by zhengdy-soc
		//---------------------
		// call member function
		//---------------------
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
	ifs >> dft[0] >> dft[1] >> dft[2] >> dft[3];
	
	//dft[i](i=0-3) gives the four components of xc functional:
	//local X, local C, semilocal X, semilocal C
	//dft_tot is the name of the combination
	string dft_tot;
	READ_VALUE(ifs, dft_tot);
	
	// dft functional enforced to modify
	// mohan add 2010-07-15
	if(GlobalV::DFT_FUNCTIONAL!="none")
	{
		/*xiaohui modify 2015-03-24
		dft[0] = GlobalV::DFT_FUNCTIONAL;
		dft[1] = GlobalV::DFT_FUNCTIONAL;
		dft[2] = GlobalV::DFT_FUNCTIONAL;
		dft[3] = GlobalV::DFT_FUNCTIONAL;
		xiaohui modify 2015-03-24*/

		//xiaohui add 2015-03-23
		
		string dft_functional;
		if(dft[1] == "PZ")
		{
			dft_functional = "lda";
		}
		else if(dft[1] == "PBE")
		{
			dft_functional = "pbe";
		}
		else if(dft[1] == "SCAN")
		{
			dft_functional = "scan";
		}
		
		if(dft_tot != GlobalV::DFT_FUNCTIONAL)
		{
			functional_error = 1;

			cout << " dft_functional readin is: " << GlobalV::DFT_FUNCTIONAL << endl;
			cout << " dft_functional in pseudopot file is: " << dft_tot << endl;
			GlobalV::ofs_warning << " dft_functional readin is: " << GlobalV::DFT_FUNCTIONAL << endl;
			GlobalV::ofs_warning << " dft_functional in pseudopot file is: " << dft_tot << endl;
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
	
	return;
}

void Pseudopot_upf::read_pseudo_nl(ifstream &ifs)
{
//	int nb, mb, n, ir, idum, ldum, lp, i, ikk;
	int nb, mb, ir, idum;

	if (nbeta == 0)
	{
		delete[] kkbeta;
		delete[] lll;
		this->kkbeta = new int[1];
		this->lll = new int[1];
		this->beta.create(1, 1);
		this->dion.create(1, 1);
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
		}
		else // not tvanp
		{
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

	ofs << " End of pseudopot_upf." << endl;

	return;

}
