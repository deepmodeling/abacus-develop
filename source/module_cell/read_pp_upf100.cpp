#include "read_pp.h"

//  read pseudopot_upf potential "upf" in the Unified
//  Pseudopot_upfpotential Format
int Pseudopot_upf::read_pseudo_upf(std::ifstream &ifs)
{
    std::string dummy;
    this->has_so = false;
    this->q_with_l = false;

    // addinfo_loop
    ifs.rdstate();

    while (ifs.good())
    {
        ifs >> dummy;
        if (dummy == "<PP_ADDINFO>")
        {
            this->has_so = true;
        }
        else if (dummy == "<PP_QIJ_WITH_L>")
        {
            this->q_with_l = true;
        }

        if (has_so && q_with_l)
        {
            break;
        }
        ifs.rdstate();
    }

    // Search for Header
    // This version doesn't use the new routine SCAN_BEGIN
    // because this search must set extra flags for
    // compatibility with other pp format reading

    int ierr = 0;

    ifs.clear();
    ifs.seekg(0);
    ifs.rdstate();

    // header_loop:
    while (ifs.good())
    {
        ifs >> dummy;
        if (dummy == "<PP_HEADER>")
        {
            ierr = 1;
            //---------------------
            // call member function
            //---------------------
            read_pseudo_header(ifs);
            ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_HEADER>");
            break;
        }
    }

    if (ierr == 0)
    {
        // 2: something in pseudopotential file not match.
        return 2;
    }

    // Search for mesh information
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_MESH>"))
    {
        //---------------------
        // call member function
        //---------------------
        read_pseudo_mesh(ifs);
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_MESH>");
    }

    // If  present, search for nlcc
    if (this->nlcc)
    {
		ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NLCC>"); 
		//---------------------
		// call member function
		//---------------------
		read_pseudo_nlcc(ifs);
		ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NLCC>");
	}

    if (!this->coulomb_potential)
    {
        // Search for Local potential
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_LOCAL>");
        //---------------------
        // call member function
        //---------------------
        read_pseudo_local(ifs);
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_LOCAL>");
    }

	// Search for Nonlocal potential
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NONLOCAL>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_nl(ifs);
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NONLOCAL>");

	// Search for atomic wavefunctions
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_PSWFC>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_pswfc(ifs);
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_PSWFC>");

	// Search for atomic charge
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RHOATOM>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_rhoatom(ifs);
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RHOATOM>");

	// Search for add_info
	if (has_so)
	{
		ModuleBase::GlobalFunc::SCAN_BEGIN (ifs,"<PP_ADDINFO>");//added by zhengdy-soc
		//---------------------
		// call member function
		//---------------------
		read_pseudo_so (ifs);
		ModuleBase::GlobalFunc::SCAN_END (ifs,"</PP_ADDINFO>");
	}
    if (mesh%2 == 0)
	{
		mesh -= 1;
	}
	ifs.clear();
	ifs.seekg(0);

	// return 0: read in sucessfully.
	return 0;
}// end subroutine read_pseudopot_upf


void Pseudopot_upf::read_pseudo_header(std::ifstream &ifs)
{
	ModuleBase::GlobalFunc::READ_VALUE(ifs, this->nv);// Version number
	ModuleBase::GlobalFunc::READ_VALUE(ifs, this->psd);// Element label

	// Type of pseudo : NC or US
	ModuleBase::GlobalFunc::READ_VALUE(ifs, this->pp_type);
	if(pp_type=="US")
	{
		this->tvanp = true;
        this->coulomb_potential = false;
    }
    else if (pp_type == "NC")
    {
        this->tvanp = false;
        this->coulomb_potential = false;
    }
    else if (pp_type == "1/r")
    {
        this->tvanp = false;
        this->coulomb_potential = true;
    }
    else
    {
        // A bug here!!! can't quit together.
        std::cout << " pp_type=" << pp_type << std::endl;
        ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_header", "unknown pseudo type");
    }

    // If use nlcc
    std::string nlc;
    ModuleBase::GlobalFunc::READ_VALUE(ifs, nlc);

    if (nlc == "T")
    {
		this->nlcc = true;
	}
	else
	{
		this->nlcc = false;
	}

	// mohan modify 2009-12-15
	std::string junk;
	ifs >> junk >> junk >> junk >> junk;
	ModuleBase::GlobalFunc::READ_VALUE(ifs, xc_func);

    ModuleBase::GlobalFunc::READ_VALUE(ifs, this->zp);
	ModuleBase::GlobalFunc::READ_VALUE(ifs, this->etotps);

	ifs >> this->ecutwfc >> this->ecutrho;
	ifs.ignore(75, '\n');

	ModuleBase::GlobalFunc::READ_VALUE(ifs, this->lmax);
	ModuleBase::GlobalFunc::READ_VALUE(ifs, this->mesh);

	ifs >> this->nwfc >> this->nbeta ;
	ifs.ignore(75, '\n');
	ifs.ignore(75, '\n');

	delete[] els;
	delete[] lchi;
	delete[] oc;
	this->els = new std::string[nwfc];
	this->lchi = new int[nwfc];
	this->oc = new double[nwfc];

	ModuleBase::GlobalFunc::ZEROS(lchi, nwfc); // angular momentum of each orbital
	ModuleBase::GlobalFunc::ZEROS(oc, nwfc);//occupation of each orbital

	for(int i=0;i<nwfc;i++)
	{
		ifs >> els[i] >> this->lchi[i] >> this->oc[i];
	}
	if (this->coulomb_potential)
	{
		this->nbeta = 0;
        this->lmax = 0;
        this->lloc = 0;
	}
	return;
}

void Pseudopot_upf::read_pseudo_mesh(std::ifstream &ifs)
{
	assert(mesh>0);

	delete[] r;
	delete[] rab;
	this->r = new double[mesh];
	this->rab = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(r,mesh);
	ModuleBase::GlobalFunc::ZEROS(rab,mesh);

	int ir = 0;

	if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_R>", false) )
	{
		for (ir = 0;ir < mesh;ir++)
		{
			ifs >> this->r[ir];
		}
		ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_R>");
	}

	if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RAB>", false) )
	{
		for (ir = 0;ir < mesh;ir++)
		{
			ifs >> this->rab[ir];
		}
		ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RAB>");
	}
	return;
}

void Pseudopot_upf::read_pseudo_nlcc(std::ifstream &ifs)
{
	assert(mesh>0);
	delete[] rho_atc;
	this->rho_atc = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(rho_atc, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rho_atc[ir];
	}
	return;
}

void Pseudopot_upf::read_pseudo_local(std::ifstream &ifs)
{
	assert(mesh>0);
	delete[] vloc;
	this->vloc = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(vloc, mesh);

	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->vloc[ir];
	}
	
	return;
}

void Pseudopot_upf::read_pseudo_nl(std::ifstream &ifs)
{
    //	int nb, mb, n, ir, idum, ldum, lp, i, ikk;
    int nb = 0;
    int mb = 0;
    int ir = 0;
    int idum = 0;
    int ldum = 0;

    if (nbeta == 0)
    {
        this->nqf = 0;
        this->nqlc = 0;
        this->kkbeta = 0;
        return;
    }
    else
    {
        delete[] kbeta;
        delete[] lll;
        this->kbeta = new int[nbeta];
        this->lll = new int[nbeta];
        this->beta.create(nbeta, mesh);
        this->dion.create(nbeta, nbeta);
        kkbeta = 0;

        for (int i = 0; i < nbeta; i++)
        {
            ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_BETA>", false);
            ifs >> idum;
            ModuleBase::GlobalFunc::READ_VALUE(ifs, this->lll[i]);// nl_1
            ModuleBase::GlobalFunc::READ_VALUE(ifs, this->kbeta[i]); // nl_2
            // number of mesh points for projectors

            for (ir = 0; ir < kbeta[i]; ir++)
            {
				ifs >> this->beta(i, ir);// nl_3
			}
			ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_BETA>");
            kkbeta = (kbeta[i] > kkbeta) ? kbeta[i] : kkbeta;
        }

        // DIJ
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_DIJ>", false);
        ModuleBase::GlobalFunc::READ_VALUE(ifs, this->nd); // nl_4
        for (int i = 0; i < this->nd; i++)
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
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_DIJ>");

		// QIJ
		if (tvanp)
		{
            if (!ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_QIJ>", false))
            {
                ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_QIJ_WITH_L>", false);
            }
            // If nqf is not zero, Qij's inside rinner are computed using qfcoef's
            ModuleBase::GlobalFunc::READ_VALUE(ifs, this->nqf);
            this->nqlc = 2 * this->lmax + 1;
            delete[] rinner;
            this->rinner = new double[nqlc];
            qqq.create(nbeta, nbeta);
            if (q_with_l)
            {
                this->qfuncl.create(2 * lmax + 1, nbeta * (nbeta + 1) / 2, mesh);
            }
            else
            {
                this->qfunc.create(nbeta * (nbeta + 1) / 2, mesh);
            }

            if (nqf <= 0)
            {
                ModuleBase::GlobalFunc::ZEROS(rinner, nqlc);
                this->qfcoef.create(1, 1, 1, 1);
            }
            else
            {
                this->qfcoef.create(nbeta, nbeta, nqlc, nqf);
                ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RINNER>", false);
                for (int i = 0; i < nqlc; i++)
                {
                    ifs >> idum >> rinner[i];
                }
                ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RINNER>");
            }

            for (int nb = 0; nb < nbeta; nb++)
            {
                int ln = lll[nb];
                for (int mb = nb; mb < nbeta; mb++)
                {
                    int lm = lll[mb];
                    int nmb = mb * (mb + 1) / 2 + nb;
                    ifs >> idum >> idum >> ldum; // i  j  (l(j))
                    ifs.ignore(75, '\n');

                    if (ldum != lm)
                    {
                        ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_nl",
                                                 "inconsistent angular momentum for Q_ij");
                    }

                    ModuleBase::GlobalFunc::READ_VALUE(ifs, this->qqq(nb, mb));
                    this->qqq(mb, nb) = this->qqq(nb, mb);

                    if (q_with_l)
                    {
                        for (int l = std::abs(ln - lm); l <= ln + lm; l += 2)
                        {
                            for (int ir = 0; ir < mesh; ir++)
                            {
                                ifs >> qfuncl(l, nmb, ir);
                            }
                        }
                    }
                    else
                    {
                        for (int ir = 0; ir < mesh; ir++)
                        {
                            ifs >> qfunc(nmb, ir);
                        }
                    }

                    if (this->nqf > 0)
                    {
                        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_QFCOEF>", false);
                        for (int k = 0; k < nqlc; k++)
                        {
                            for (int l = 0; l < nqf; l++)
                            {
                                ifs >> qfcoef(nb, mb, k, l);
                                qfcoef(mb, nb, k, l) = qfcoef(nb, mb, k, l);
                            }
                        }
                        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_QFCOEF>");
                    }
                }
            }
            ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_QIJ>");
        }
        else // not tvanp
        {
        }
    }
    return;
}

void Pseudopot_upf::read_pseudo_pswfc(std::ifstream &ifs)
{
	this->chi.create(this->nwfc, this->mesh);
	for (int i=0;i<nwfc;i++)
	{
		std::string OrbitalName;
		int BelongToL = 0;
		double occupation = 0.0;
		std::string dummy;
		ifs >> OrbitalName >> BelongToL >> occupation >> dummy;
		for (int ir = 0;ir < mesh;ir++)
		{
			ifs >> this->chi(i, ir);
		}
	}
	return;
}

void Pseudopot_upf::read_pseudo_rhoatom(std::ifstream &ifs)
{
	delete[] rho_at;
	this->rho_at = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(rho_at, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rho_at[ir];
	}
	return;
}

void Pseudopot_upf::read_pseudo_so(std::ifstream &ifs)
{
       //read soc info from upf, added by zhengdy-soc
       if(!this->has_so) return;
       delete[] nn;
       delete[] jchi;
       delete[] jjj;
       this->nn = new int[nwfc];
       this->jchi = new double[nwfc];
       this->jjj = new double[nbeta];
       ModuleBase::GlobalFunc::ZEROS(nn,nwfc);
       ModuleBase::GlobalFunc::ZEROS(jchi,nwfc);
       ModuleBase::GlobalFunc::ZEROS(jjj,nbeta);
       //RELWFC
       for(int nw=0;nw< nwfc;nw++)
       {
             ifs >> this->els[nw] >>this->nn[nw] >> this->lchi[nw] >> this->jchi[nw] >> this->oc[nw];
             if(this->lchi[nw]-this->jchi[nw]-0.5>1e-7 && this->lchi[nw]-this->jchi[nw]-0.5<1e-7)
             {
                  std::cout<<"Ignore ADDINFO section"<<std::endl;
                  this->has_so = 0;
             }
       }
       //RELBETA
       for(int nb = 0;nb < nbeta;nb++)
       {
             ifs >> this->lll[nb] >> this->jjj[nb];
             if(this->lll[nb]-this->jjj[nb]-0.5>1e-7 && this->lll[nb]-this->jjj[nb]-0.5<1e-7)
             {
                  std::cout<<"Ignore ADDINFO section"<<std::endl;
                  this->has_so = 0;
             }
       }
       return;
}


void Pseudopot_upf::print_pseudo_upf(std::ofstream &ofs)
{
	ModuleBase::TITLE("Pseudopot_upf","print_pseudo_upf");
	ofs << " ==== read_pseudo_upf === " << std::endl;

	// print header
	ofs << " has_so: " << has_so << std::endl;
	ofs << " Version number : " << nv << std::endl;
	ofs << " Element label : " << psd << std::endl;
	ofs << " pp_type: " << pp_type << std::endl;
	ofs << " tvanp: " << tvanp << std::endl;
	ofs << " nlcc: " << nlcc << std::endl; 
	ofs << " dft: " << xc_func << std::endl;
	ofs << " zp: " << zp << std::endl;
	ofs << " etotps: " << etotps << std::endl;
	ofs << " ecutwfc: " << ecutwfc << std::endl;
	ofs << " ecutrho: " << ecutrho << std::endl;
	ofs << " lmax: " << lmax << std::endl;
	ofs << " mesh: " << mesh << std::endl;
	ofs << " nwfc: " << nwfc << std::endl;
	ofs << " nbeta: " << nbeta << std::endl;
	for(int i=0; i<nwfc; ++i)
	{
		ofs << " iw=" << i << " els=" << els[i] << " lchi=" << lchi[i] << " oc=" << oc[i] << std::endl;
	}

	ofs << " End of pseudopot_upf." << std::endl;

	return;

}
