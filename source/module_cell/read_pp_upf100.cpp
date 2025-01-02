#include "read_pp.h"

//  read pseudopot_upf potential "upf" in the Unified
//  Pseudopot_upfpotential Format
int Pseudopot_upf::read_pseudo_upf(std::ifstream &ifs, Atom_pseudo& pp)
{
    std::string dummy;
    pp.has_so = false;
    this->q_with_l = false;
    this->mesh_changed = false;

    // addinfo_loop
    ifs.rdstate();

    while (ifs.good())
    {
        ifs >> dummy;
        if (dummy == "<PP_ADDINFO>")
        {
            pp.has_so = true;
        }
        else if (dummy == "<PP_QIJ_WITH_L>")
        {
            this->q_with_l = true;
        }

        if (pp.has_so && q_with_l)
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
            read_pseudo_header(ifs, pp);
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
        read_pseudo_mesh(ifs, pp);
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_MESH>");
    }

    // If  present, search for nlcc
    if (pp.nlcc)
    {
		ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NLCC>"); 
		//---------------------
		// call member function
		//---------------------
		read_pseudo_nlcc(ifs, pp);
		ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NLCC>");
	}

    if (!this->coulomb_potential)
    {
        // Search for Local potential
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_LOCAL>");
        //---------------------
        // call member function
        //---------------------
        read_pseudo_local(ifs, pp);
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_LOCAL>");
    }

	// Search for Nonlocal potential
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NONLOCAL>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_nl(ifs, pp);
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NONLOCAL>");

	// Search for atomic wavefunctions
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_PSWFC>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_pswfc(ifs, pp);
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_PSWFC>");

	// Search for atomic charge
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RHOATOM>");
	//---------------------
	// call member function
	//---------------------
	read_pseudo_rhoatom(ifs, pp);
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RHOATOM>");

	// Search for add_info
	if (pp.has_so)
	{
		ModuleBase::GlobalFunc::SCAN_BEGIN (ifs,"<PP_ADDINFO>");//added by zhengdy-soc
		//---------------------
		// call member function
		//---------------------
		read_pseudo_so(ifs, pp);
		ModuleBase::GlobalFunc::SCAN_END (ifs,"</PP_ADDINFO>");
	}

	ifs.clear();
	ifs.seekg(0);

	// return 0: read in sucessfully.
	return 0;
}// end subroutine read_pseudopot_upf


void Pseudopot_upf::read_pseudo_header(std::ifstream &ifs, Atom_pseudo& pp)
{
	ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.nv);// Version number
	ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.psd);// Element label

	// Type of pseudo : NC or US
	ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.pp_type);
	if(pp.pp_type=="US")
	{
		pp.tvanp = true;
        this->coulomb_potential = false;
    }
    else if (pp.pp_type == "NC")
    {
        pp.tvanp = false;
        this->coulomb_potential = false;
    }
    else if (pp.pp_type == "1/r")
    {
        pp.tvanp = false;
        this->coulomb_potential = true;
    }
    else
    {
        // A bug here!!! can't quit together.
        std::cout << " pp_type=" << pp.pp_type << std::endl;
        ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_header", "unknown pseudo type");
    }

    // If use nlcc
    std::string nlc;
    ModuleBase::GlobalFunc::READ_VALUE(ifs, nlc);

    if (nlc == "T")
    {
		pp.nlcc = true;
	}
	else
	{
		pp.nlcc = false;
	}

	// mohan modify 2009-12-15
	std::string junk;
	ifs >> junk >> junk >> junk >> junk;
	ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.xc_func);

    ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.zv);
	ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.etotps);

	ifs >> pp.ecutwfc >> pp.ecutrho;
	ifs.ignore(75, '\n');

	ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.lmax);
	ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.mesh);
    if (pp.mesh%2 == 0)
	{
		pp.mesh -= 1;
        this->mesh_changed = true;
	}

	ifs >> pp.nchi >> pp.nbeta ;
	ifs.ignore(75, '\n');
	ifs.ignore(75, '\n');

	pp.els = std::vector<std::string>(pp.nchi, "");
	pp.lchi = std::vector<int>(pp.nchi, 0);
	pp.oc = std::vector<double>(pp.nchi, 0.0);

	for(int i=0;i<pp.nchi;i++)
	{
		ifs >> pp.els[i] >> pp.lchi[i] >> pp.oc[i];
	}
	if (this->coulomb_potential)
	{
		pp.nbeta = 0;
        pp.lmax = 0;
        this->lloc = 0;
	}
	return;
}

void Pseudopot_upf::read_pseudo_mesh(std::ifstream &ifs, Atom_pseudo& pp)
{
	assert(pp.mesh>0);

	pp.r = std::vector<double>(pp.mesh, 0.0);
	pp.rab = std::vector<double>(pp.mesh, 0.0);

	int ir = 0;

	if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_R>", false) )
	{
		for (ir = 0;ir < pp.mesh;ir++)
		{
			ifs >> pp.r[ir];
		}
        this->skip_number(ifs, this->mesh_changed);
		ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_R>");
	}

	if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RAB>", false) )
	{
		for (ir = 0;ir < pp.mesh;ir++)
		{
			ifs >> pp.rab[ir];
		}
        this->skip_number(ifs, this->mesh_changed);
		ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RAB>");
	}
	return;
}

void Pseudopot_upf::read_pseudo_nlcc(std::ifstream &ifs, Atom_pseudo& pp)
{
	assert(pp.mesh>0);
	pp.rho_atc = std::vector<double>(pp.mesh, 0.0);
	for (int ir = 0;ir < pp.mesh;ir++)
	{
		ifs >> pp.rho_atc[ir];
	}
    this->skip_number(ifs, this->mesh_changed);
	return;
}

void Pseudopot_upf::read_pseudo_local(std::ifstream &ifs, Atom_pseudo& pp)
{
	assert(pp.mesh>0);
	pp.vloc_at = std::vector<double>(pp.mesh, 0.0);

	for (int ir = 0;ir < pp.mesh;ir++)
	{
		ifs >> pp.vloc_at[ir];
	}
    this->skip_number(ifs, this->mesh_changed);
	return;
}

void Pseudopot_upf::read_pseudo_nl(std::ifstream &ifs, Atom_pseudo& pp)
{
    //	int nb, mb, n, ir, idum, ldum, lp, i, ikk;
    int nb = 0;
    int mb = 0;
    int ir = 0;
    int idum = 0;
    int ldum = 0;

    if (pp.nbeta == 0)
    {
        this->nqf = 0;
        pp.nqlc = 0;
        pp.kkbeta = 0;
        return;
    }
    else
    {
        this->kbeta = std::vector<int>(pp.nbeta, 0);
        pp.lll = std::vector<int>(pp.nbeta, 0);
        pp.betar.create(pp.nbeta, pp.mesh);
        pp.dion.create(pp.nbeta, pp.nbeta);
        pp.kkbeta = 0;

        for (int i = 0; i < pp.nbeta; i++)
        {
            ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_BETA>", false);
            ifs >> idum;
            ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.lll[i]);// nl_1
            ModuleBase::GlobalFunc::READ_VALUE(ifs, this->kbeta[i]); // nl_2
            if (this->kbeta[i] > pp.mesh)
            {
                this->kbeta[i] = pp.mesh;
            }
            // number of mesh points for projectors

            for (ir = 0; ir < this->kbeta[i]; ir++)
            {
				ifs >> pp.betar(i, ir);// nl_3
			}
			ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_BETA>");
            pp.kkbeta = (this->kbeta[i] > pp.kkbeta) ? this->kbeta[i] : pp.kkbeta;
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
			pp.dion(mb, nb) = swap;// nl_5
			pp.dion(nb, mb) = swap;
        }
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_DIJ>");

		// QIJ
		if (pp.tvanp)
		{
            if (!ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_QIJ>", false))
            {
                ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_QIJ_WITH_L>", false);
            }
            // If nqf is not zero, Qij's inside rinner are computed using qfcoef's
            ModuleBase::GlobalFunc::READ_VALUE(ifs, this->nqf);
            pp.nqlc = 2 * pp.lmax + 1;
            this->rinner = std::vector<double>(pp.nqlc, 0.0);
            pp.qqq.create(pp.nbeta, pp.nbeta);
            if (q_with_l)
            {
                pp.qfuncl.create(2 * pp.lmax + 1, pp.nbeta * (pp.nbeta + 1) / 2, pp.mesh);
            }
            else
            {
                this->qfunc.create(pp.nbeta * (pp.nbeta + 1) / 2, pp.mesh);
            }

            if (nqf <= 0)
            {
                this->qfcoef.create(1, 1, 1, 1);
            }
            else
            {
                this->qfcoef.create(pp.nbeta, pp.nbeta, pp.nqlc, nqf);
                ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RINNER>", false);
                for (int i = 0; i < pp.nqlc; i++)
                {
                    ifs >> idum >> rinner[i];
                }
                ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RINNER>");
            }

            for (int nb = 0; nb < pp.nbeta; nb++)
            {
                int ln = pp.lll[nb];
                for (int mb = nb; mb < pp.nbeta; mb++)
                {
                    int lm = pp.lll[mb];
                    int nmb = mb * (mb + 1) / 2 + nb;
                    ifs >> idum >> idum >> ldum; // i  j  (l(j))
                    ifs.ignore(75, '\n');

                    if (ldum != lm)
                    {
                        ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_nl",
                                                 "inconsistent angular momentum for Q_ij");
                    }

                    ModuleBase::GlobalFunc::READ_VALUE(ifs, pp.qqq(nb, mb));
                    pp.qqq(mb, nb) = pp.qqq(nb, mb);

                    if (q_with_l)
                    {
                        for (int l = std::abs(ln - lm); l <= ln + lm; l += 2)
                        {
                            for (int ir = 0; ir < pp.mesh; ir++)
                            {
                                ifs >> pp.qfuncl(l, nmb, ir);
                            }
                        }
                    }
                    else
                    {
                        for (int ir = 0; ir < pp.mesh; ir++)
                        {
                            ifs >> qfunc(nmb, ir);
                        }
                    }
                    this->skip_number(ifs, this->mesh_changed);

                    if (this->nqf > 0)
                    {
                        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_QFCOEF>", false);
                        for (int k = 0; k < pp.nqlc; k++)
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

void Pseudopot_upf::read_pseudo_pswfc(std::ifstream &ifs, Atom_pseudo& pp)
{
	pp.chi.create(pp.nchi, pp.mesh);
	for (int i=0;i<pp.nchi;i++)
	{
		std::string OrbitalName;
		int BelongToL = 0;
		double occupation = 0.0;
		std::string dummy;
		ifs >> OrbitalName >> BelongToL >> occupation >> dummy;
		for (int ir = 0;ir < pp.mesh;ir++)
		{
			ifs >> pp.chi(i, ir);
		}
        this->skip_number(ifs, this->mesh_changed);
	}
	return;
}

void Pseudopot_upf::read_pseudo_rhoatom(std::ifstream &ifs, Atom_pseudo& pp)
{
	pp.rho_at = std::vector<double>(pp.mesh, 0.0);
	for (int ir = 0;ir < pp.mesh;ir++)
	{
		ifs >> pp.rho_at[ir];
	}
    this->skip_number(ifs, this->mesh_changed);
	return;
}

void Pseudopot_upf::read_pseudo_so(std::ifstream &ifs, Atom_pseudo& pp)
{
       //read soc info from upf, added by zhengdy-soc
       if(!pp.has_so) { return;}
       pp.nn = std::vector<int>(pp.nchi, 0);
       pp.jchi = std::vector<double>(pp.nchi, 0.0);
       pp.jjj = std::vector<double>(pp.nbeta, 0.0);
       //RELWFC
       for(int nw=0;nw< pp.nchi;nw++)
       {
             ifs >> pp.els[nw] >>pp.nn[nw] >> pp.lchi[nw] >> pp.jchi[nw] >> pp.oc[nw];
             if(pp.lchi[nw]-pp.jchi[nw]-0.5>1e-7 && pp.lchi[nw]-pp.jchi[nw]-0.5<1e-7)
             {
                  std::cout<<"Ignore ADDINFO section"<<std::endl;
                  pp.has_so = false;
             }
       }
       //RELBETA
       for(int nb = 0;nb < pp.nbeta;nb++)
       {
             ifs >> pp.lll[nb] >> pp.jjj[nb];
             if(pp.lll[nb]-pp.jjj[nb]-0.5>1e-7 && pp.lll[nb]-pp.jjj[nb]-0.5<1e-7)
             {
                  std::cout<<"Ignore ADDINFO section"<<std::endl;
                  pp.has_so = false;
             }
       }
       return;
}


void Pseudopot_upf::print_pseudo_upf(std::ofstream &ofs, Atom_pseudo& pp)
{
	ModuleBase::TITLE("Pseudopot_upf","print_pseudo_upf");
	ofs << " ==== read_pseudo_upf === " << std::endl;

	// print header
	ofs << " has_so: " << pp.has_so << std::endl;
	ofs << " Version number : " << pp.nv << std::endl;
	ofs << " Element label : " << pp.psd << std::endl;
	ofs << " pp_type: " << pp.pp_type << std::endl;
	ofs << " tvanp: " << pp.tvanp << std::endl;
	ofs << " nlcc: " << pp.nlcc << std::endl; 
	ofs << " dft: " << pp.xc_func << std::endl;
	ofs << " zp: " << pp.zv << std::endl;
	ofs << " etotps: " << pp.etotps << std::endl;
	ofs << " ecutwfc: " << pp.ecutwfc << std::endl;
	ofs << " ecutrho: " << pp.ecutrho << std::endl;
	ofs << " lmax: " << pp.lmax << std::endl;
	ofs << " mesh: " << pp.mesh << std::endl;
	ofs << " nwfc: " << pp.nchi << std::endl;
	ofs << " nbeta: " << pp.nbeta << std::endl;
	for(int i=0; i<pp.nchi; ++i)
	{
		ofs << " iw=" << i << " els=" << pp.els[i] << " lchi=" << pp.lchi[i] << " oc=" << pp.oc[i] << std::endl;
	}

	ofs << " End of pseudopot_upf." << std::endl;

	return;

}
