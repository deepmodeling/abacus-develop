#include "read_pseudo.h"

namespace elecstate
{
    void cal_nwfc(std::ofstream& log,UnitCell& ucell,Atom* atoms) 
    {
        ModuleBase::TITLE("UnitCell", "cal_nwfc");
        const int ntype = ucell.ntype;
        const int nat   = ucell.nat;
        assert(ntype > 0);
        assert(nat > 0);

        //===========================
        // (1) set iw2l, iw2n, iw2m
        //===========================
        for (int it = 0; it < ntype; it++) 
        {
            ucell.atoms[it].set_index();
        }
        //===========================
        // (2) set namax and nwmax
        //===========================
        ucell.namax = 0;
        ucell.nwmax = 0;
        for (int it = 0; it < ntype; it++) {
            ucell.namax = std::max(atoms[it].na, ucell.namax);
            ucell.nwmax = std::max(atoms[it].nw, ucell.nwmax);
        }

        //===========================
        // (3) set nwfc and stapos_wf
        //===========================
        int nlocal_tmp = 0;
        for (int it = 0; it < ntype; it++) {
            atoms[it].stapos_wf = nlocal_tmp;
            const int nlocal_it = atoms[it].nw * atoms[it].na;
            if (PARAM.inp.nspin != 4) {
                nlocal_tmp += nlocal_it;
            } else {
                nlocal_tmp += nlocal_it * 2; // zhengdy-soc
            }
        }

//        log << " " << std::setw(40) << "NLOCAL"
  //          << " = " << nlocal_tmp << std::endl;
        //========================================================
        // (4) set index for itia2iat, itiaiw2iwt
        //========================================================

        // mohan add 2010-09-26
        assert(nlocal_tmp > 0);
        assert(nlocal_tmp == PARAM.globalv.nlocal);
        delete[] ucell.iwt2iat;
        delete[] ucell.iwt2iw;
        ucell.iwt2iat = new int[nlocal_tmp];
        ucell.iwt2iw = new int[nlocal_tmp];

        ucell.itia2iat.create(ntype, ucell.namax);
        ucell.set_iat2iwt(PARAM.globalv.npol);
        int iat = 0;
        int iwt = 0;
        for (int it = 0; it < ntype; it++) {
            for (int ia = 0; ia < atoms[it].na; ia++) {
                ucell.itia2iat(it, ia) = iat;
                for (int iw = 0; iw < atoms[it].nw * PARAM.globalv.npol; iw++) {
                    ucell.iwt2iat[iwt] = iat;
                    ucell.iwt2iw[iwt] = iw;
                    ++iwt;
                }
                ++iat;
            }
        }

        //========================
        // (5) set lmax and nmax
        //========================
        ucell.lmax = 0;
        ucell.nmax = 0;
        ucell.nmax_total = 0;
        for (int it = 0; it < ntype; it++) {
            ucell.lmax = std::max(ucell.lmax, atoms[it].nwl);
            for (int l = 0; l < atoms[it].nwl + 1; l++) {
                ucell.nmax = std::max(ucell.nmax, atoms[it].l_nchi[l]);
            }

            int nchi = 0;
            for (int l = 0; l < atoms[it].nwl + 1; l++) {
                nchi += atoms[it].l_nchi[l];
            }
            ucell.nmax_total = std::max(ucell.nmax_total, nchi);
        }

        //=======================
        // (6) set lmax_ppwf
        //=======================
        ucell.lmax_ppwf = 0;
        for (int it = 0; it < ntype; it++) {
            for (int ic = 0; ic < atoms[it].ncpp.nchi; ic++) {
                if (ucell.lmax_ppwf < atoms[it].ncpp.lchi[ic]) {
                    ucell.lmax_ppwf = atoms[it].ncpp.lchi[ic];
                }
            }
        }
        //=====================
        // Use localized basis
        //=====================
        if ((PARAM.inp.basis_type == "lcao") || (PARAM.inp.basis_type == "lcao_in_pw")
            || ((PARAM.inp.basis_type == "pw") && (PARAM.inp.init_wfc.substr(0, 3) == "nao")
                && (PARAM.inp.esolver_type == "ksdft"))) // xiaohui add 2013-09-02
        {
            ModuleBase::GlobalFunc::AUTO_SET("NBANDS", PARAM.inp.nbands);
        } else // plane wave basis
        {
            // if(winput::after_iter && winput::sph_proj)
            //{
            //	if(PARAM.inp.nbands < PARAM.globalv.nlocal)
            //	{
            //		ModuleBase::WARNING_QUIT("cal_nwfc","NBANDS must > PARAM.globalv.nlocal
            //!");
            //	}
            // }
        }

        return;
    }

    void cal_meshx(int& meshx,const Atom* atoms, const int ntype) 
    {
        meshx = 0;
        for (int it = 0; it < ntype; it++) {
            const int mesh = atoms[it].ncpp.msh;
            if (mesh > meshx) 
            {
                meshx = mesh;
            }
        }
    }


    void cal_natomwfc(std::ofstream& log,int& natomwfc,const int ntype,const Atom* atoms) 
    {
        natomwfc = 0;
		for (int it = 0; it < ntype; it++) 
		{
			//============================
			// Use pseudo-atomic orbitals
			//============================
			int tmp = 0;
			for (int l = 0; l < atoms[it].ncpp.nchi; l++) 
			{
				if (atoms[it].ncpp.oc[l] >= 0) 
				{
					if (PARAM.inp.nspin == 4) 
					{
						if (atoms[it].ncpp.has_so) 
						{
							tmp += 2 * atoms[it].ncpp.lchi[l];
							if (fabs(atoms[it].ncpp.jchi[l] - atoms[it].ncpp.lchi[l] - 0.5)< 1e-6) 
							{
								tmp += 2;
							}
						} else 
						{
							tmp += 2 * (2 * atoms[it].ncpp.lchi[l] + 1);
						}
					} else 
					{
						tmp += 2 * atoms[it].ncpp.lchi[l] + 1;
					}
				}
			}
			natomwfc += tmp * atoms[it].na;
		}
		ModuleBase::GlobalFunc::OUT(log, "Number of pseudo atomic orbitals", natomwfc);
		return;
    }
}
