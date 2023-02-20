#include "esolver_fp.h"
#include "../module_base/global_variable.h"
#include "../module_base/element_elec_config.h"
#include "../module_hamilt_pw/hamilt_pwdft/global.h"
namespace ModuleESolver
{   ESolver_FP::ESolver_FP()
    {
        // pw_rho = new ModuleBase::PW_Basis();
        
        pw_rho = new ModulePW::PW_Basis_Big(); 
        GlobalC::rhopw = this->pw_rho; //Temporary
        //temporary, it will be removed
        GlobalC::bigpw = static_cast<ModulePW::PW_Basis_Big*>(pw_rho);
        GlobalC::bigpw->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);
    }
    ESolver_FP::~ESolver_FP()
    {
        delete pw_rho;
        delete this->pelec;
    }
    void ESolver_FP::Init(Input& inp, UnitCell& cell)
    {
        this->read_pseudo(cell, GlobalV::ofs_running);

#ifdef __MPI
            this->pw_rho->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
        if (this->classname == "ESolver_OF") this->pw_rho->setfullpw(inp.of_full_pw, inp.of_full_pw_dim);
        // Initalize the plane wave basis set
        if (inp.nx * inp.ny * inp.nz == 0)
            this->pw_rho->initgrids(cell.lat0, cell.latvec, inp.ecutrho);
	    else
            this->pw_rho->initgrids(cell.lat0, cell.latvec, inp.nx, inp.ny, inp.nz);
        
        this->pw_rho->initparameters(false, inp.ecutrho);
        this->pw_rho->setuptransform();
        this->pw_rho->collect_local_pw(); 
        this->pw_rho->collect_uniqgg();
        this->print_rhofft(inp, GlobalV::ofs_running);
        
    }

    void ESolver_FP::read_pseudo(UnitCell &cell, ofstream &ofs)
    {
        // read in non-local pseudopotential and ouput the projectors.
        ofs << "\n\n\n\n";
        ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        ofs << " |                                                                    |" << std::endl;
        ofs << " | Reading pseudopotentials files:                                    |" << std::endl;
        ofs << " | The pseudopotential file is in UPF format. The 'NC' indicates that |" << std::endl;
        ofs << " | the type of pseudopotential is 'norm conserving'. Functional of    |" << std::endl;
        ofs << " | exchange and correlation is decided by 4 given parameters in UPF   |" << std::endl;
        ofs << " | file.  We also read in the 'core correction' if there exists.      |" << std::endl;
        ofs << " | Also we can read the valence electrons number and the maximal      |" << std::endl;
        ofs << " | angular momentum used in this pseudopotential. We also read in the |" << std::endl;
        ofs << " | trail wave function, trail atomic density and local-pseudopotential|" << std::endl;
        ofs << " | on logrithmic grid. The non-local pseudopotential projector is also|" << std::endl;
        ofs << " | read in if there is any.                                           |" << std::endl;
        ofs << " |                                                                    |" << std::endl;
        ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        ofs << "\n\n\n\n";

        cell.read_cell_pseudopots(GlobalV::global_pseudo_dir, ofs);

        if(GlobalV::MY_RANK == 0 && GlobalV::out_element_info)
        {
            for(int it=0; it<cell.ntype; it++)
            {
                Atom* atom = &cell.atoms[it];
                std::stringstream ss;
                ss << GlobalV::global_out_dir << atom->label 
                << "/" << atom->label
                << ".NONLOCAL";
                std::ofstream ofs(ss.str().c_str());
                
                ofs << "<HEADER>" << std::endl;
                ofs << std::setw(10) << atom->label << "\t" << "label" << std::endl;
                ofs << std::setw(10) << atom->ncpp.pp_type << "\t" << "pseudopotential type" << std::endl;
                ofs << std::setw(10) << atom->ncpp.lmax << "\t" << "lmax" << std::endl;
                ofs << "</HEADER>" << std::endl;
                
                ofs << "\n<DIJ>" << std::endl;
                ofs << std::setw(10) << atom->ncpp.nbeta << "\t" << "nummber of projectors." << std::endl;
                for(int ib=0; ib<atom->ncpp.nbeta; ib++)
                {
                    for(int ib2=0; ib2<atom->ncpp.nbeta; ib2++)
                    {
                        ofs<<std::setw(10) << atom->ncpp.lll[ib] 
                        << " " << atom->ncpp.lll[ib2]
                        << " " << atom->ncpp.dion(ib,ib2)<<std::endl;
                    }
                }
                ofs << "</DIJ>" << std::endl;

                for(int i=0; i<atom->ncpp.nbeta; i++)
                {
                    ofs << "<PP_BETA>" << std::endl;
                    ofs << std::setw(10) << i << "\t" << "the index of projectors." <<std::endl;
                    ofs << std::setw(10) << atom->ncpp.lll[i] << "\t" << "the angular momentum." <<std::endl;

                    // mohan add
                    // only keep the nonzero part.
                    int cut_mesh = atom->ncpp.mesh; 
                    for(int j=atom->ncpp.mesh-1; j>=0; --j)
                    {
                        if( abs( atom->ncpp.betar(i,j) ) > 1.0e-10 )
                        {
                            cut_mesh = j; 
                            break;
                        }
                    }
                    if(cut_mesh %2 == 0) ++cut_mesh;

                    ofs << std::setw(10) << cut_mesh << "\t" << "the number of mesh points." << std::endl;

                    for(int j=0; j<cut_mesh; ++j)
                    {
                        ofs << std::setw(15) << atom->ncpp.r[j]
                            << std::setw(15) << atom->ncpp.betar(i, j)
                            << std::setw(15) << atom->ncpp.rab[j] << std::endl;
                    }
                    ofs << "</PP_BETA>" << std::endl;
                }

                ofs.close();
            }
        }

#ifdef __MPI
        cell.bcast_unitcell_pseudo2();
#endif

        for(int it=0; it<cell.ntype; it++)
        {
            if(cell.atoms[0].ncpp.xc_func !=cell.atoms[it].ncpp.xc_func)
            {
                GlobalV::ofs_warning << "\n type " << cell.atoms[0].label << " functional is " 
                    << cell.atoms[0].ncpp.xc_func;

                GlobalV::ofs_warning << "\n type " << cell.atoms[it].label << " functional is " 
                    << cell.atoms[it].ncpp.xc_func << std::endl;

                ModuleBase::WARNING_QUIT("setup_cell","All DFT functional must consistent.");
            }
        }

        cell.check_structure(GlobalV::MIN_DIST_COEF);

        // setup the total number of PAOs
        cell.cal_natomwfc(ofs);

        // setup GlobalV::NLOCAL
        cell.cal_nwfc(ofs);

        // Check whether the number of valence is minimum 
        if(GlobalV::MY_RANK==0)
        {
            int abtype = 0;
            for(int it=0; it<cell.ntype; it++)
            {
                if(ModuleBase::MinZval.find(cell.atoms[it].ncpp.psd) != ModuleBase::MinZval.end())
                {
                    if(cell.atoms[it].ncpp.zv > ModuleBase::MinZval.at(cell.atoms[it].ncpp.psd))
                    {
                        abtype += 1;
                        if(abtype == 1)
                        {
                            std::cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
                            ofs << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
                        }
                        std::cout<<" Warning: number valence electrons > " << ModuleBase::MinZval.at(cell.atoms[it].ncpp.psd);
                        std::cout<<" for " << cell.atoms[it].ncpp.psd << ": " << ModuleBase::EleConfig.at(cell.atoms[it].ncpp.psd) << std::endl;
                        ofs << " Warning: number valence electrons > " << ModuleBase::MinZval.at(cell.atoms[it].ncpp.psd);
                        ofs << " for " << cell.atoms[it].ncpp.psd << ": " << ModuleBase::EleConfig.at(cell.atoms[it].ncpp.psd) << std::endl;
                    }
                }
            }
            if(abtype>0)
            {
                std::cout<< " Please make sure the pseudopotential file is what you need"<<std::endl;
                std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"<<std::endl;
                ofs << " Please make sure the pseudopential file is what you need"<<std::endl;
                ofs << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
                ModuleBase::GlobalFunc::OUT(ofs,"");
            }
        }

        cell.cal_meshx();
    }

    void ESolver_FP::print_rhofft(Input&inp, ofstream &ofs)
    {
        std::cout << " UNIFORM GRID DIM     : " << GlobalC::rhopw->nx << " * " << GlobalC::rhopw->ny << " * " << GlobalC::rhopw->nz << std::endl;
        std::cout << " UNIFORM GRID DIM(BIG): " << GlobalC::bigpw->nbx << " * " << GlobalC::bigpw->nby << " * " << GlobalC::bigpw->nbz << std::endl;

        ofs << "\n\n\n\n";
	    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	    ofs << " |                                                                    |" << std::endl;
	    ofs << " | Setup plane waves of charge/potential:                             |" << std::endl;
	    ofs << " | Use the energy cutoff and the lattice vectors to generate the      |" << std::endl;
	    ofs << " | dimensions of FFT grid. The number of FFT grid on each processor   |" << std::endl;
	    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space is   |" << std::endl;
	    ofs << " | different for charege/potential and wave functions. We also set    |" << std::endl;
	    ofs << " | the 'sticks' for the parallel of FFT. The number of plane waves    |" << std::endl;
	    ofs << " | is 'npw' in each processor.                                        |" << std::endl;
	    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	    ofs << "\n\n\n\n";
	    ofs << "\n SETUP THE PLANE WAVE BASIS" << std::endl;
        double ecut = INPUT.ecutrho;
        if(inp.nx * inp.ny * inp.nz > 0)
        {
            ecut = this->pw_rho->gridecut_lat * this->pw_rho->tpiba2;
            ofs << "use input fft dimensions for wave functions." << std::endl;
            ofs << "calculate energy cutoff from nx, ny, nz:" << std::endl;

        }
        ModuleBase::GlobalFunc::OUT(ofs,"energy cutoff for charge/potential (unit:Ry)", ecut);
            
	    ModuleBase::GlobalFunc::OUT(ofs,"fft grid for charge/potential", this->pw_rho->nx,this->pw_rho->ny,this->pw_rho->nz);
	    ModuleBase::GlobalFunc::OUT(ofs,"fft grid division",GlobalC::bigpw->bx,GlobalC::bigpw->by,GlobalC::bigpw->bz);
	    ModuleBase::GlobalFunc::OUT(ofs,"big fft grid for charge/potential",GlobalC::bigpw->nbx,GlobalC::bigpw->nby,GlobalC::bigpw->nbz);
        ModuleBase::GlobalFunc::OUT(ofs,"nbxx",GlobalC::bigpw->nbxx);
	    ModuleBase::GlobalFunc::OUT(ofs,"nrxx",this->pw_rho->nrxx);

        ofs << "\n SETUP PLANE WAVES FOR CHARGE/POTENTIAL" << std::endl;
        ModuleBase::GlobalFunc::OUT(ofs,"number of plane waves",this->pw_rho->npwtot);
	    ModuleBase::GlobalFunc::OUT(ofs,"number of sticks", this->pw_rho->nstot);

        ofs << "\n PARALLEL PW FOR CHARGE/POTENTIAL" << std::endl;
        ofs <<" "<< std::setw(8)  << "PROC"<< std::setw(15) << "COLUMNS(POT)"<< std::setw(15) << "PW" << std::endl;
        for (int i = 0; i < GlobalV::NPROC_IN_POOL ; ++i)
        {
            ofs <<" "<<std::setw(8)<< i+1 << std::setw(15) << this->pw_rho->nst_per[i] << std::setw(15) << this->pw_rho->npw_per[i] << std::endl;
        }
        ofs << " --------------- sum -------------------" << std::endl;
        ofs << " " << std::setw(8)  << GlobalV::NPROC_IN_POOL << std::setw(15) << this->pw_rho->nstot << std::setw(15) << this->pw_rho->npwtot << std::endl;
        
        ModuleBase::GlobalFunc::OUT(ofs,"number of |g|", this->pw_rho->ngg);
        ModuleBase::GlobalFunc::OUT(ofs,"max |g|", this->pw_rho->gg_uniq[ this->pw_rho->ngg-1]);
	    ModuleBase::GlobalFunc::OUT(ofs,"min |g|", this->pw_rho->gg_uniq[0]);
    }
}