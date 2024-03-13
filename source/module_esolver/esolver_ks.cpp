#include "esolver_ks.h"

#include <time.h>
#ifdef __MPI
#include <mpi.h>
#else
#include <chrono>
#endif

#include <iostream>

#include "module_io/print_info.h"
#include "module_base/timer.h"
#include "module_io/input.h"
#include "module_io/json_output/init_info.h"


//--------------Temporary----------------
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//---------------------------------------

#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#include "module_base/parallel_common.h"
#endif

namespace ModuleESolver
{

    template<typename T, typename Device>
    ESolver_KS<T, Device>::ESolver_KS()
    {
        classname = "ESolver_KS";
        basisname = "PLEASE ADD BASISNAME FOR CURRENT ESOLVER.";
        scf_thr = GlobalV::SCF_THR;
        drho = 0.0;
        maxniter = GlobalV::SCF_NMAX;
        niter = maxniter;
        out_freq_elec = GlobalV::OUT_FREQ_ELEC;

        // pw_rho = new ModuleBase::PW_Basis();
        //temporary, it will be removed
        pw_wfc = new ModulePW::PW_Basis_K_Big(GlobalV::device_flag, GlobalV::precision_flag);
        ModulePW::PW_Basis_K_Big* tmp = static_cast<ModulePW::PW_Basis_K_Big*>(pw_wfc);
        tmp->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);

        ///----------------------------------------------------------
        /// charge mixing
        ///----------------------------------------------------------
        p_chgmix = new Charge_Mixing();
        p_chgmix->set_rhopw(this->pw_rho, this->pw_rhod);

        ///----------------------------------------------------------
        /// wavefunc
        ///----------------------------------------------------------
        this->wf.init_wfc = INPUT.init_wfc;
        this->wf.mem_saver = INPUT.mem_saver;
        this->wf.out_wfc_pw = INPUT.out_wfc_pw;
        this->wf.out_wfc_r = INPUT.out_wfc_r;
    }

    template<typename T, typename Device>
    ESolver_KS<T, Device>::~ESolver_KS()
    {
        delete this->psi;
        delete this->pw_wfc;
        delete this->p_hamilt;
        delete this->phsol;
        delete this->p_chgmix;
    }

    template<typename T, typename Device>
    void ESolver_KS<T, Device>::Init(Input& inp, UnitCell& ucell)
    {
        ESolver_FP::Init(inp,ucell);

        //------------------Charge Mixing------------------
        p_chgmix->set_mixing(GlobalV::MIXING_MODE,
                             GlobalV::MIXING_BETA,
                             GlobalV::MIXING_NDIM,
                             GlobalV::MIXING_GG0,
                             GlobalV::MIXING_TAU,
                             GlobalV::MIXING_BETA_MAG,
                             GlobalV::MIXING_GG0_MAG,
                             GlobalV::MIXING_GG0_MIN,
                             GlobalV::MIXING_ANGLE,
                             GlobalV::MIXING_DMR);
        
#ifdef USE_PAW
        if(GlobalV::use_paw)
        {
            int * atom_type;
            double ** atom_coord;
            std::vector<std::string> filename_list;

            atom_type = new int [GlobalC::ucell.nat];
            atom_coord = new double * [GlobalC::ucell.nat];
            filename_list.resize(GlobalC::ucell.ntype);

            for(int ia = 0; ia < GlobalC::ucell.nat; ia ++)
            {
                atom_coord[ia] = new double [3];
            }

            int iat = 0;
            for(int it = 0; it < GlobalC::ucell.ntype; it ++)
            {
                for(int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia ++)
                {
                    atom_type[iat] = it;
                    atom_coord[iat][0] = GlobalC::ucell.atoms[it].taud[ia].x;
                    atom_coord[iat][1] = GlobalC::ucell.atoms[it].taud[ia].y;
                    atom_coord[iat][2] = GlobalC::ucell.atoms[it].taud[ia].z;
                    iat ++;
                }
            }            

            if(GlobalV::MY_RANK == 0)
            {
                std::ifstream ifa(GlobalV::stru_file.c_str(), std::ios::in);
                if (!ifa)
                {
                    ModuleBase::WARNING_QUIT("set_libpaw_files", "can not open stru file");
                }

                std::string line;
                while(!ifa.eof())
                {
                    getline(ifa,line);
                    if (line.find("PAW_FILES") != std::string::npos) break;
                }

                for(int it = 0; it < GlobalC::ucell.ntype; it++)
                {
                    ifa >> filename_list[it];
                }
            }
#ifdef __MPI
            for(int it = 0; it < GlobalC::ucell.ntype; it++)
            {
                Parallel_Common::bcast_string(filename_list[it]);
            }
#endif
            
            GlobalC::paw_cell.init_paw_cell(INPUT.ecutwfc, INPUT.cell_factor,
                GlobalC::ucell.omega,GlobalC::ucell.nat,GlobalC::ucell.ntype,
                atom_type,(const double **) atom_coord,
                filename_list);

            for(int iat = 0; iat < GlobalC::ucell.nat; iat ++)
            {
                delete [] atom_coord[iat];
            }
            delete [] atom_coord;
            delete [] atom_type;
        }
#endif
        ucell.cal_nelec(GlobalV::nelec);

        /* it has been established that that
         xc_func is same for all elements, therefore
         only the first one if used*/
        if(GlobalV::use_paw)
        {
            XC_Functional::set_xc_type(GlobalV::DFT_FUNCTIONAL);
        }
        else
        {
            XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
        }
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

        // symmetry analysis should be performed every time the cell is changed
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
        }

        // Setup the k points according to symmetry.
        this->kv.set(ucell.symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

        // print information
        // mohan add 2021-01-30
        Print_Info::setup_parameters(ucell, this->kv);

        //if(GlobalV::BASIS_TYPE=="pw" || GlobalV::CALCULATION=="get_wf")
        //{
            //Envelope function is calculated as lcao_in_pw
            //new plane wave basis
#ifdef __MPI
        this->pw_wfc->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
        this->pw_wfc->initgrids(inp.ref_cell_factor * ucell.lat0,
                                ucell.latvec,
                                this->pw_rho->nx,
                                this->pw_rho->ny,
                                this->pw_rho->nz);
        this->pw_wfc->initparameters(false, inp.ecutwfc, this->kv.nks, this->kv.kvec_d.data());
#ifdef __MPI
        if (INPUT.pw_seed > 0)
            MPI_Allreduce(MPI_IN_PLACE, &this->pw_wfc->ggecut, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            // qianrui add 2021-8-13 to make different kpar parameters can get the same results
#endif
        this->pw_wfc->ft.fft_mode = inp.fft_mode;
        this->pw_wfc->setuptransform();
        for (int ik = 0; ik < this->kv.nks; ++ik)
            this->kv.ngk[ik] = this->pw_wfc->npwk[ik];
        this->pw_wfc->collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);
        this->print_wfcfft(inp, GlobalV::ofs_running);
        //}
        // initialize the real-space uniform grid for FFT and parallel
        // distribution of plane waves
        GlobalC::Pgrid.init(this->pw_rhod->nx,
                            this->pw_rhod->ny,
                            this->pw_rhod->nz,
                            this->pw_rhod->nplane,
                            this->pw_rhod->nrxx,
                            pw_big->nbz,
                            pw_big->bz); // mohan add 2010-07-22, update 2011-05-04

        // Calculate Structure factor
        this->sf.setup_structure_factor(&GlobalC::ucell, this->pw_rhod);

        // Initialize charge extrapolation
        CE.Init_CE(GlobalC::ucell.nat);

#ifdef USE_PAW
        if(GlobalV::use_paw)
        {
            GlobalC::paw_cell.set_libpaw_ecut(INPUT.ecutwfc/2.0,INPUT.ecutwfc/2.0); //in Hartree
            GlobalC::paw_cell.set_libpaw_fft(this->pw_wfc->nx,this->pw_wfc->ny,this->pw_wfc->nz,
                                            this->pw_wfc->nx,this->pw_wfc->ny,this->pw_wfc->nz,
                                            this->pw_wfc->startz,this->pw_wfc->numz);
#ifdef __MPI
            if(GlobalV::RANK_IN_POOL == 0) GlobalC::paw_cell.prepare_paw();
#else
            GlobalC::paw_cell.prepare_paw();
#endif
            GlobalC::paw_cell.set_sij();

            GlobalC::paw_cell.set_eigts(
                this->pw_wfc->nx,this->pw_wfc->ny,this->pw_wfc->nz,
                this->sf.eigts1.c,this->sf.eigts2.c,this->sf.eigts3.c);

            std::vector<std::vector<double>> rhoijp;
            std::vector<std::vector<int>> rhoijselect;
            std::vector<int> nrhoijsel;
#ifdef __MPI
            if(GlobalV::RANK_IN_POOL == 0)
            {
                GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

                for(int iat = 0; iat < GlobalC::ucell.nat; iat ++)
                {
                    GlobalC::paw_cell.set_rhoij(iat,nrhoijsel[iat],rhoijselect[iat].size(),rhoijselect[iat].data(),rhoijp[iat].data());
                }  
            }
#else
            GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

            for(int iat = 0; iat < GlobalC::ucell.nat; iat ++)
            {
                GlobalC::paw_cell.set_rhoij(iat,nrhoijsel[iat],rhoijselect[iat].size(),rhoijselect[iat].data(),rhoijp[iat].data());
            }
#endif
        }
#endif
    }

    template<typename T, typename Device>
    void ESolver_KS<T, Device>::init_after_vc(Input& inp, UnitCell& ucell)
    {
        ModuleBase::TITLE("ESolver_KS", "init_after_vc");

        ESolver_FP::init_after_vc(inp, ucell);

        if (GlobalV::md_prec_level == 2)
        {
            // initialize the real-space uniform grid for FFT and parallel
            // distribution of plane waves
            GlobalC::Pgrid.init(this->pw_rhod->nx,
                                this->pw_rhod->ny,
                                this->pw_rhod->nz,
                                this->pw_rhod->nplane,
                                this->pw_rhod->nrxx,
                                pw_big->nbz,
                                pw_big->bz); // mohan add 2010-07-22, update 2011-05-04

            // Calculate Structure factor
            this->sf.setup_structure_factor(&ucell, this->pw_rhod);
        }
    }

    template<typename T, typename Device>
    void ESolver_KS<T, Device>::hamilt2density(const int istep, const int iter, const double ethr)
    {
        ModuleBase::timer::tick(this->classname, "hamilt2density");
        //Temporarily, before HSolver is constructed, it should be overrided by
        //LCAO, PW, SDFT and TDDFT.
        //After HSolver is constructed, LCAO, PW, SDFT should delete their own
        //hamilt2density() and use:
        //this->phsol->solve(this->phamilt, this->pes, this->wf, ETHR);
        ModuleBase::timer::tick(this->classname, "hamilt2density");
    }

    template<typename T, typename Device>
    void ESolver_KS<T, Device>::print_wfcfft(Input& inp, std::ofstream &ofs)
    {
        ofs << "\n\n\n\n";
	    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	    ofs << " |                                                                    |" << std::endl;
	    ofs << " | Setup plane waves of wave functions:                               |" << std::endl;
	    ofs << " | Use the energy cutoff and the lattice vectors to generate the      |" << std::endl;
	    ofs << " | dimensions of FFT grid. The number of FFT grid on each processor   |" << std::endl;
	    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space is   |" << std::endl;
	    ofs << " | different for charege/potential and wave functions. We also set    |" << std::endl;
	    ofs << " | the 'sticks' for the parallel of FFT. The number of plane wave of  |" << std::endl;
	    ofs << " | each k-point is 'npwk[ik]' in each processor                       |" << std::endl;
	    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	    ofs << "\n\n\n\n";
        ofs << "\n SETUP PLANE WAVES FOR WAVE FUNCTIONS" << std::endl;
        double ecut = INPUT.ecutwfc;
        if(std::abs(ecut-this->pw_wfc->gk_ecut * this->pw_wfc->tpiba2) > 1e-6)
        {
            ecut = this->pw_wfc->gk_ecut * this->pw_wfc->tpiba2;
            ofs<<"Energy cutoff for wavefunc is incompatible with nx, ny, nz and it will be reduced!"<<std::endl;

        }  
        ModuleBase::GlobalFunc::OUT(ofs,"energy cutoff for wavefunc (unit:Ry)", ecut);
        ModuleBase::GlobalFunc::OUT(ofs,"fft grid for wave functions", this->pw_wfc->nx,this->pw_wfc->ny,this->pw_wfc->nz);
        ModuleBase::GlobalFunc::OUT(ofs,"number of plane waves",this->pw_wfc->npwtot);
	    ModuleBase::GlobalFunc::OUT(ofs,"number of sticks", this->pw_wfc->nstot);

        ofs << "\n PARALLEL PW FOR WAVE FUNCTIONS" << std::endl;
        ofs <<" "<< std::setw(8)  << "PROC"<< std::setw(15) << "COLUMNS(POT)"<< std::setw(15) << "PW" << std::endl;
        for (int i = 0; i < GlobalV::NPROC_IN_POOL ; ++i)
        {
            ofs <<" "<<std::setw(8)<< i+1 << std::setw(15) << this->pw_wfc->nst_per[i] << std::setw(15) << this->pw_wfc->npw_per[i] << std::endl;
        }
        ofs << " --------------- sum -------------------" << std::endl;
        ofs << " " << std::setw(8)  << GlobalV::NPROC_IN_POOL << std::setw(15) << this->pw_wfc->nstot << std::setw(15) << this->pw_wfc->npwtot << std::endl;
        ModuleBase::GlobalFunc::DONE(ofs, "INIT PLANEWAVE");
    }

    template<typename T, typename Device>
    void ESolver_KS<T, Device>::Run(const int istep, UnitCell& ucell)
    {
        if (!(GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md"
            || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax"))
        {
            this->othercalculation(istep);
        }
        else
        {
            ModuleBase::timer::tick(this->classname, "Run");

            this->beforescf(istep); //Something else to do before the iter loop
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT SCF");
            if(this->maxniter > 0)  this->printhead(); //print the headline on the screen.

            bool firstscf = true;
            this->conv_elec = false;
            this->niter = this->maxniter;
            for (int iter = 1; iter <= this->maxniter; ++iter)
            {
                writehead(GlobalV::ofs_running, istep, iter);
#ifdef __MPI
                auto iterstart = MPI_Wtime();
#else
                auto iterstart = std::chrono::system_clock::now();
#endif
                double diag_ethr = this->phsol->set_diagethr(istep, iter, drho);
                eachiterinit(istep, iter);
                this->hamilt2density(istep, iter, diag_ethr);

                //<Temporary> It may be changed when more clever parallel algorithm is put forward.
                //When parallel algorithm for bands are adopted. Density will only be treated in the first group.
                //(Different ranks should have abtained the same, but small differences always exist in practice.)
                //Maybe in the future, density and wavefunctions should use different parallel algorithms, in which 
                //they do not occupy all processors, for example wavefunctions uses 20 processors while density uses 10.
                if(GlobalV::MY_STOGROUP == 0)
                {
                    // double drho = this->estate.caldr2(); 
                    // EState should be used after it is constructed.

                    drho = p_chgmix->get_drho(pelec->charge, GlobalV::nelec);
                    double hsolver_error = 0.0;
                    if (firstscf)
                    {
                        firstscf = false;
                        hsolver_error = this->phsol->cal_hsolerror();
                        // The error of HSolver is larger than drho, so a more precise HSolver should be excuconv_elected.
                        if (hsolver_error > drho)
                        {
                            diag_ethr = this->phsol->reset_diagethr(GlobalV::ofs_running, hsolver_error, drho);
                            this->hamilt2density(istep, iter, diag_ethr);
                            drho = p_chgmix->get_drho(pelec->charge, GlobalV::nelec);
                            hsolver_error = this->phsol->cal_hsolerror();
                        }
                    }
                    // mixing will restart at this->p_chgmix->mixing_restart steps
                    if (drho <= GlobalV::MIXING_RESTART && GlobalV::MIXING_RESTART > 0.0 && this->p_chgmix->mixing_restart > iter)
                    {
                        this->p_chgmix->mixing_restart = iter + 1;
                    }
                    // drho will be 0 at this->p_chgmix->mixing_restart step, which is not ground state
                    this->conv_elec = (drho < this->scf_thr && !(iter==this->p_chgmix->mixing_restart && GlobalV::MIXING_RESTART > 0.0));

                    // If drho < hsolver_error in the first iter or drho < scf_thr, we do not change rho.
                    if (drho < hsolver_error || this->conv_elec)
                    {
                        if (drho < hsolver_error)    GlobalV::ofs_warning << " drho < hsolver_error, keep charge density unchanged." << std::endl;
                    }
                    else
                    {
                        //----------charge mixing---------------
                        // mixing will restart after this->p_chgmix->mixing_restart steps
                        if (GlobalV::MIXING_RESTART > 0 && iter == this->p_chgmix->mixing_restart - 1)
                        {
                            // do not mix charge density
                        }
                        else
                        {
                            p_chgmix->mix_rho(pelec->charge); // update chr->rho by mixing
                        }
                        if (GlobalV::SCF_THR_TYPE == 2) pelec->charge->renormalize_rho(); // renormalize rho in R-space would induce a error in K-space
                        //----------charge mixing done-----------  
                    }
                }
#ifdef __MPI
		        MPI_Bcast(&drho, 1, MPI_DOUBLE , 0, PARAPW_WORLD);
		        MPI_Bcast(&this->conv_elec, 1, MPI_DOUBLE , 0, PARAPW_WORLD);
                MPI_Bcast(pelec->charge->rho[0], this->pw_rhod->nrxx, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif

                // Hamilt should be used after it is constructed.
                // this->phamilt->update(conv_elec);
                updatepot(istep, iter);
                eachiterfinish(iter);
#ifdef __MPI
                double duration = (double)(MPI_Wtime() - iterstart);
#else
                double duration = (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - iterstart)).count() / static_cast<double>(1e6);
#endif
                /*
                    SCF print: G1    -3.435545e+03  0.000000e+00   3.607e-01  2.862e-01
                */
                printiter(iter, drho, duration, diag_ethr);
                if (this->conv_elec)
                {
                    this->niter = iter;
                    bool stop = this->do_after_converge(iter);
                    if(stop) break;
                }
                // notice for restart
                if (GlobalV::MIXING_RESTART > 0 && iter == this->p_chgmix->mixing_restart - 1)
                {
                    std::cout<<"SCF restart after this step!"<<std::endl;
                }
            }
            afterscf(istep);

            ModuleBase::timer::tick(this->classname, "Run");
        }       
        // add nkstot,nkstot_ibz to output json
        #ifdef __RAPIDJSON
                int Jnkstot = this->pelec->klist->nkstot;
                int Jnkstot_ibz = this->pelec->klist->nkstot_ibz;
                Json::add_nkstot(Jnkstot,Jnkstot_ibz);
        #endif //__RAPIDJSON 
        return;
    };

    template<typename T, typename Device>
    void ESolver_KS<T, Device>::printhead()
    {
        std::cout << " " << std::setw(7) << "ITER";
        if (GlobalV::NSPIN == 2)
        {
            std::cout << std::setw(10) << "TMAG";
            std::cout << std::setw(10) << "AMAG";
        }
        std::cout << std::setw(15) << "ETOT(eV)";
        std::cout << std::setw(15) << "EDIFF(eV)";
        std::cout << std::setw(11) << "DRHO";
        std::cout << std::setw(11) << "TIME(s)" << std::endl;
    }

    template<typename T, typename Device>
    void ESolver_KS<T, Device>::printiter(const int iter, const double drho, const double duration, const double ethr)
    {
        this->pelec->print_etot(this->conv_elec, iter, drho, duration, INPUT.printe, ethr);
    }

    template<typename T, typename Device>
    void ESolver_KS<T, Device>::writehead(std::ofstream& ofs_running, const int istep, const int iter)
    {
        ofs_running
            << "\n "
            << this->basisname
            << " ALGORITHM --------------- ION=" << std::setw(4) << istep + 1
            << "  ELEC=" << std::setw(4) << iter
            << "--------------------------------\n";
    }

    template<typename T, typename Device>
    int ESolver_KS<T, Device>::getniter()
    {
        return this->niter;
    }

    template<typename T, typename Device>
    ModuleIO::Output_Rho ESolver_KS<T, Device>::create_Output_Rho(int is, int iter, const std::string& prefix)
    {
        int precision = 3;
        std::string tag = "CHG";
        return ModuleIO::Output_Rho(this->pw_big,
                                    this->pw_rhod,
                                    is,
                                    GlobalV::NSPIN,
                                    pelec->charge->rho_save[is],
                                    iter,
                                    this->pelec->eferm.get_efval(is),
                                    &(GlobalC::ucell),
                                    GlobalV::global_out_dir,
                                    precision,
                                    tag,
                                    prefix);
    }

    template<typename T, typename Device>
    ModuleIO::Output_Rho ESolver_KS<T, Device>::create_Output_Kin(int is, int iter, const std::string& prefix)
    {
        int precision = 11;
        std::string tag = "TAU";
        return ModuleIO::Output_Rho(this->pw_big,
                                    this->pw_rhod,
                                    is,
                                    GlobalV::NSPIN,
                                    pelec->charge->kin_r_save[is],
                                    iter,
                                    this->pelec->eferm.get_efval(is),
                                    &(GlobalC::ucell),
                                    GlobalV::global_out_dir,
                                    precision,
                                    tag,
                                    prefix);
    }

    template<typename T, typename Device>
    ModuleIO::Output_Potential ESolver_KS<T, Device>::create_Output_Potential(int iter, const std::string& prefix)
    {
        int precision = 3;
        std::string tag = "POT";
        return ModuleIO::Output_Potential(this->pw_big,
                                          this->pw_rhod,
                                          GlobalV::NSPIN,
                                          iter,
                                          GlobalV::out_pot,
                                          this->pelec->pot->get_effective_v(),
                                          this->pelec->pot->get_fixed_v(),
                                          &(GlobalC::ucell),
                                          pelec->charge,
                                          precision,
                                          GlobalV::global_out_dir,
                                          tag,
                                          prefix);
    }



template class ESolver_KS<std::complex<float>, psi::DEVICE_CPU>;
template class ESolver_KS<std::complex<double>, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS<std::complex<float>, psi::DEVICE_GPU>;
template class ESolver_KS<std::complex<double>, psi::DEVICE_GPU>;
#endif
#ifdef __LCAO
template class ESolver_KS<double, psi::DEVICE_CPU>;
#endif
}