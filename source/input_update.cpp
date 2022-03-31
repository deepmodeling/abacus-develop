//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-11-6
//==========================================================
#include "input_update.h"
#include "src_pw/global.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "input.h"
#include "src_ions/ions_move_basic.h"
#include "src_io/optical.h"
#ifdef __LCAO
#include "src_lcao/FORCE_STRESS.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/global_fp.h" // mohan update 2021-01-30
#endif

Update_input::Update_input() {}
Update_input::~Update_input() {}

void Update_input::init(const std::string &fn)
{
	GlobalV::ofs_warning << "\n CHECK UPDATE INPUT INFORMATION" << std::endl;
    this->Read(fn);
#ifdef __MPI
    Bcast();
#endif

    return;
}

bool Update_input::Read(const std::string &fn)
{
    if (GlobalV::MY_RANK!=0) return false;

    std::ifstream ifs(fn.c_str(), ios::in);	// "in_datas/input_parameters"

    if (!ifs) 
	{
		return false;
	}

    ifs.clear();
    ifs.seekg(0);

    char word[80];
    char word1[80];
    int ierr = 0;

    //ifs >> std::setiosflags(ios::uppercase);
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word;
        ifs.ignore(150, '\n');
        if (strcmp(word , "INPUT_PARAMETERS") == 0)
        {
            ierr = 1;
            break;
        }
        ifs.rdstate();
    }

    if (ierr == 0)
    {
		std::cout << " Error parameter list." << std::endl;
		return false;// return error : false
    }

    ifs.rdstate();

    while (ifs.good())
    {
        ifs >> word1;
        INPUT.strtolower(word1, word);
		
		// 1
        if (strcmp("cal_force", word) == 0)
        {
            read_value(ifs, cal_force);
			if(cal_force!=GlobalV::CAL_FORCE)
			{
				this->change(GlobalV::ofs_warning,"CAL_FORCE",GlobalV::CAL_FORCE, cal_force);
    			GlobalV::CAL_FORCE = this->cal_force;
			}
        }
		// 2
        else if (strcmp("force_thr", word) == 0)
        {
            read_value(ifs, force_thr);
			if(force_thr!=GlobalV::FORCE_THR)
			{
				this->change(GlobalV::ofs_warning,"forth_thr(Ry/Bohr)",GlobalV::FORCE_THR,force_thr);
    			GlobalV::FORCE_THR = this->force_thr;
			}
        }
		// 2
        else if (strcmp("force_thr_ev", word) == 0)
        {
            read_value(ifs, force_thr);
			force_thr = force_thr / 13.6058 * 0.529177;
			if(force_thr!=GlobalV::FORCE_THR)
			{
				this->change(GlobalV::ofs_warning,"force_thr(Ry/Bohr)",GlobalV::FORCE_THR,force_thr);
    			GlobalV::FORCE_THR = this->force_thr;
			}
        }
		// 3
#ifdef __FP
        else if (strcmp("force_thr_ev2", word) == 0)
        {
            read_value(ifs, force_thr_ev2);
			if(force_thr_ev2!=Force_Stress_LCAO::force_invalid_threshold_ev)
			{
				this->change(GlobalV::ofs_warning,"force_thr threshold(Ry/Bohr)",Force_Stress_LCAO::force_invalid_threshold_ev,force_thr_ev2);
    			Force_Stress_LCAO::force_invalid_threshold_ev = this->force_thr_ev2;
			}
        }
#endif
		// 4
        else if (strcmp("scf_thr_rho", word) == 0)
        {
            read_value(ifs, scf_thr_rho);
			if(scf_thr_rho!=GlobalV::SCF_THR_RHO)
			{
				this->change(GlobalV::ofs_warning,"scf_thr_rho",GlobalV::SCF_THR_RHO,scf_thr_rho);
    			GlobalV::SCF_THR_RHO = this->scf_thr_rho;
			}
        }
		// 5
        else if (strcmp("scf_nmax", word) == 0)
        {
            read_value(ifs, scf_nmax);
			if(scf_nmax!=GlobalV::SCF_NMAX)
			{
				this->change(GlobalV::ofs_warning,"scf_nmax",GlobalV::SCF_NMAX,scf_nmax);
    			GlobalV::SCF_NMAX = this->scf_nmax;
			}
        }
		// 6
        else if (strcmp("relax_nmax", word) == 0)
        {
            read_value(ifs, relax_nmax);
			if(relax_nmax!=GlobalV::RELAX_NMAX)
			{
				this->change(GlobalV::ofs_warning,"relax_nmax",GlobalV::RELAX_NMAX,relax_nmax);
    			GlobalV::RELAX_NMAX = this->relax_nmax;
			}
        }
        else if (strcmp("md_nstep", word) == 0)
        {
            read_value(ifs, md_nstep);
			if(md_nstep!=GlobalV::MD_NSTEP)
			{
				this->change(GlobalV::ofs_warning, "md_nstep", GlobalV::MD_NSTEP, md_nstep);
    			GlobalV::MD_NSTEP = this->md_nstep;
			}
        }
		// 7
        else if (strcmp("mixing_beta", word) == 0)
        {
            read_value(ifs, mixing_beta);
			if(mixing_beta!=GlobalC::CHR.mixing_beta)
			{
				this->change(GlobalV::ofs_warning,"mixing_beta",GlobalC::CHR.mixing_beta,mixing_beta);
    			GlobalC::CHR.mixing_beta = mixing_beta;
			}
        }
		// 8
        else if (strcmp("printe", word) == 0)
        {
            read_value(ifs, printe);
			if(printe!=GlobalC::en.printe)
			{
				this->change(GlobalV::ofs_warning,"printe",GlobalC::en.printe,printe);
				GlobalC::en.printe = this->printe; // mohan add 2011-03-16
			}
        }
		// 9 
        else if (strcmp("chg_extrap", word) == 0)//xiaohui modify 2015-02-01
        {
            read_value(ifs, chg_extrap);

			if(chg_extrap!=GlobalC::pot.chg_extrap)
			{
				this->change(GlobalV::ofs_warning,"chg_extrap",GlobalC::pot.chg_extrap,chg_extrap);
				GlobalC::pot.chg_extrap = this->chg_extrap;
				if(out_dm==0) 
				{
					out_dm = 10000;
#ifdef __FP
					this->change(GlobalV::ofs_warning,"out_dm",Local_Orbital_Charge::out_dm, out_dm);
					Local_Orbital_Charge::out_dm=out_dm;
#endif
				}
			}

        }
		// 10
        else if (strcmp("out_chg", word) == 0)
        {
            read_value(ifs, out_chg);
			if(out_chg!=GlobalC::CHR.out_chg)
			{
				this->change(GlobalV::ofs_warning,"out_chg",GlobalC::CHR.out_chg,out_chg);
    			GlobalC::CHR.out_chg = this->out_chg;
			}
        }
		// 11
        else if (strcmp("out_dm", word) == 0)
        {
            read_value(ifs, out_dm);
#ifdef __FP
			if(out_dm!=Local_Orbital_Charge::out_dm)
			{
				this->change(GlobalV::ofs_warning,"out_dm",Local_Orbital_Charge::out_dm,out_dm);
				Local_Orbital_Charge::out_dm = this->out_dm;
			}
#endif
        }
		// 12
        else if (strcmp("out_dos", word) == 0)
        {
            read_value(ifs, out_dos);
			if(out_dos!=GlobalC::en.out_dos)
			{
				this->change(GlobalV::ofs_warning,"out_dos",GlobalC::en.out_dos,out_dos);
				GlobalC::en.out_dos = this->out_dos;
			}
        }
		// 13
        else if (strcmp("out_wfc_lcao", word) == 0)
        {
            read_value(ifs, out_wfc_lcao);
#ifdef __FP
			//if(out_wfc_lcao!=out_wfc_lcao)
			if(out_wfc_lcao!=Pdiag_Double::out_wfc_lcao)		// Peize Lin change out_wfc_lcao to GlobalC::ParaO.out_wfc_lcao at 2020.01.31
			{
				this->change(GlobalV::ofs_warning,"out_wfc_lcao",Pdiag_Double::out_wfc_lcao,out_wfc_lcao);
				Pdiag_Double::out_wfc_lcao = this->out_wfc_lcao;
			}
#endif
        }
        else
        {
            //std::cout << " THE PARAMETER NAME '" << word
            //     << "' IS NOT USED!" << std::endl;
            ifs.ignore(150, '\n');
        }

        ifs.rdstate();
        if (ifs.eof() != 0)
        {
			break;
        }
        else if (ifs.bad() != 0)
        {
			std::cout << " Bad input parameters. " << std::endl;
            return false;
        }
        else if (ifs.fail() != 0)
        {
			std::cout << " word = " << word << std::endl;
			std::cout << " Fail to read parameters. " << std::endl; 
            ifs.clear();
			return false;
        }
        else if (ifs.good() == 0)
        {
			break;
        }
    }

    return true;
}//end read_parameters

#include "src_parallel/parallel_common.h"
#ifdef __MPI
void Update_input::Bcast()
{
    Parallel_Common::bcast_int( GlobalV::CAL_FORCE );
    Parallel_Common::bcast_double( GlobalV::FORCE_THR);
#ifdef __LCAO
    Parallel_Common::bcast_double( Force_Stress_LCAO::force_invalid_threshold_ev);
    Parallel_Common::bcast_int( Local_Orbital_Charge::out_dm );
    Parallel_Common::bcast_int( Pdiag_Double::out_wfc_lcao );
#endif
    Parallel_Common::bcast_double( GlobalV::SCF_THR_RHO );
    Parallel_Common::bcast_int( GlobalV::SCF_NMAX );
    Parallel_Common::bcast_int( GlobalV::RELAX_NMAX );
    Parallel_Common::bcast_double( GlobalC::CHR.mixing_beta );
    Parallel_Common::bcast_int( GlobalC::en.printe );
    Parallel_Common::bcast_string( GlobalC::pot.chg_extrap );//xiaohui modify 2015-02-01
    Parallel_Common::bcast_int( GlobalC::CHR.out_chg );
	Parallel_Common::bcast_int( GlobalC::en.out_dos );
    Parallel_Common::bcast_double( GlobalC::CHR.nelec );
	
    Parallel_Common::bcast_int( GlobalV::MD_NSTEP );
    return;
}
#endif

