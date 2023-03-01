/**********************************************************************
  Mulliken_Charge.cpp:

     Mulliken_Charge.cpp is a subrutine to calculate Mulliken charge.
 
  Log of Mulliken_Charge.cpp:

     12/Oct/2018  Released by Feng Qi

***********************************************************************/

#include "mulliken_charge.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_base/scalapack_connector.h"
#include "write_orb_info.h"
#ifdef __LCAO
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_gen_fixedH.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/global_fp.h"
#endif
#include "module_cell/module_neighbor/sltk_atom_arrange.h"

ModuleBase::matrix Mulliken_Charge::cal_mulliken(const std::vector<ModuleBase::matrix> &dm,
    LCAO_Hamilt &uhm
)
{
    ModuleBase::TITLE("Mulliken_Charge", "cal_mulliken");

    const int nspin = (GlobalV::NSPIN == 2) ? 2 : 1;
    const int nlocal = (GlobalV::NSPIN == 4) ? GlobalV::NLOCAL/2 : GlobalV::NLOCAL;
    // std::vector<std::vector<double>> MecMulP(GlobalV::NSPIN, std::vector<double>(nlocal, 0));
    // std::vector<std::vector<double>> orbMulP(GlobalV::NSPIN, std::vector<double>(nlocal, 0));
    ModuleBase::matrix MecMulP, orbMulP;
    MecMulP.create(GlobalV::NSPIN, nlocal);
    orbMulP.create(GlobalV::NSPIN, nlocal);

    for(size_t is=0; is!=nspin; ++is)
    {
        ModuleBase::matrix mud;
        mud.create(uhm.LM->ParaV->nrow, uhm.LM->ParaV->ncol);
#ifdef __MPI
        const char T_char = 'T';
        const char N_char = 'N';
        const int one_int = 1;
        const double one_float[2] = {1.0, 0.0}, zero_float[2] = {0.0, 0.0};        
        pdgemm_(&T_char,
                &N_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                dm[is].c,
                &one_int,
                &one_int,
                uhm.LM->ParaV->desc,
                uhm.LM->Sloc.data(),
                &one_int,
                &one_int,
                uhm.LM->ParaV->desc,
                &zero_float[0],
                mud.c,
                &one_int,
                &one_int,
                uhm.LM->ParaV->desc);
        if(GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
        {
            for(size_t i=0; i!=GlobalV::NLOCAL; ++i)
                if(uhm.LM->ParaV->in_this_processor(i, i))
                {
                    const int ir = uhm.LM->ParaV->trace_loc_row[i];
                    const int ic = uhm.LM->ParaV->trace_loc_col[i];
                    MecMulP(is, i) += mud(ir, ic);
                }
        }
        else if(GlobalV::NSPIN == 4)
        {
            for(size_t i=0; i!=GlobalV::NLOCAL; ++i)
            {
                const int index = i%2;
                if(!index)
                {
                    const int j = i/2;
                    const int k1 = 2*j;
                    const int k2 = 2*j+1;
                    if(uhm.LM->ParaV->in_this_processor(k1, k1))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k1];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k1];
                        MecMulP(0, j) += mud(ir, ic);
                    }
                    if(uhm.LM->ParaV->in_this_processor(k1, k2))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k1];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k2];
                        MecMulP(1, j) += mud(ir, ic);
                    }
                    if(uhm.LM->ParaV->in_this_processor(k2, k1))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k2];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k1];
                        MecMulP(2, j) += mud(ir, ic);
                    }
                    if(uhm.LM->ParaV->in_this_processor(k2, k2))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k2];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k2];
                        MecMulP(3, j) += mud(ir, ic);
                    }
                }
            }
        }
#endif
    }
#ifdef __MPI 
    MPI_Reduce(MecMulP.c, orbMulP.c, GlobalV::NSPIN*nlocal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif 

    return orbMulP;
}

ModuleBase::matrix Mulliken_Charge::cal_mulliken_k(const std::vector<ModuleBase::ComplexMatrix> &dm,
    LCAO_Hamilt &uhm
)
{
    ModuleBase::TITLE("Mulliken_Charge", "cal_mulliken");

    const int nspin = (GlobalV::NSPIN == 2) ? 2 : 1;
    const int nlocal = (GlobalV::NSPIN == 4) ? GlobalV::NLOCAL/2 : GlobalV::NLOCAL;
    // std::vector<std::vector<double>> MecMulP(GlobalV::NSPIN, std::vector<double>(nlocal, 0));
    // std::vector<std::vector<double>> orbMulP(GlobalV::NSPIN, std::vector<double>(nlocal, 0));
    ModuleBase::matrix MecMulP, orbMulP;
    MecMulP.create(GlobalV::NSPIN, nlocal);
    orbMulP.create(GlobalV::NSPIN, nlocal);

    GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
		GlobalV::ofs_running,
		GlobalV::OUT_LEVEL,
		GlobalC::ORB.get_rcutmax_Phi(), 
		GlobalC::ucell.infoNL.get_rcutmax_Beta(), 
		GlobalV::GAMMA_ONLY_LOCAL);
	atom_arrange::search(
		GlobalV::SEARCH_PBC,
		GlobalV::ofs_running,
		GlobalC::GridD, 
		GlobalC::ucell, 
		GlobalV::SEARCH_RADIUS,
		GlobalV::test_atom_input);
	uhm.LM->allocate_HS_R(uhm.LM->ParaV->nnr);
	uhm.LM->zeros_HSR('S');
	uhm.genH.calculate_S_no(uhm.LM->SlocR.data());
	uhm.genH.build_ST_new('S', false, GlobalC::ucell, uhm.LM->SlocR.data());

    for(size_t ik = 0; ik != GlobalC::kv.nks; ++ik)
    {
        uhm.LM->allocate_HS_k(uhm.LM->ParaV->nloc);
        uhm.LM->zeros_HSk('S');
		uhm.LM->folding_fixedH(ik);

        ModuleBase::ComplexMatrix mud;
        mud.create(uhm.LM->ParaV->nrow, uhm.LM->ParaV->ncol);

#ifdef __MPI
        const char T_char = 'T';
        const char N_char = 'N';
        const int one_int = 1;
        const double one_float[2] = {1.0, 0.0}, zero_float[2] = {0.0, 0.0};        
        pzgemm_(&T_char,
                &N_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                dm[ik].c,
                &one_int,
                &one_int,
                uhm.LM->ParaV->desc,
                uhm.LM->Sloc2.data(),
                &one_int,
                &one_int,
                uhm.LM->ParaV->desc,
                &zero_float[0],
                mud.c,
                &one_int,
                &one_int,
                uhm.LM->ParaV->desc);
        if(GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
        {
            const int spin = GlobalC::kv.isk[ik];
            for(size_t i=0; i!=GlobalV::NLOCAL; ++i)
                if(uhm.LM->ParaV->in_this_processor(i, i))
                {
                    const int ir = uhm.LM->ParaV->trace_loc_row[i];
                    const int ic = uhm.LM->ParaV->trace_loc_col[i];
                    MecMulP(spin, i) += mud(ir, ic).real();
                }
        }
        else if(GlobalV::NSPIN == 4)
        {
            for(size_t i=0; i!=GlobalV::NLOCAL; ++i)
            {
                const int index = i%2;
                if(!index)
                {
                    const int j = i/2;
                    const int k1 = 2*j;
                    const int k2 = 2*j+1;
                    if(uhm.LM->ParaV->in_this_processor(k1, k1))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k1];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k1];
                        MecMulP(0, j) += mud(ir, ic).real();
                    }
                    if(uhm.LM->ParaV->in_this_processor(k1, k2))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k1];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k2];
                        MecMulP(1, j) += mud(ir, ic).real();
                    }
                    if(uhm.LM->ParaV->in_this_processor(k2, k1))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k2];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k1];
                        MecMulP(2, j) += mud(ir, ic).real();
                    }
                    if(uhm.LM->ParaV->in_this_processor(k2, k2))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k2];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k2];
                        MecMulP(3, j) += mud(ir, ic).real();
                    }
                }
            }
        }
#endif
    }
#ifdef __MPI
    atom_arrange::delete_vector(
		GlobalV::ofs_running,
		GlobalV::SEARCH_PBC,
		GlobalC::GridD, 
		GlobalC::ucell, 
		GlobalV::SEARCH_RADIUS, 
		GlobalV::test_atom_input); 
    MPI_Reduce(MecMulP.c, orbMulP.c, GlobalV::NSPIN*nlocal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif 

    return orbMulP;
}

std::vector<std::vector<std::vector<double>>> Mulliken_Charge::convert(const ModuleBase::matrix &orbMulP)
{
    std::vector<std::vector<std::vector<double>>> AorbMulP;
    AorbMulP.resize(GlobalV::NSPIN);
    
    for (size_t is=0; is!=GlobalV::NSPIN; ++is)
    {
        int num=0;
        AorbMulP[is].resize(GlobalC::ucell.nat);
        for (size_t i=0; i!=GlobalC::ucell.nat; ++i)
        {   
            const int a = GlobalC::ucell.iat2ia[i];
            const int t = GlobalC::ucell.iat2it[i];
            AorbMulP[is][i].resize(GlobalC::ucell.atoms[t].nw);
            for (size_t j=0; j!=GlobalC::ucell.atoms[t].nw; ++j)
            {
                AorbMulP[is][i][j] = orbMulP(is, num);
                num++;
            }
        }
    }
    
    return AorbMulP;
}

void Mulliken_Charge::out_mulliken(LCAO_Hamilt &uhm, Local_Orbital_Charge &loc)
{
    ModuleBase::TITLE("Mulliken_Charge", "out_mulliken");

    ModuleBase::matrix orbMulP;
    if(GlobalV::GAMMA_ONLY_LOCAL)
        orbMulP = this->cal_mulliken(loc.dm_gamma, uhm);
    else
        orbMulP = this->cal_mulliken_k(loc.dm_k, uhm);

    std::vector<std::vector<std::vector<double>>> AorbMulP = this->convert(orbMulP);

    if(GlobalV::MY_RANK == 0)
    {
        const int nlocal = (GlobalV::NSPIN == 4) ? GlobalV::NLOCAL/2 : GlobalV::NLOCAL;
        std::stringstream as;
        as << GlobalV::global_out_dir << "mulliken.txt";
        std::ofstream os(as.str().c_str());
        os << "CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM" << std::endl;

		double sch = 0.0;
		os << std::setprecision(8);
		for(size_t is=0; is!=GlobalV::NSPIN; ++is)
		{
			double sss = 0.0;
			for(size_t iw=0; iw!=nlocal; ++iw)
			{
				sch += orbMulP(is, iw);
				sss += orbMulP(is, iw);
			}
			os << " Total charge of spin " << is+1 << ":\t" << sss << std::endl;
		}
		os << " Total charge:\t" << sch << std::endl;
		os << "Decomposed Mulliken populations" << std::endl;

        for (size_t i = 0; i != GlobalC::ucell.nat; ++i)
        {
            double total_charge = 0.0, atom_mag = 0.0;
            std::vector<double> total_charge_soc(GlobalV::NSPIN);
            const int t = GlobalC::ucell.iat2it[i];
            int num = 0;
            if (GlobalV::NSPIN==1)
                os << i << std::setw(25) << "Zeta of " << GlobalC::ucell.atoms[t].label << std::setw(30) << "Spin 1" << std::endl;
            else if (GlobalV::NSPIN==2)
                os << i << std::setw(25) << "Zeta of " << GlobalC::ucell.atoms[t].label << std::setw(30) << "Spin 1" << std::setw(30) << "Spin 2" << std::setw(30) << "Sum" << std::setw(30) << "Diff" << std::endl;
            else if (GlobalV::NSPIN==4)
                os << i << std::setw(25) << "Zeta of " << GlobalC::ucell.atoms[t].label << std::setw(30) << "Spin 1" << std::setw(30) << "Spin 2" << std::setw(30) << "Spin 3" << std::setw(30) << "Spin 4" << std::endl;
        
            for (size_t L=0; L!=GlobalC::ucell.atoms[t].nwl+1; ++L)
            {
                std::vector<double> sum_l(GlobalV::NSPIN, 0.0);
                for (size_t Z=0; Z!=GlobalC::ucell.atoms[t].l_nchi[L]; ++Z)
                {
                    std::vector<double> sum_m(GlobalV::NSPIN, 0.0);
                    for (size_t M=0; M!=(2*L+1); ++M)
                    {
                        if (GlobalV::NSPIN==1)
                        {
                            double spin1 = AorbMulP[0][i][num];
                            os << GlobalC::en.Name_Angular[L][M] << std::setw(25) << Z << std::setw(32) << spin1 << std::endl;
                            sum_m[0] += spin1;
                        }
                        else if (GlobalV::NSPIN==2)
                        {
                            double spin1 = AorbMulP[0][i][num]; 
                            double spin2 = AorbMulP[1][i][num];
                            double sum = spin1 + spin2;
                            double diff = spin1 - spin2;
                            os << GlobalC::en.Name_Angular[L][M] << std::setw(25) << Z << std::setw(32) << spin1 << std::setw(30) << spin2 << std::setw(30) << sum << std::setw(30) << diff << std::endl;
                            sum_m[0] += AorbMulP[0][i][num];
                            sum_m[1] += AorbMulP[1][i][num];
                        }
                        else if (GlobalV::NSPIN==4)
                        {
                            double spin1 = AorbMulP[0][i][num]; 
                            double spin2 = AorbMulP[1][i][num]; 
                            double spin3 = AorbMulP[2][i][num]; 
                            double spin4 = AorbMulP[3][i][num]; 
                            os << GlobalC::en.Name_Angular[L][M] << std::setw(25) << Z << std::setw(32) << spin1 << std::setw(30) << spin2 << std::setw(30) << spin3 << std::setw(30) << spin4 << std::endl;
                            sum_m[0] += spin1;
                            sum_m[1] += spin2;
                            sum_m[2] += spin3;
                            sum_m[3] += spin4;
                        }
                        num++;
                    }

                    if (GlobalV::NSPIN==1)
                    {
                        double spin1 = sum_m[0];
                        os << "  sum over m "<< std::setw(45) << spin1 << std::endl;
                        sum_l[0] += spin1;
                    }
                    else if (GlobalV::NSPIN==2)
                    {
                        double spin1 = sum_m[0];
                        double spin2 = sum_m[1];
                        double sum = spin1 + spin2;
                        double diff = spin1 - spin2;
                        os << "  sum over m "<< std::setw(45) << spin1 << std::setw(30) << spin2 << std::setw(35) << sum << std::setw(25) << diff << std::endl;
                        sum_l[0] += spin1;
                        sum_l[1] += spin2;
                    }
                    else if (GlobalV::NSPIN==4)
                    {
                        double spin1 = sum_m[0];
                        double spin2 = sum_m[1];
                        double spin3 = sum_m[2];
                        double spin4 = sum_m[3];
                        os << "  sum over m "<< std::setw(45) << spin1 << std::setw(30) << spin2 << std::setw(30) << spin3 << std::setw(30) << spin4 << std::endl;
                        sum_l[0] += spin1;
                        sum_l[1] += spin2;
                        sum_l[2] += spin3;
                        sum_l[3] += spin4;
                    }
                }
            
                if(GlobalC::ucell.atoms[t].l_nchi[L])
                {
                    if (GlobalV::NSPIN==1)
                    {
                        double spin1 = sum_l[0];
                        os << "  sum over m+zeta "<< std::setw(40) << spin1 << std::endl;
                        total_charge += sum_l[0];
                    }
                    else if (GlobalV::NSPIN==2)
                    {
                        double spin1 = sum_l[0];
                        double spin2 = sum_l[1];
                        double sum = spin1 + spin2;
                        double diff = spin1 - spin2;
                        os << "  sum over m+zeta "<< std::setw(40) << spin1 << std::setw(30) << spin2 << std::setw(35) << sum << std::setw(25) << diff << std::endl;
                        total_charge += sum_l[0] + sum_l[1];
                        atom_mag += sum_l[0] - sum_l[1];
                    }
                    else if (GlobalV::NSPIN==4)
                    {
                        double spin1 = sum_l[0];
                        double spin2 = sum_l[1];
                        double spin3 = sum_l[2];
                        double spin4 = sum_l[3];
                        os << "  sum over m+zeta "<< std::setw(40) << spin1 << std::setw(30) << spin2 << std::setw(30) << spin3 << std::setw(30) << spin4 << std::endl;
                        total_charge_soc[0] += spin1;
                        total_charge_soc[1] += spin2;
                        total_charge_soc[2] += spin3;
                        total_charge_soc[3] += spin4;
                    }
                }
            }

            if (GlobalV::NSPIN==1)
                os << "Total Charge on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(20) << total_charge <<std::endl;
            else if (GlobalV::NSPIN==2)
            {
                os << "Total Charge on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(20) << total_charge <<std::endl;
                os << "Total Magnetism on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(20) << atom_mag <<std::endl;
            }
            else if (GlobalV::NSPIN==4)
            {
                double spin1 = total_charge_soc[0];
                double spin2 = total_charge_soc[1];
                double spin3 = total_charge_soc[2];
                double spin4 = total_charge_soc[3];
                os << "Total Charge on atom in four components: " << GlobalC::ucell.atoms[t].label <<  std::setw(20) 
                << "(" << spin1 << ", " << spin2 << ", " << spin3 << ", " << spin4 << ")" 
                <<std::endl;
            }
            os << std::endl <<std::endl;
        }
        os.close();
        ModuleIO::write_orb_info();
    }
}