/**********************************************************************
  Mulliken_Charge.cpp:

     Mulliken_Charge.cpp is a subrutine to calculate Mulliken charge.
 
  Log of Mulliken_Charge.cpp:

     12/Oct/2018  Released by Feng Qi
     03/2023/     Refactored by Yuyang Ji
     18/04/2023   Convert to namespace by Yu Liu

***********************************************************************/

#include "mulliken_charge.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/name_angular.h"
#include "module_base/scalapack_connector.h"
#include "write_orb_info.h"
#include "module_elecstate/elecstate_lcao.h"

template<>
ModuleBase::matrix ModuleIO::cal_mulliken(const std::vector<std::vector<double>> &dm,
    LCAO_Matrix* LM, const K_Vectors& kv, hamilt::Hamilt<double>* ham_in
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
        mud.create(LM->ParaV->ncol, LM->ParaV->nrow);
#ifdef __MPI
        const char T_char = 'T';
        const char N_char = 'N';
        const int one_int = 1;
        const double one_float = 1.0, zero_float = 0.0;        
        pdgemm_(&T_char,
                &T_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float,
                dm[is].data(),
                &one_int,
                &one_int,
                LM->ParaV->desc,
                LM->Sloc.data(),
                &one_int,
                &one_int,
                LM->ParaV->desc,
                &zero_float,
                mud.c,
                &one_int,
                &one_int,
                LM->ParaV->desc);
        if(GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
        {
            for(size_t i=0; i!=GlobalV::NLOCAL; ++i)
                if(LM->ParaV->in_this_processor(i, i))
                {
                    const int ir = LM->ParaV->global2local_row(i);
                    const int ic = LM->ParaV->global2local_col(i);
                    MecMulP(is, i) += mud(ic, ir);
                }
        }
        //NSPIN=4 is forbidden for gamma_only case
        /*else if(GlobalV::NSPIN == 4)
        {
            for(size_t i=0; i!=GlobalV::NLOCAL; ++i)
            {
                const int index = i%2;
                if(!index)
                {
                    const int j = i/2;
                    const int k1 = 2*j;
                    const int k2 = 2*j+1;
                    if(LM->ParaV->in_this_processor(k1, k1))
                    {
                        const int ir = LM->ParaV->global2local_row(k1);
                        const int ic = LM->ParaV->global2local_col(k1);
                        MecMulP(0, j) += mud(ir, ic);
                    }
                    if(LM->ParaV->in_this_processor(k1, k2))
                    {
                        const int ir = LM->ParaV->global2local_row(k1);
                        const int ic = LM->ParaV->global2local_col(k2);
                        MecMulP(1, j) += mud(ir, ic);
                    }
                    if(LM->ParaV->in_this_processor(k2, k1))
                    {
                        const int ir = LM->ParaV->global2local_row(k2);
                        const int ic = LM->ParaV->global2local_col(k1);
                        MecMulP(2, j) += mud(ir, ic);
                    }
                    if(LM->ParaV->in_this_processor(k2, k2))
                    {
                        const int ir = LM->ParaV->global2local_row(k2);
                        const int ic = LM->ParaV->global2local_col(k2);
                        MecMulP(3, j) += mud(ir, ic);
                    }
                }
            }
        }*/
#endif
    }
#ifdef __MPI 
    MPI_Reduce(MecMulP.c, orbMulP.c, GlobalV::NSPIN*nlocal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif 

    return orbMulP;
}

template<>
ModuleBase::matrix ModuleIO::cal_mulliken(const std::vector<std::vector<std::complex<double>>>& dm,
    LCAO_Matrix* LM, const K_Vectors& kv, hamilt::Hamilt<std::complex<double>>* ham_in
)
{
    ModuleBase::TITLE("Mulliken_Charge", "cal_mulliken_k");

    const int nspin = (GlobalV::NSPIN == 2) ? 2 : 1;
    const int nlocal = (GlobalV::NSPIN == 4) ? GlobalV::NLOCAL/2 : GlobalV::NLOCAL;
    // std::vector<std::vector<double>> MecMulP(GlobalV::NSPIN, std::vector<double>(nlocal, 0));
    // std::vector<std::vector<double>> orbMulP(GlobalV::NSPIN, std::vector<double>(nlocal, 0));
    ModuleBase::matrix MecMulP, orbMulP;
    MecMulP.create(GlobalV::NSPIN, nlocal, true);
    orbMulP.create(GlobalV::NSPIN, nlocal, true);

    for(size_t ik = 0; ik != kv.nks; ++ik)
    {
        // calculate SK for current k point
        // the target matrix is LM->Sloc2 with collumn-major
        if(GlobalV::NSPIN == 4)
        {
            dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(ham_in)->updateSk(ik, LM, 1);
        }
        else
        {
            dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(ham_in)->updateSk(ik, LM, 1);
        }

        ModuleBase::ComplexMatrix mud;
        mud.create(LM->ParaV->ncol, LM->ParaV->nrow, true);

#ifdef __MPI
        const char T_char = 'T';
        const char N_char = 'N';
        const int one_int = 1;
        const std::complex<double> one_float = {1.0, 0.0}, zero_float = {0.0, 0.0};        
        pzgemm_(&T_char,
                &T_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float,
                dm[ik].data(),
                &one_int,
                &one_int,
                LM->ParaV->desc,
                LM->Sloc2.data(),
                &one_int,
                &one_int,
                LM->ParaV->desc,
                &zero_float,
                mud.c,
                &one_int,
                &one_int,
                LM->ParaV->desc);
        if(GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
        {
            const int spin = kv.isk[ik];
            for(size_t i=0; i!=GlobalV::NLOCAL; ++i)
                if(LM->ParaV->in_this_processor(i, i))
                {
                    const int ir = LM->ParaV->global2local_row(i);
                    const int ic = LM->ParaV->global2local_col(i);
                    MecMulP(spin, i) += mud(ic, ir).real();
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
                    if(LM->ParaV->in_this_processor(k1, k1))
                    {
                        const int ir = LM->ParaV->global2local_row(k1);
                        const int ic = LM->ParaV->global2local_col(k1);
                        MecMulP(0, j) += mud(ic, ir).real();
                        MecMulP(3, j) += mud(ic, ir).real();
                    }
                    if(LM->ParaV->in_this_processor(k1, k2))
                    {
                        const int ir = LM->ParaV->global2local_row(k1);
                        const int ic = LM->ParaV->global2local_col(k2);
                        MecMulP(1, j) += mud(ic, ir).real();
                        MecMulP(2, j) += mud(ic, ir).imag();
                    }
                    if(LM->ParaV->in_this_processor(k2, k1))
                    {
                        const int ir = LM->ParaV->global2local_row(k2);
                        const int ic = LM->ParaV->global2local_col(k1);
                        MecMulP(1, j) += mud(ic, ir).real();
                        MecMulP(2, j) -= mud(ic, ir).imag();
                    }
                    if(LM->ParaV->in_this_processor(k2, k2))
                    {
                        const int ir = LM->ParaV->global2local_row(k2);
                        const int ic = LM->ParaV->global2local_col(k2);
                        MecMulP(0, j) += mud(ic, ir).real();
                        MecMulP(3, j) -= mud(ic, ir).real();
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

std::vector<std::vector<std::vector<double>>> ModuleIO::convert(const ModuleBase::matrix &orbMulP)
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

template<typename T>
void ModuleIO::out_mulliken(const int& step, LCAO_Matrix* LM, const elecstate::ElecState* pelec, const K_Vectors& kv, hamilt::Hamilt<T>* ham_in)
{
    ModuleBase::TITLE("Mulliken_Charge", "out_mulliken");

    ModuleBase::matrix orbMulP;
    const std::vector<std::vector<T>>& dm =
        dynamic_cast<const elecstate::ElecStateLCAO<T>*> (pelec)->get_DM()->get_DMK_vector();
    orbMulP = ModuleIO::cal_mulliken(dm, LM, kv, ham_in);

    std::vector<std::vector<std::vector<double>>> AorbMulP = ModuleIO::convert(orbMulP);

    if(GlobalV::MY_RANK == 0)
    {
        const int nlocal = (GlobalV::NSPIN == 4) ? GlobalV::NLOCAL/2 : GlobalV::NLOCAL;
        std::stringstream as;
        as << GlobalV::global_out_dir << "mulliken.txt";
        std::ofstream os;
        if (step == 0)
        {
            os.open(as.str().c_str());
        }
        else
        {
            os.open(as.str().c_str(), std::ios::app);
        }
        os << "STEP: " << step << std::endl;
        os << "CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM" << std::endl;

		double sch = 0.0;
		os << std::setprecision(4);
		for(size_t is=0; is!=GlobalV::NSPIN; ++is)
		{
            if(GlobalV::NSPIN == 4 && is>0) continue;
			double sss = 0.0;
			for(size_t iw=0; iw!=nlocal; ++iw)
			{
				sch += orbMulP(is, iw);
				sss += orbMulP(is, iw);
			}
            if(GlobalV::NSPIN == 2)
            {
			    os << " Total charge of spin " << is+1 << ":\t" << sss << std::endl;
            }
        }
		os << " Total charge:\t" << sch << std::endl;
		os << "Decomposed Mulliken populations" << std::endl;
        GlobalV::ofs_running << std::endl <<std::endl;

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
                            double spin1 = ModuleIO::output_cut(AorbMulP[0][i][num]);
                            os << ModuleBase::Name_Angular[L][M] << std::setw(25) << Z << std::setw(32) << spin1 << std::endl;
                            sum_m[0] += AorbMulP[0][i][num];
                        }
                        else if (GlobalV::NSPIN==2)
                        {
                            double spin1 = ModuleIO::output_cut(AorbMulP[0][i][num]); 
                            double spin2 = ModuleIO::output_cut(AorbMulP[1][i][num]);
                            double sum = ModuleIO::output_cut(spin1 + spin2);
                            double diff = ModuleIO::output_cut(spin1 - spin2);
                            os << ModuleBase::Name_Angular[L][M] << std::setw(25) << Z << std::setw(32) << spin1 << std::setw(30) << spin2 << std::setw(30) << sum << std::setw(30) << diff << std::endl;
                            sum_m[0] += AorbMulP[0][i][num];
                            sum_m[1] += AorbMulP[1][i][num];
                        }
                        else if (GlobalV::NSPIN==4)
                        {
                            double spin[4];
                            for(int j=0;j<4;j++)
                            {
                                spin[j] = ModuleIO::output_cut(AorbMulP[j][i][num]);
                                sum_m[j] += AorbMulP[j][i][num];
                            } 
                            os << ModuleBase::Name_Angular[L][M] << std::setw(25) << Z << std::setw(32) << spin[0] << std::setw(30) << spin[1] << std::setw(30) << spin[2] << std::setw(30) << spin[3] << std::endl;
                        }
                        num++;
                    }

                    if (GlobalV::NSPIN==1)
                    {
                        double spin1 = ModuleIO::output_cut(sum_m[0]);
                        os << "  sum over m "<< std::setw(45) << spin1 << std::endl;
                        sum_l[0] += sum_m[0];
                    }
                    else if (GlobalV::NSPIN==2)
                    {
                        double spin1 = ModuleIO::output_cut(sum_m[0]);
                        double spin2 = ModuleIO::output_cut(sum_m[1]);
                        double sum = ModuleIO::output_cut(spin1 + spin2);
                        double diff = ModuleIO::output_cut(spin1 - spin2);
                        os << "  sum over m "<< std::setw(45) << spin1 << std::setw(30) << spin2 << std::setw(35) << sum << std::setw(25) << diff << std::endl;
                        sum_l[0] += sum_m[0];
                        sum_l[1] += sum_m[1];
                    }
                    else if (GlobalV::NSPIN==4)
                    {
                        double spin[4];
                        for(int j=0;j<4;j++)
                        {
                            spin[j] = ModuleIO::output_cut(sum_m[j]);
                            sum_l[j] += sum_m[j];
                        }
                        os << "  sum over m "<< std::setw(45) << spin[0] << std::setw(30) << spin[1] << std::setw(30) << spin[2] << std::setw(30) << spin[3] << std::endl;
                    }
                }
            
                if(GlobalC::ucell.atoms[t].l_nchi[L])
                {
                    if (GlobalV::NSPIN==1)
                    {
                        double spin1 = ModuleIO::output_cut(sum_l[0]);
                        os << "  sum over m+zeta "<< std::setw(40) << spin1 << std::endl;
                        total_charge += sum_l[0];
                    }
                    else if (GlobalV::NSPIN==2)
                    {
                        double spin1 = ModuleIO::output_cut(sum_l[0]);
                        double spin2 = ModuleIO::output_cut(sum_l[1]);
                        double sum = ModuleIO::output_cut(spin1 + spin2);
                        double diff = ModuleIO::output_cut(spin1 - spin2);
                        os << "  sum over m+zeta "<< std::setw(40) << spin1 << std::setw(30) << spin2 << std::setw(35) << sum << std::setw(25) << diff << std::endl;
                        total_charge += sum_l[0] + sum_l[1];
                        atom_mag += sum_l[0] - sum_l[1];
                    }
                    else if (GlobalV::NSPIN==4)
                    {
                        double spin[4];
                        for(int j=0;j<4;j++)
                        {
                            spin[j] = ModuleIO::output_cut(sum_l[j]);
                            total_charge_soc[j] += sum_l[j];
                        }
                        os << "  sum over m+zeta "<< std::setw(40) << spin[0] << std::setw(30) << spin[1] << std::setw(30) << spin[2] << std::setw(30) << spin[3] << std::endl;
                    }
                }
            }

            if (GlobalV::NSPIN==1)
                os << "Total Charge on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(20) << total_charge <<std::endl;
            else if (GlobalV::NSPIN==2)
            {
                os << "Total Charge on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(20) << total_charge <<std::endl;
                os << "Total Magnetism on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(20) << ModuleIO::output_cut(atom_mag) <<std::endl;
                GlobalV::ofs_running << "Total Magnetism on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(20) << std::setprecision(10) << atom_mag <<std::endl;
            }
            else if (GlobalV::NSPIN==4)
            {
                double spin1 = ModuleIO::output_cut(total_charge_soc[0]);
                double spin2 = ModuleIO::output_cut(total_charge_soc[1]);
                double spin3 = ModuleIO::output_cut(total_charge_soc[2]);
                double spin4 = ModuleIO::output_cut(total_charge_soc[3]);
                os << "Total Charge on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(20) 
                << spin1 <<std::endl;
                os << "Total Magnetism on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(20) 
                << "("  << spin2 << ", " << spin3 << ", " << spin4 << ")" 
                <<std::endl;
                GlobalV::ofs_running << "Total Magnetism on atom:  " << GlobalC::ucell.atoms[t].label <<  std::setw(10)
                << "("  << std::setprecision(10) << spin2 << ", " << std::setprecision(10) << spin3 << ", " << std::setprecision(10) << spin4 << ")"
                <<std::endl;
            }
            os << std::endl <<std::endl;
        }
        GlobalV::ofs_running << std::endl <<std::endl;
        os.close();
        ModuleIO::write_orb_info(&(GlobalC::ucell));
    }
}
template void ModuleIO::out_mulliken(const int& step, LCAO_Matrix* LM, const elecstate::ElecState* pelec, const K_Vectors& kv, hamilt::Hamilt<double>* ham_in);
template void ModuleIO::out_mulliken(const int& step, LCAO_Matrix* LM, const elecstate::ElecState* pelec, const K_Vectors& kv, hamilt::Hamilt<std::complex<double>>* ham_in);