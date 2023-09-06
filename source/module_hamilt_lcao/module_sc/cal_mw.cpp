#include "spin_constrain.h"

#include <iostream>

#include "module_base/matrix.h"
#include "module_base/name_angular.h"
#include "module_base/tool_title.h"
#include "module_base/scalapack_connector.h"

void SpinConstrain::cal_MW(const int& step,
                        LCAO_Hamilt& uhm,
                        const std::vector<ModuleBase::ComplexMatrix> &dm,
                        const K_Vectors& kv,
                        const UnitCell& ucell)
{
    ModuleBase::TITLE("module_sc", "cal_MW");

    ModuleBase::matrix orbMulP;
    orbMulP = this->cal_MW_k(dm, uhm, kv);
    
    std::vector<std::vector<std::vector<double>>> AorbMulP = this->convert(orbMulP);

    size_t nw = this->get_nw();
    
    if (GlobalV::MY_RANK == 0)
    {
        const int nlocal = nw / 2;
        std::stringstream as;
        as << GlobalV::global_out_dir << "atomic_magnetization.txt";
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
        os << "CALCULATE THE MULLIkEN ANALYSIS OF MAGNETIZATION FOR EACH ATOM" << std::endl;

        double sch = 0.0;
        os << std::setprecision(4);
        for (int is = 0; is < this->nspin_; ++is)
        {
            if (this->nspin_ == 4 && is > 0)
                continue;
            double sss = 0.0;
            for (size_t iw = 0; iw != nlocal; ++iw)
            {
                sch += orbMulP(is, iw);
                sss += orbMulP(is, iw);
            }
        }
        os << " Total charge:\t" << sch << std::endl;
        os << "Decomposed Mulliken populations" << std::endl;

        for (size_t i = 0; i != ucell.nat; ++i)
        {
            double total_charge = 0.0, atom_mag = 0.0;
            std::vector<double> total_charge_soc(this->nspin_);
            const int t = ucell.iat2it[i];
            int num = 0;
            os << i << std::setw(25) << "Zeta of " << ucell.atoms[t].label << std::setw(30) << "Spin 1" << std::setw(30)
               << "Spin 2" << std::setw(30) << "Spin 3" << std::setw(30) << "Spin 4" << std::endl;

            for (size_t L = 0; L != ucell.atoms[t].nwl + 1; ++L)
            {
                std::vector<double> sum_l(this->nspin_, 0.0);
                for (size_t Z = 0; Z != ucell.atoms[t].l_nchi[L]; ++Z)
                {
                    std::vector<double> sum_m(this->nspin_, 0.0);
                    for (size_t M = 0; M != (2 * L + 1); ++M)
                    {
                        double spin[4];
                        for (int j = 0; j < 4; j++)
                        {
                            spin[j] = this->output_cut(AorbMulP[j][i][num]);
                            sum_m[j] += AorbMulP[j][i][num];
                        }
                        os << ModuleBase::Name_Angular[L][M] << std::setw(25) << Z << std::setw(32) << spin[0]
                           << std::setw(30) << spin[1] << std::setw(30) << spin[2] << std::setw(30) << spin[3]
                           << std::endl;
                        num++;
                    }

                    double spin[4];
                    for (int j = 0; j < 4; j++)
                    {
                        spin[j] = this->output_cut(sum_m[j]);
                        sum_l[j] += sum_m[j];
                    }
                    os << "  sum over m " << std::setw(45) << spin[0] << std::setw(30) << spin[1] << std::setw(30)
                       << spin[2] << std::setw(30) << spin[3] << std::endl;
                }

                if (ucell.atoms[t].l_nchi[L])
                {
                    double spin[4];
                    for (int j = 0; j < 4; j++)
                    {
                        spin[j] = this->output_cut(sum_l[j]);
                        total_charge_soc[j] += sum_l[j];
                    }
                    os << "  sum over m+zeta " << std::setw(40) << spin[0] << std::setw(30) << spin[1] << std::setw(30)
                       << spin[2] << std::setw(30) << spin[3] << std::endl;
                }
            }

            double spin1 = this->output_cut(total_charge_soc[0]);
            double spin2 = this->output_cut(total_charge_soc[1]);
            double spin3 = this->output_cut(total_charge_soc[2]);
            double spin4 = this->output_cut(total_charge_soc[3]);
            os << "Total Charge on atom:  " << ucell.atoms[t].label << std::setw(20) << spin1 << std::endl;
            os << "Total Magnetism on atom:  " << ucell.atoms[t].label << std::setw(20) << "(" << spin2 << ", " << spin3
               << ", " << spin4 << ")" << std::endl;
            os << std::endl << std::endl;
        }
        os.close();
    }
}

ModuleBase::matrix SpinConstrain::cal_MW_k(const std::vector<ModuleBase::ComplexMatrix> &dm,
    LCAO_Hamilt &uhm, const K_Vectors& kv
)
{
    ModuleBase::TITLE("module_sc", "cal_MW_k");
    int nw = this->get_nw();
    const int nlocal = nw/2;
    ModuleBase::matrix MecMulP, orbMulP;
    MecMulP.create(this->nspin_, nlocal);
    orbMulP.create(this->nspin_, nlocal);

    for(size_t ik = 0; ik != kv.nks; ++ik)
    {
        uhm.LM->zeros_HSk('S');
		uhm.LM->folding_fixedH(ik, kv.kvec_d);

        ModuleBase::ComplexMatrix mud;
        mud.create(this->ParaV->ncol, this->ParaV->nrow);

#ifdef __MPI
        const char T_char = 'T';
        const char N_char = 'N';
        const int one_int = 1;
        const std::complex<double> one_float = {1.0, 0.0}, zero_float = {0.0, 0.0};        
        pzgemm_(&T_char,
                &T_char,
                &nw,
                &nw,
                &nw,
                &one_float,
                dm[ik].c,
                &one_int,
                &one_int,
                this->ParaV->desc,
                uhm.LM->Sloc2.data(),
                &one_int,
                &one_int,
                this->ParaV->desc,
                &zero_float,
                mud.c,
                &one_int,
                &one_int,
                this->ParaV->desc);

        for(size_t i=0; i < nw; ++i)
        {
            const int index = i%2;
            if(!index)
            {
                const int j = i/2;
                const int k1 = 2*j;
                const int k2 = 2*j+1;
                if(this->ParaV->in_this_processor(k1, k1))
                {
                    const int ir = this->ParaV->global2local_row(k1);
                    const int ic = this->ParaV->global2local_col(k1);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) += mud(ic, ir).real();
                }
                if(this->ParaV->in_this_processor(k1, k2))
                {
                    const int ir = this->ParaV->global2local_row(k1);
                    const int ic = this->ParaV->global2local_col(k2);
                    MecMulP(1, j) += mud(ic, ir).real();
                    MecMulP(2, j) += mud(ic, ir).imag();
                }
                if(this->ParaV->in_this_processor(k2, k1))
                {
                    const int ir = this->ParaV->global2local_row(k2);
                    const int ic = this->ParaV->global2local_col(k1);
                    MecMulP(1, j) += mud(ic, ir).real();
                    MecMulP(2, j) -= mud(ic, ir).imag();
                }
                if(this->ParaV->in_this_processor(k2, k2))
                {
                    const int ir = this->ParaV->global2local_row(k2);
                    const int ic = this->ParaV->global2local_col(k2);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) -= mud(ic, ir).real();
                }
            }
        }
#endif
    }
#ifdef __MPI
    MPI_Reduce(MecMulP.c, orbMulP.c, this->nspin_*nlocal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif 

    return orbMulP;
}

std::vector<std::vector<std::vector<double>>> SpinConstrain::convert(const ModuleBase::matrix &orbMulP)
{
    std::vector<std::vector<std::vector<double>>> AorbMulP;
    AorbMulP.resize(this->nspin_);
    int nat = this->get_nat();
    for (int is=0; is < this->nspin_; ++is)
    {
        int num=0;
        AorbMulP[is].resize(nat);
        for (const auto& sc_elem : this->get_atomCounts())
        {
            int it = sc_elem.first;
            int nat_it = sc_elem.second;
            int nw_it = this->get_orbitalCounts().at(it);
            for (int ia = 0; ia < nat_it; ia++)
            {
                int iat = this->get_iat(it, ia);
                AorbMulP[is][iat].resize(nw_it);
                for (int iw = 0; iw < nw_it; iw++)
                {
                    AorbMulP[is][iat][iw] = orbMulP(is, num);
                    num++;
                }
            }
        }
    }
    
    return AorbMulP;
}