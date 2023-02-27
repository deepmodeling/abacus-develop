/**********************************************************************
  Mulliken_Charge.cpp:

     Mulliken_Charge.cpp is a subrutine to calculate Mulliken charge.
 
  Log of Mulliken_Charge.cpp:

     12/Oct/2018  Released by Feng Qi

***********************************************************************/

#include "mulliken_charge.h"
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
        const char N_char = 'N';
        const int one_int = 1;
        const double one_float[2] = {1.0, 0.0}, zero_float[2] = {0.0, 0.0};        
        pdgemm_(&N_char,
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
                        MecMulP(0, i) += mud(ir, ic);
                    }
                    else if(uhm.LM->ParaV->in_this_processor(k1, k2))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k1];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k2];
                        MecMulP(1, i) += mud(ir, ic);
                    }
                    else if(uhm.LM->ParaV->in_this_processor(k2, k1))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k2];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k1];
                        MecMulP(2, i) += mud(ir, ic);
                    }
                    else if(uhm.LM->ParaV->in_this_processor(k2, k2))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k2];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k2];
                        MecMulP(3, i) += mud(ir, ic);
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
        ModuleBase::ComplexMatrix mud;
        mud.create(uhm.LM->ParaV->nrow, uhm.LM->ParaV->ncol);

#ifdef __MPI
        const char N_char = 'N';
        const int one_int = 1;
        const double one_float[2] = {1.0, 0.0}, zero_float[2] = {0.0, 0.0};        
        pzgemm_(&N_char,
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
                        MecMulP(0, i) += mud(ir, ic).real();
                    }
                    else if(uhm.LM->ParaV->in_this_processor(k1, k2))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k1];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k2];
                        MecMulP(1, i) += mud(ir, ic).real();
                    }
                    else if(uhm.LM->ParaV->in_this_processor(k2, k1))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k2];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k1];
                        MecMulP(2, i) += mud(ir, ic).real();
                    }
                    else if(uhm.LM->ParaV->in_this_processor(k2, k2))
                    {
                        const int ir = uhm.LM->ParaV->trace_loc_row[k2];
                        const int ic = uhm.LM->ParaV->trace_loc_col[k2];
                        MecMulP(3, i) += mud(ir, ic).real();
                    }
                }
            }
        }
    atom_arrange::delete_vector(
		GlobalV::ofs_running,
		GlobalV::SEARCH_PBC,
		GlobalC::GridD, 
		GlobalC::ucell, 
		GlobalV::SEARCH_RADIUS, 
		GlobalV::test_atom_input);
#endif
    }
#ifdef __MPI 
    MPI_Reduce(MecMulP.c, orbMulP.c, GlobalV::NSPIN*nlocal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif 

    return orbMulP;
}

void Mulliken_Charge::out_mulliken(LCAO_Hamilt &uhm, Local_Orbital_Charge &loc)
{
    ModuleBase::matrix orbMulP;
    if(GlobalV::GAMMA_ONLY_LOCAL)
        orbMulP = this->cal_mulliken(loc.dm_gamma, uhm);
    else
        orbMulP = this->cal_mulliken_k(loc.dm_k, uhm);

    if(GlobalV::MY_RANK == 0)
    {
        const int nlocal = (GlobalV::NSPIN == 4) ? GlobalV::NLOCAL/2 : GlobalV::NLOCAL;
        std::stringstream as;
        as << GlobalV::global_out_dir << "Mulliken.dat";
        std::ofstream os(as.str().c_str());
        for(size_t i = 0; i != GlobalV::NLOCAL; ++i)
        {
            if(GlobalV::NSPIN == 1)
                os << std::setw(13) << orbMulP(0, i) << std::endl;
            else if(GlobalV::NSPIN == 2)
                os << std::setw(13) << orbMulP(0, i) << std::setw(13) << orbMulP(1, i) << std::endl;
            else if(GlobalV::NSPIN == 4)
                os << std::setw(13) << orbMulP(0, i) << std::setw(13) << orbMulP(1, i) << std::setw(13) << orbMulP(2, i) << std::setw(13) << orbMulP(3, i) << std::endl;
        }
        os.close();
        ModuleIO::write_orb_info();
    }
}