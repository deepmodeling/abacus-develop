#include "density_matrix.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/libm/libm.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/klist.h"

#include <cstddef>
#include <memory>
#include <stdexcept>

namespace module_dm
{

// initialize density matrix DMR from UnitCell (mainly used in UnitTest)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_dmr(const Grid_Driver* GridD_in, const UnitCell* ucell)
{
    ModuleBase::TITLE("DensityMatrix", "init_dmr");
    this->clear_dmr();
    // construct a new DMR
    std::unique_ptr<hamilt::HContainer<TR>> tmp_DMR(new hamilt::HContainer<TR>(this->pv));
    // set up a HContainer
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        ModuleBase::Vector3<double> tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD_in->Find_atom(*ucell, tau1, T1, I1, &adjs);
        // std::cout << "adjs.adj_num: " <<adjs.adj_num << std::endl;
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell->itia2iat(T2, I2);
            if (this->pv->is_invalid_atom_pair(iat1, iat2))
            {
                continue;
            }
            ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            // std::cout << "R_index: " << R_index.x << " " << R_index.y << " " << R_index.z << std::endl;
            hamilt::AtomPair<TR> tmp_ap(iat1, iat2, R_index, this->pv);
            tmp_DMR->insert_pair(tmp_ap);
        }
    }
    // allocate the memory of BaseMatrix in SR, and set the new values to zero
    if (std::is_same<TK, double>::value)
    {
        tmp_DMR->fix_gamma();
    }
    tmp_DMR->allocate(nullptr, true);
    this->dmr.push_back(tmp_DMR.release());
    // add another DMR if nspin==2
    if (this->spin_mult == 2)
    {
        std::unique_ptr<hamilt::HContainer<TR>> tmp_DMR1(new hamilt::HContainer<TR>(*this->dmr[0]));
        this->dmr.push_back(tmp_DMR1.release());
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->dmr.size() * this->dmr[0]->get_memory_size());
}

/// initialize density matrix DMR from UnitCell and RA (mainly used in UnitTest)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_dmr(Record_adj& ra, const UnitCell* ucell)
{
    ModuleBase::TITLE("DensityMatrix", "init_dmr");
    this->clear_dmr();
    // construct a new DMR
    std::unique_ptr<hamilt::HContainer<TR>> tmp_DMR(new hamilt::HContainer<TR>(this->pv));
    // set up a HContainer
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        ModuleBase::Vector3<double> tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        for (int ad = 0; ad < ra.na_each[iat1]; ++ad)
        {
            const int T2 = ra.get_info(iat1, ad)[3];
            const int I2 = ra.get_info(iat1, ad)[4];
            int iat2 = ucell->itia2iat(T2, I2);
            if (this->pv->is_invalid_atom_pair(iat1, iat2))
            {
                continue;
            }
            hamilt::AtomPair<TR> tmp_ap(iat1,
                                        iat2,
                                        ra.get_info(iat1, ad)[0],
                                        ra.get_info(iat1, ad)[1],
                                        ra.get_info(iat1, ad)[2],
                                        this->pv);
            tmp_DMR->insert_pair(tmp_ap);
        }
    }
    if (std::is_same<TK, double>::value)
    {
        tmp_DMR->fix_gamma();
    }
    tmp_DMR->allocate(nullptr, true);
    this->dmr.push_back(tmp_DMR.release());
    // add another DMR if nspin==2
    if (this->spin_mult == 2)
    {
        std::unique_ptr<hamilt::HContainer<TR>> tmp_DMR1(new hamilt::HContainer<TR>(*this->dmr[0]));
        this->dmr.push_back(tmp_DMR1.release());
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->dmr.size() * this->dmr[0]->get_memory_size());
}

// initialize density matrix DMR from another HContainer (mainly used)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_dmr(const hamilt::HContainer<TR>& DMR_in)
{
    ModuleBase::TITLE("DensityMatrix", "init_dmr");
    this->clear_dmr();
    // set up a HContainer using another one
    for (int is = 0; is < this->spin_mult; ++is) // loop over spin
    {
        std::unique_ptr<hamilt::HContainer<TR>> tmp_DMR(new hamilt::HContainer<TR>(DMR_in));
        // zero.out
        tmp_DMR->set_zero();
        this->dmr.push_back(tmp_DMR.release());
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->dmr.size() * this->dmr[0]->get_memory_size());
}

template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_dmr(const hamilt::HContainer<TRShift>& DMR_in)
{
    ModuleBase::TITLE("DensityMatrix", "init_dmr");
    this->clear_dmr();
    // set up a HContainer using another one
    int size_ap = DMR_in.size_atom_pairs();
    if (size_ap > 0)
    {
        const Parallel_Orbitals* paraV_ = DMR_in.get_atom_pair(0).get_paraV();
        std::unique_ptr<hamilt::HContainer<TR>> tmp_DMR(new hamilt::HContainer<TR>(paraV_));
        for (int iap = 0; iap < size_ap; iap++)
        {
            const int iat1 = DMR_in.get_atom_pair(iap).get_atom_i();
            const int iat2 = DMR_in.get_atom_pair(iap).get_atom_j();
            for (int ir = 0; ir < DMR_in.get_atom_pair(iap).get_R_size(); ir++)
            {
                const ModuleBase::Vector3<int> R_index = DMR_in.get_atom_pair(iap).get_R_index(ir);
                hamilt::AtomPair<TR> tmp_ap(iat1, iat2, R_index, paraV_);
                tmp_DMR->insert_pair(tmp_ap);
            }
        }
        tmp_DMR->allocate(nullptr, true);
        this->dmr.push_back(tmp_DMR.release());
        if (this->spin_mult == 2)
        {
            std::unique_ptr<hamilt::HContainer<TR>> tmp_DMR1(new hamilt::HContainer<TR>(*this->dmr[0]));
            this->dmr.push_back(tmp_DMR1.release());
        }
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->dmr.size() * this->dmr[0]->get_memory_size());
}

// T of HContainer can be double or std::complex<double>
template class DensityMatrix<double, double>;               // Gamma-Only case
template class DensityMatrix<std::complex<double>, double>; // Multi-k case
template class DensityMatrix<std::complex<double>, std::complex<double>>; // For EXX in future

} // namespace module_dm
