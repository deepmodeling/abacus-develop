#include "density_matrix.h"

#include "module_parameter/parameter.h"
#include "module_base/libm/libm.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_cell/klist.h"

namespace elecstate
{

// initialize density matrix DMR from UnitCell (mainly used in UnitTest)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(Grid_Driver* GridD_in, const UnitCell* ucell)
{
    ModuleBase::TITLE("DensityMatrix", "init_DMR");
    // ensure _DMR is empty
    for (auto& it: this->_DMR)
    {
        delete it;
    }
    this->_DMR.clear();
    // construct a new DMR
    hamilt::HContainer<TR>* tmp_DMR;
    tmp_DMR = new hamilt::HContainer<TR>(this->_paraV);
    // set up a HContainer
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
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
            if (this->_paraV->get_row_size(iat1) <= 0 || this->_paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            // std::cout << "R_index: " << R_index.x << " " << R_index.y << " " << R_index.z << std::endl;
            hamilt::AtomPair<TR> tmp_ap(iat1, iat2, R_index, this->_paraV);
            tmp_DMR->insert_pair(tmp_ap);
        }
    }
    // allocate the memory of BaseMatrix in SR, and set the new values to zero
    if (std::is_same<TK, double>::value)
    {
        tmp_DMR->fix_gamma();
    }
    tmp_DMR->allocate(nullptr, true);
    this->_DMR.push_back(tmp_DMR);
    // add another DMR if nspin==2
    if (this->_nspin == 2)
    {
        hamilt::HContainer<TR>* tmp_DMR1;
        tmp_DMR1 = new hamilt::HContainer<TR>(*tmp_DMR);
        this->_DMR.push_back(tmp_DMR1);
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->_DMR.size() * this->_DMR[0]->get_memory_size());
}

/// initialize density matrix DMR from UnitCell and RA (mainly used in UnitTest)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(Record_adj& ra, const UnitCell* ucell)
{
    ModuleBase::TITLE("DensityMatrix", "init_DMR");
    // ensure _DMR is empty
    for (auto& it: this->_DMR)
    {
        delete it;
    }
    this->_DMR.clear();
    // construct a new DMR
    hamilt::HContainer<TR>* tmp_DMR;
    tmp_DMR = new hamilt::HContainer<TR>(this->_paraV);
    // set up a HContainer
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        for (int ad = 0; ad < ra.na_each[iat1]; ++ad)
        {
            const int T2 = ra.info[iat1][ad][3];
            const int I2 = ra.info[iat1][ad][4];
            int iat2 = ucell->itia2iat(T2, I2);
            if (this->_paraV->get_row_size(iat1) <= 0 || this->_paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            hamilt::AtomPair<TR> tmp_ap(iat1,
                                        iat2,
                                        ra.info[iat1][ad][0],
                                        ra.info[iat1][ad][1],
                                        ra.info[iat1][ad][2],
                                        this->_paraV);
            tmp_DMR->insert_pair(tmp_ap);
        }
    }
    // allocate the memory of BaseMatrix in SR, and set the new values to zero
    if (std::is_same<TK, double>::value)
    {
        tmp_DMR->fix_gamma();
    }
    tmp_DMR->allocate(nullptr, true);
    this->_DMR.push_back(tmp_DMR);
    // add another DMR if nspin==2
    if (this->_nspin == 2)
    {
        hamilt::HContainer<TR>* tmp_DMR1;
        tmp_DMR1 = new hamilt::HContainer<TR>(*tmp_DMR);
        this->_DMR.push_back(tmp_DMR1);
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->_DMR.size() * this->_DMR[0]->get_memory_size());
}

// initialize density matrix DMR from another HContainer (mainly used)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(const hamilt::HContainer<TR>& DMR_in)
{
    ModuleBase::TITLE("DensityMatrix", "init_DMR");
    // ensure _DMR is empty
    for (auto& it: this->_DMR)
    {
        delete it;
    }
    this->_DMR.clear();
    // set up a HContainer using another one
    for (int is = 0; is < this->_nspin; ++is) // loop over spin
    {
        hamilt::HContainer<TR>* tmp_DMR;
        tmp_DMR = new hamilt::HContainer<TR>(DMR_in);
        // zero.out
        tmp_DMR->set_zero();
        this->_DMR.push_back(tmp_DMR);
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->_DMR.size() * this->_DMR[0]->get_memory_size());
}

template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(const hamilt::HContainer<TRShift>& DMR_in)
{
    ModuleBase::TITLE("DensityMatrix", "init_DMR");
    // ensure _DMR is empty
    for (auto& it: this->_DMR)
    {
        delete it;
    }
    this->_DMR.clear();
    // set up a HContainer using another one
    int size_ap = DMR_in.size_atom_pairs();
    if (size_ap > 0)
    {
        const Parallel_Orbitals* paraV_ = DMR_in.get_atom_pair(0).get_paraV();
        hamilt::HContainer<TR>* tmp_DMR = new hamilt::HContainer<TR>(paraV_);
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
        this->_DMR.push_back(tmp_DMR);
        if (this->_nspin == 2)
        {
            hamilt::HContainer<TR>* tmp_DMR1 = new hamilt::HContainer<TR>(*tmp_DMR);
            this->_DMR.push_back(tmp_DMR1);
        }
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->_DMR.size() * this->_DMR[0]->get_memory_size());
}

// get _DMR pointer
template <typename TK, typename TR>
hamilt::HContainer<TR>* DensityMatrix<TK, TR>::get_DMR_pointer(const int ispin) const
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    return this->_DMR[ispin - 1];
}

// get _DMK[ik] pointer
template <typename TK, typename TR>
TK* DensityMatrix<TK, TR>::get_DMK_pointer(const int ik) const
{
#ifdef __DEBUG
    assert(ik < this->_nk * this->_nspin);
#endif
    return const_cast<TK*>(this->_DMK[ik].data());
}

// set DMK using a pointer
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK_pointer(const int ik, TK* DMK_in)
{
#ifdef __DEBUG
    assert(ik < this->_nk * this->_nspin);
#endif
    this->_DMK[ik].assign(DMK_in, DMK_in + this->_paraV->nrow * this->_paraV->ncol);
}

// set _DMK element
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK(const int ispin, const int ik, const int i, const int j, const TK value)
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
    assert(ik >= 0 && ik < this->_nk);
#endif
    // consider transpose col=>row
    this->_DMK[ik + this->_nk * (ispin - 1)][i * this->_paraV->nrow + j] = value;
}

// set _DMK element
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK_zero()
{
    for (int ik = 0; ik < _nspin * _nk; ik++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->_DMK[ik].data(),
                                      this->_paraV->get_row_size() * this->_paraV->get_col_size());
    }
}

// get a matrix element of density matrix dm(k)
template <typename TK, typename TR>
TK DensityMatrix<TK, TR>::get_DMK(const int ispin, const int ik, const int i, const int j) const
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    // consider transpose col=>row
    return this->_DMK[ik + this->_nk * (ispin - 1)][i * this->_paraV->nrow + j];
}

// get _DMK nks, nrow, ncol
template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_nks() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() == _nk * _nspin);
#endif
    return _nk * _nspin;
}

template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_size() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->_DMK.size();
}

template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_nrow() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->_paraV->nrow;
}

template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_ncol() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->_paraV->ncol;
}

template <typename TK, typename TR>
void DensityMatrix<TK, TR>::save_DMR()
{
    ModuleBase::TITLE("DensityMatrix", "save_DMR");
    ModuleBase::timer::tick("DensityMatrix", "save_DMR");

    const int nnr = this->_DMR[0]->get_nnr();
    // allocate if _DMR_save is empty
    if (_DMR_save.size() == 0)
    {
        _DMR_save.resize(this->_DMR.size());
    }
    // resize if _DMR_save[is].size is not equal to _DMR.size
    for (int is = 0; is < _DMR_save.size(); is++)
    {
        if (_DMR_save[is].size() != nnr)
        {
            _DMR_save[is].resize(nnr);
        }
    }
    // save _DMR to _DMR_save
    for (int is = 0; is < this->_DMR.size(); is++)
    {
        TR* DMR_pointer = this->_DMR[is]->get_wrapper();
        TR* DMR_save_pointer = _DMR_save[is].data();
        // set to zero
        ModuleBase::GlobalFunc::ZEROS(DMR_save_pointer, nnr);
        for (int i = 0; i < nnr; i++)
        {
            DMR_save_pointer[i] = DMR_pointer[i];
        }
    }

    ModuleBase::timer::tick("DensityMatrix", "save_DMR");
}

// read *.dmk into density matrix dm(k)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::read_DMK(const std::string directory, const int ispin, const int ik)
{
    ModuleBase::TITLE("DensityMatrix", "read_DMK");
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    // read
    std::string fn;
    fn = directory + "SPIN" + std::to_string(ispin) + "_" + std::to_string(ik) + ".dmk";
    //
    bool quit_abacus = false;

    std::ifstream ifs;

    ifs.open(fn.c_str());
    if (!ifs)
    {
        quit_abacus = true;
    }
    else
    {
        // if the number is not match,
        // quit the program or not.
        bool quit = false;

        ModuleBase::CHECK_DOUBLE(ifs, this->_kvec_d[ik].x, quit);
        ModuleBase::CHECK_DOUBLE(ifs, this->_kvec_d[ik].y, quit);
        ModuleBase::CHECK_DOUBLE(ifs, this->_kvec_d[ik].z, quit);
        ModuleBase::CHECK_INT(ifs, this->_paraV->nrow);
        ModuleBase::CHECK_INT(ifs, this->_paraV->ncol);
    } // If file exist, read in data.
    // Finish reading the first part of density matrix.

    for (int i = 0; i < this->_paraV->nrow; ++i)
    {
        for (int j = 0; j < this->_paraV->ncol; ++j)
        {
            ifs >> this->_DMK[ik + this->_nk * (ispin - 1)][i * this->_paraV->ncol + j];
        }
    }
    ifs.close();
}

// output density matrix dm(k) into *.dmk
template <>
void DensityMatrix<double, double>::write_DMK(const std::string directory, const int ispin, const int ik)
{
    ModuleBase::TITLE("DensityMatrix", "write_DMK");
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    // write
    std::string fn;
    fn = directory + "SPIN" + std::to_string(ispin) + "_" + std::to_string(ik) + ".dmk";
    std::ofstream ofs;
    ofs.open(fn.c_str());
    if (!ofs)
    {
        ModuleBase::WARNING("elecstate::write_dmk", "Can't create DENSITY MATRIX File!");
    }
    ofs << this->_kvec_d[ik].x << " " << this->_kvec_d[ik].y << " " << this->_kvec_d[ik].z << std::endl;
    ofs << "\n  " << this->_paraV->nrow << " " << this->_paraV->ncol << std::endl;

    ofs << std::setprecision(3);
    ofs << std::scientific;

    for (int i = 0; i < this->_paraV->nrow; ++i)
    {
        for (int j = 0; j < this->_paraV->ncol; ++j)
        {
            if (j % 8 == 0)
            {
                ofs << "\n";
            }
            ofs << " " << this->_DMK[ik + this->_nk * (ispin - 1)][i * this->_paraV->ncol + j];
        }
    }

    ofs.close();
}

template <>
void DensityMatrix<std::complex<double>, double>::write_DMK(const std::string directory, const int ispin, const int ik)
{
    ModuleBase::TITLE("DensityMatrix", "write_DMK");
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    // write
    std::string fn;
    fn = directory + "SPIN" + std::to_string(ispin) + "_" + std::to_string(ik) + ".dmk";
    std::ofstream ofs;
    ofs.open(fn.c_str());
    if (!ofs)
    {
        ModuleBase::WARNING("elecstate::write_dmk", "Can't create DENSITY MATRIX File!");
    }
    ofs << this->_kvec_d[ik].x << " " << this->_kvec_d[ik].y << " " << this->_kvec_d[ik].z << std::endl;
    ofs << "\n  " << this->_paraV->nrow << " " << this->_paraV->ncol << std::endl;

    ofs << std::setprecision(3);
    ofs << std::scientific;

    for (int i = 0; i < this->_paraV->nrow; ++i)
    {
        for (int j = 0; j < this->_paraV->ncol; ++j)
        {
            if (j % 8 == 0)
            {
                ofs << "\n";
            }
            ofs << " " << this->_DMK[ik + this->_nk * (ispin - 1)][i * this->_paraV->ncol + j].real();
        }
    }

    ofs.close();
}

// T of HContainer can be double or complex<double>
template class DensityMatrix<double, double>;               // Gamma-Only case
template class DensityMatrix<std::complex<double>, double>; // Multi-k case
template class DensityMatrix<std::complex<double>, std::complex<double>>; // For EXX in future

} // namespace elecstate