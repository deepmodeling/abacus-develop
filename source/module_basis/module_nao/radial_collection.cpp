#include "module_basis/module_nao/radial_collection.h"
#include <memory>

#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/atomic_radials.h"
#include "module_basis/module_nao/beta_radials.h"
#include "module_basis/module_nao/sphbes_radials.h"

#include "module_base/parallel_common.h"
#include "module_base/tool_quit.h"
#include "module_base/global_variable.h"

RadialCollection::RadialCollection(const RadialCollection& other) :
    ntype_(other.ntype_),
    lmax_(other.lmax_),
    nchi_(other.nchi_),
    nzeta_max_(other.nzeta_max_),
    rcut_max_(other.rcut_max_),
    radset_(nullptr),
    iter_(nullptr),
    nl_(nullptr)
{
    if (ntype_ == 0)
    {
        return;
    }

    radset_ = new RadialSet*[ntype_];
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = other.radset_[itype]->clone();
    }

    iter_build();
}

RadialCollection& RadialCollection::operator=(const RadialCollection& rhs)
{
    if (&rhs == this)
    {
        return *this;
    }

    cleanup();

    ntype_ = rhs.ntype_;
    lmax_ = rhs.lmax_;
    nchi_ = rhs.nchi_;
    nzeta_max_ = rhs.nzeta_max_;
    rcut_max_ = rhs.rcut_max_;

    radset_ = new RadialSet*[ntype_];
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = rhs.radset_[itype]->clone();
    }

    iter_build();

    return *this;
}

RadialCollection::~RadialCollection()
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        delete radset_[itype];
    }
    delete[] radset_;
    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    delete[] nl_;
}

void RadialCollection::set_rcut_max()
{
    rcut_max_ = 0.0;
    for (int itype = 0; itype < ntype_; ++itype)
    {
        rcut_max_ = std::max(rcut_max_, radset_[itype]->rcut_max());
    }
}

void RadialCollection::cleanup()
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        delete radset_[itype];
    }
    delete[] radset_;
    radset_ = nullptr;

    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    iter_ = nullptr;

    delete[] nl_;
    nl_ = nullptr;

    ntype_ = 0;
    lmax_ = -1;
    nchi_ = 0;
    nzeta_max_ = 0;
}

void RadialCollection::iter_build()
{
    /*
     * collect the addresses of NumericalRadial objects from different RadialSet objects
     * so that all NumericalRadial objects can be iterated over in a single loop
     *
     * objects are sorted by l first, by itype next, by izeta last.
     *                                                                                      */
    delete[] iter_; // iterator does not control memory; simply delete the pointer array
    delete[] nl_;

    nl_ = new int[lmax_ + 1];
    iter_ = new const NumericalRadial*[nchi_];

    int i = 0;
    std::fill(nl_, nl_ + lmax_ + 1, 0);
    for (int l = 0; l <= lmax_; ++l)
    {
        for (int itype = 0; itype != ntype_; ++itype)
        {
            for (int izeta = 0; izeta < radset_[itype]->nzeta(l); ++izeta)
            {
                iter_[i] = &radset_[itype]->chi(l, izeta);
                ++i;
                ++nl_[l];
            }
        }
    }
}

void RadialCollection::build(const int ntype, Numerical_Nonlocal* const nls)
{
    cleanup();
    ntype_ = ntype;
    radset_ = new RadialSet*[ntype_];

    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype] = new BetaRadials;
        radset_[itype]->build(nls[itype], itype);

        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    iter_build();
    set_rcut_max();
}

void RadialCollection::build(const int nfile, const std::string* const file, const char ftype)
{

    cleanup();

    ntype_ = nfile;
#ifdef __MPI
    Parallel_Common::bcast_int(ntype_);
#endif

    radset_ = new RadialSet*[ntype_];
    char* file_type = new char[ntype_];

    if (ftype)
    { // simply use the given file type if given
        std::fill(file_type, file_type + ntype_, ftype);
    }
    else
    { // otherwise check the file type
        for (int itype = 0; itype < ntype_; ++itype)
        {
            file_type[itype] = check_file_type(file[itype]);
        }
    }

    for (int itype = 0; itype < ntype_; ++itype)
    {
        switch(file_type[itype])
        {
          case 'o': // orbital file
            radset_[itype] = new AtomicRadials;
            break;
          case 'c': // coefficient file
            radset_[itype] = new SphbesRadials;
            break;
          default: // not supposed to happend
            ModuleBase::WARNING_QUIT("RadialCollection::build", "Unrecognized file: " + file[itype]);
        }
        radset_[itype]->build(file[itype], itype);
    }

    for (int itype = 0; itype < ntype_; ++itype)
    {
        lmax_ = std::max(lmax_, radset_[itype]->lmax());
        nchi_ += radset_[itype]->nchi();
        nzeta_max_ = std::max(nzeta_max_, radset_[itype]->nzeta_max());
    }

    iter_build();
    set_rcut_max();
}

void RadialCollection::set_transformer(ModuleBase::SphericalBesselTransformer sbt, const int update)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_transformer(sbt, update);
    }
}

void RadialCollection::set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_grid(for_r_space, ngrid, grid, mode);
    }
    if (for_r_space) rcut_max_ = grid[ngrid - 1];
}

void RadialCollection::set_uniform_grid(const bool for_r_space,
                                        const int ngrid,
                                        const double cutoff,
                                        const char mode,
                                        const bool enable_fft)
{
    for (int itype = 0; itype < ntype_; ++itype)
    {
        radset_[itype]->set_uniform_grid(for_r_space, ngrid, cutoff, mode, enable_fft);
    }

    if (for_r_space) rcut_max_ = cutoff;
}

char RadialCollection::check_file_type(const std::string& file) const
{
    // currently we only support ABACUS numerical atomic orbital file and
    // SIAB/PTG-generated orbital coefficient file. The latter contains a
    // <Coefficients ...> block, which is not present in the former.
    //
    // Unfortunately, the numerial atomic orbital file does not have any
    // distinguishing feature. Many keywords in the orbital file may also
    // be found in the coefficient file. Here we simply assume that if the
    // file contains a <Coefficients ...> block, it is a coefficient file;
    // otherwise it is an orbital file.

    char file_type = 'o';
    if (GlobalV::MY_RANK == 0)
    {
        std::ifstream ifs(file.c_str());
        std::string line;
        while (std::getline(ifs, line))
        {
            if (line.find("<Coefficients") != std::string::npos)
            {
                file_type = 'c';
                break;
            }
        }
        ifs.close();
    }
#ifdef __MPI
    Parallel_Common::bcast_char(&file_type, 1);
#endif
    return file_type;
}
