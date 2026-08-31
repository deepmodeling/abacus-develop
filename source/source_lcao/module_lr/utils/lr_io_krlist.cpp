#include "lr_io_krlist.h"
#include "lr_io.h"
#include "source_lcao/module_ri/ri_util.h"
#include "source_base/constants.h"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace LR_IO
{

namespace
{
const std::string FILE_COARSE = "stru_out";
const std::string FILE_FINE_UNIFORM = "band_kpath_info";
const std::string FILE_FINE_NONUNIFORM = "KPT_bse";
}

RI_kRlist::RI_kRlist(const UnitCell& ucell, K_Vectors* const pkv, const int nspin_in,
    const std::string& in_dir, const std::string& out_dir, const int use_fine_kgrid)
    : klist(pkv), nspin(nspin_in)
{
    read_kpts_coarse(in_dir + FILE_COARSE, ucell, this->klist, out_dir);
    this->klist_coarse = *this->klist;
    this->period = RI_Util::get_Born_vonKarmen_period(*klist);
    this->Rlist = RI_Util::get_Born_von_Karmen_cells(period);
    if (use_fine_kgrid == 1)
    {
        read_kpts_fine(in_dir + FILE_FINE_UNIFORM, ucell, this->klist, false, out_dir);
    }
    else if (use_fine_kgrid == 2)
    {
        read_kpts_fine(in_dir + FILE_FINE_NONUNIFORM, ucell, this->klist, true, out_dir);
    }
    else if (use_fine_kgrid != 0)
    {
        ModuleBase::WARNING_QUIT("LR_IO", "use_fine_kgrid must be 0, 1 or 2");
    }
}

void RI_kRlist::read_kpts_coarse(const std::string& file, const UnitCell& ucell,
                                 K_Vectors* const klist, const std::string& out_dir)
{
    std::ifstream ifs(file);
    if (!ifs) throw std::runtime_error(file + " not found");
    std::string tmp;
    for (int i = 0; i < 7; ++i) { std::getline(ifs, tmp); }
    const int nat = std::stoi(tmp);
    for (int i = 0; i != nat; ++i) { std::getline(ifs, tmp); }
    const int nks_original = klist->get_nks();

    ifs >> klist->nmp[0] >> klist->nmp[1] >> klist->nmp[2];
    const int nk = klist->nmp[0] * klist->nmp[1] * klist->nmp[2];
    const int nks = (this->nspin == 2) ? 2 * nk : nk;
    assert(nks == nks_original);

    for (int ik = 0; ik < nk; ++ik)
    {
        ifs >> klist->kvec_c[ik].x >> klist->kvec_c[ik].y >> klist->kvec_c[ik].z;
        klist->kvec_c[ik] *= (ucell.lat0 / ModuleBase::TWO_PI); // in unit of 2*pi/lat0
        klist->kvec_d[ik] = klist->kvec_c[ik] * ucell.latvec.Transpose();
        set_zero_if_close(klist->kvec_d[ik]);
        klist->wk[ik] = 1.0 / double(nk);
    }
    if (this->nspin == 2)
    {
        for (int ik = 0; ik < nk; ++ik)
        {
            klist->kvec_c[ik + nk] = klist->kvec_c[ik];
            klist->kvec_d[ik + nk] = klist->kvec_d[ik];
            klist->wk[ik + nk] = klist->wk[ik];
        }
    }

    std::ofstream ofs_kpts_coarse(out_dir + "kpts_coarse.dat");
    ofs_kpts_coarse << "kpts_coarse:" << nk << std::setw(16) << "( Cartesian"
        << std::setw(36) << "|                Direct )" << std::setw(15)
        << "| wk (normalized as sum = nk)" << std::endl;
    for (int ik = 0; ik < nks; ++ik)
    {
        ofs_kpts_coarse << std::setw(5) << ik << std::setw(12) << klist->kvec_c[ik].x
            << std::setw(12) << klist->kvec_c[ik].y << std::setw(12) << klist->kvec_c[ik].z
            << " | " << std::setw(12) << klist->kvec_d[ik].x << std::setw(12)
            << klist->kvec_d[ik].y << std::setw(12) << klist->kvec_d[ik].z << " | "
            << klist->wk[ik] * nk << std::endl;
    }
}

void RI_kRlist::read_kpts_fine(const std::string& file, const UnitCell& ucell,
                               K_Vectors* const klist, const bool is_weighted,
                               const std::string& out_dir)
{
    std::ifstream ifs(file);
    if (!ifs) throw std::runtime_error(file + " not found");
    int nk;
    if (is_weighted) { ifs >> nk; ifs.ignore(2048, '\n'); }
    else { ifs >> nk >> nk >> nk >> nk; }

    const int nks = (this->nspin == 2) ? 2 * nk : nk;
    klist->set_nks(nks);
    klist->set_nkstot(nks);
    klist->set_nkstot_full(nk);
    klist->kvec_c.resize(nks);
    klist->kvec_d.resize(nks);
    klist->wk.resize(nks);
    klist->isk.resize(0);
    klist->ngk.resize(0);

    for (int ik = 0; ik < nk; ++ik)
    {
        ifs >> klist->kvec_d[ik].x >> klist->kvec_d[ik].y >> klist->kvec_d[ik].z;
        if (is_weighted)
        {
            ifs >> klist->wk[ik];
            klist->wk[ik] /= double(nk);
        }
        else { klist->wk[ik] = 1.0 / double(nk); }
        klist->kvec_c[ik] = klist->kvec_d[ik] * ucell.G;
        set_zero_if_close(klist->kvec_c[ik]);
    }
    std::cout << "Read " << nk << " k-points and weights from " << file << std::endl;
    if (this->nspin == 2)
    {
        for (int ik = 0; ik < nk; ++ik)
        {
            klist->kvec_c[ik + nk] = klist->kvec_c[ik];
            klist->kvec_d[ik + nk] = klist->kvec_d[ik];
            klist->wk[ik + nk] = klist->wk[ik];
        }
    }
    std::ofstream ofs_kpts_fine(out_dir + "kpts_fine.dat");
    ofs_kpts_fine << "kpts_fine:" << nk << std::setw(18) << "( Cartesian"
        << std::setw(36) << "|                Direct )" << std::setw(15)
        << "| wk (normalized as sum = nk)" << std::endl;
    for (int ik = 0; ik < nk; ++ik)
    {
        ofs_kpts_fine << std::setw(5) << ik << std::setw(12) << klist->kvec_c[ik].x
            << std::setw(12) << klist->kvec_c[ik].y << std::setw(12) << klist->kvec_c[ik].z
            << " | " << std::setw(12) << klist->kvec_d[ik].x << std::setw(12)
            << klist->kvec_d[ik].y << std::setw(12) << klist->kvec_d[ik].z << " | "
            << klist->wk[ik] * nk << std::endl;
    }
}

} // namespace LR_IO
