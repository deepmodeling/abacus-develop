#include "spectrum_mo.h"
#include "source_cell/klist.h"
#include "source_base/constants.h"
#include <cassert>
#include <fstream>
#include <iomanip>
namespace LR_Util
{
/// @brief output the velocity matrix in KS presentation in human-read friendly format
void output_spectrum_mo(const std::vector<std::complex<double>>& out_spectrum_mo,
                        const std::string& filename,
                        const double* const eig_ks,
                        const int nk,
                        const int nspin_tmp,
                        const int KS_num,
                        const K_Vectors& kv)
{
    assert(out_spectrum_mo.size() == nspin_tmp * 3 * nk * KS_num * KS_num);
    std::ofstream ofs(filename);
    ofs << "Data Unit (Hartree * Bohr)" << std::endl;
    ofs << "NOTICE: KS_index are restricted in nocc and nvirt" << std::endl;

    int step = nk * KS_num * KS_num;
    int offset = 0;
    int ipair = 0;
    for (int is = 0; is < nspin_tmp; ++is)
    {
        ofs << "ispin: " << is << std::endl;
        for (int ik = 0;ik < nk;++ik)
        {
            ofs << "k-point: " << ik << " " << kv.kvec_d[ik] << std::endl;
            ofs << std::setw(4) << "KS1" << std::setw(10) << "E1(eV)" << std::setw(6) << "KS2" << std::setw(10) << "E2(eV)"
                << std::setw(16) << "x" << std::setw(23) << "|x|^2" << std::setw(18) << "y" << std::setw(23) <<"|y|^2"
                << std::setw(18) << "z" << std::setw(23) <<"|z|^2" << std::setw(13) << "average" << std::endl;

            offset = (is * 3 * nk + ik) * KS_num * KS_num;
            for (int i = 0; i < KS_num; ++i)
            {
                for (int j = i; j < KS_num ; ++j)
                {
                    ipair = offset + i * KS_num + j;
                    double average = (std::norm(out_spectrum_mo[ipair]) + std::norm(out_spectrum_mo[ipair + step]) + std::norm(out_spectrum_mo[ipair + 2 * step])) / 3.0;
                    ofs << std::setw(3) << i << std::setw(12) << std::setprecision(6) << eig_ks[i + ik*KS_num] * ModuleBase::Ry_to_eV
                    << std::setw(4) << j << std::setw(12) << eig_ks[j + ik*KS_num] * ModuleBase::Ry_to_eV
                    << std::setw(28) << out_spectrum_mo[ipair] << std::setw(13) << std::norm(out_spectrum_mo[ipair])
                    << std::setw(28) << out_spectrum_mo[ipair + step] << std::setw(13) << std::norm(out_spectrum_mo[ipair + step])
                    << std::setw(28) << out_spectrum_mo[ipair + 2 * step]  << std::setw(13) << std::norm(out_spectrum_mo[ipair + 2 * step] )
                    << std::setw(13) << average << std::endl;
                }
            }
        }
    }
    ofs.close();
}

/// @brief output the velocity matrix in KS presentation in LibRPA format
void output_spectrum_mo_librpa(const std::vector<std::complex<double>>& out_spectrum_mo,
                        const std::string& filename,
                        const int nk, const int nspin_tmp, const int nbands, const int nbasis,
                        const int KS_num/*= nocc+nvirt, can be different from nbands*/)
{
    assert(out_spectrum_mo.size() == nspin_tmp * 3 * nk * KS_num * KS_num);
    std::ofstream ofs(filename);
    ofs << std::scientific << nk << std::endl;
    ofs << nspin_tmp << std::endl;
    ofs << nbands << std::endl;
    ofs << nbasis << std::endl;
    double HaBohrToEvAng = ModuleBase::Hartree_to_eV * ModuleBase::BOHR_TO_A; // 27.211396 * 0.5291770
    int offset = 0;
    int ipair = 0;
    for (int is = 0; is < nspin_tmp; ++is)
    {
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int id = 0; id < 3; ++id)
            {
                offset = (is * 3 * nk + id * nk + ik) * KS_num * KS_num;
                ofs <<"  " << id + 1 << "  " << ik + 1 << "  " << is + 1 << std::endl;
                for (int ib1 = 0; ib1 < KS_num; ++ib1)
                {
                    for (int ib2 = 0; ib2 < KS_num; ++ib2)
                    {
                        ipair = offset + ib1 * KS_num + ib2;
                        ofs << std::setw(26) << std::fixed << std::setprecision(16) << out_spectrum_mo[ipair].real() * HaBohrToEvAng
                            << std::setw(26) << out_spectrum_mo[ipair].imag() * HaBohrToEvAng << std::endl;
                    }
                }
            }

        }
    }
}}
