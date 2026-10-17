#include "source_pw/module_pwdft/dftu_base_io.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_cell/unitcell.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

namespace
{

// Consume the label of one SOC block in the compact layout, e.g.
// "spin 1 nelec 1.07922327" or "spin 12 re". The optional
// "nelec <value>" token pair is skipped when present.
void skip_compact_soc_label(std::ifstream& ifs, char* word)
{
    ifs >> word; // "spin"
    ifs >> word; // block index: "1", "2" or "12"
    ifs >> word; // "nelec", "re" or "im"
    if (strcmp(word, "nelec") == 0)
    {
        double nelec = 0.0;
        ifs >> nelec;
    }
}

// Read n*n doubles of one compact matrix block row by row. When out is
// non-null the values are stored, otherwise they are read and discarded.
void read_matrix_block(std::ifstream& ifs, const int n, double* out)
{
    for (int i = 0; i < n * n; ++i)
    {
        double value = 0.0;
        ifs >> value;
        if (ifs.fail())
        {
            ModuleBase::WARNING_QUIT("read_matrix_block",
                                     "failed to read a matrix element from the occupation matrix file");
        }
        if (out != nullptr)
        {
            out[i] = value;
        }
    }
}

} // namespace

namespace DFTU_BASE
{

void read_occup_m(const UnitCell& ucell,
                  OccupationMatrix& occ,
                  const std::vector<int>& l_channel,
                  const int init_occ_mat,
                  const std::string& fn,
                  const std::string& init_chg,
                  int nspin,
                  int npol,
                  OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "read_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    std::ifstream ifdftu(fn.c_str(), std::ios::in);

    if (!ifdftu)
    {
        if (init_occ_mat > 0)
        {
            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "Can not find the file dm_onsite_ini.txt. Please check your dm_onsite_ini.txt");
        }
        else
        {
            if (init_chg == "file")
            {
                const std::string not_found_msg = "Can not find the file " + fn
                                                  + ". Please do scf calculation first";
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", not_found_msg);
            }
        }
        const std::string open_failed_msg = "Can not open " + fn + " file";
        ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", open_failed_msg);
    }

    GlobalV::ofs_running << " DFT+U: read occupation matrix from " << fn << std::endl;
    std::cout << " DFT+U: read occupation matrix from " << fn << std::endl;

    ifdftu.clear();
    ifdftu.seekg(0);

    char word[20];

    int T = 0;
    int iat = 0;
    int spin = 0;
    int L = 0;
    int zeta = 0;

    while (true)
    {
        if (!(ifdftu >> word))
        {
            break;
        }

        // Comment lines ("# ...") carry the provenance header and the
        // electronic-step marker in the compact layout; skip the line.
        if (word[0] == '#')
        {
            ifdftu.ignore(150, '\n');
            continue;
        }

        // The legacy "Atom=" token and the compact "Atom" label both
        // start one atom block; every other token (element prefix, log
        // sentences) is ignored.
        const bool legacy_atom = (strcmp(word, "Atom=") == 0);
        const bool compact_atom = (strcmp(word, "Atom") == 0);
        if (!legacy_atom && !compact_atom)
        {
            continue;
        }
        const bool legacy_layout = legacy_atom;

        ifdftu >> iat;
        iat -= 1;

        if (legacy_layout)
        {
            ifdftu >> word;

            // Accept both the compact "L=2" and the legacy split "L= 2".
            if (strncmp("L=", word, 2) != 0)
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            if (strlen(word) > 2)
            {
                L = atoi(word + 2);
            }
            else
            {
                ifdftu >> L;
            }
            ifdftu >> word;

            // The ORBITAL= token is optional in newly written files; when
            // absent the radial channel defaults to 0.
            if (strncmp("ORBITAL=", word, 8) == 0)
            {
                if (strlen(word) > 8)
                {
                    zeta = atoi(word + 8);
                }
                else
                {
                    ifdftu >> zeta;
                }
            }
            else
            {
                zeta = 0;
            }
            ifdftu.ignore(150, '\n');

            if (zeta != 0)
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m",
                                         "only the first radial channel (ORBITAL=0) is supported");
            }
        }
        else
        {
            // Compact header, e.g. " Fe Atom 1 L 2 mag 0.03..."
            ifdftu >> word;
            if (strcmp(word, "L") != 0)
            {
                ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            ifdftu >> L;
            // The optional "mag ..." tail is redundant for a restart.
            ifdftu.ignore(150, '\n');
            zeta = 0;
        }

        T = ucell.iat2it[iat];
        const int NL = ucell.atoms[T].nwl + 1;

        for (int l = 0; l < NL; l++)
        {
            if (l != l_channel[T])
            {
                continue;
            }

            const int nm = 2 * L + 1;

            if (nspin == 1 || nspin == 2)
            {
                for (int is = 0; is < 2; is++)
                {
                    if (legacy_layout)
                    {
                        ifdftu >> word;
                        if (strcmp("spin=", word) != 0)
                        {
                            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
                        }
                        ifdftu >> spin;
                        spin -= 1;
                        ifdftu.ignore(150, '\n');

                        double value = 0.0;
                        for (int m0 = 0; m0 < nm; m0++)
                        {
                            for (int m1 = 0; m1 < nm; m1++)
                            {
                                ifdftu >> value;
                                occ.set(iat, L, spin, m0, m1, value);
                            }
                            ifdftu.ignore(150, '\n');
                        }
                    }
                    else
                    {
                        // Compact label: "spin 1 nelec 0.468..."
                        ifdftu >> word;
                        if (strcmp("spin", word) != 0)
                        {
                            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
                        }
                        int spin_idx = 0;
                        ifdftu >> spin_idx;
                        if (spin_idx != is + 1)
                        {
                            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG SPIN INDEX IN Plus_U FILE");
                        }
                        ifdftu >> word;
                        if (strcmp("nelec", word) != 0)
                        {
                            ModuleBase::WARNING_QUIT("DFTU_BASE::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
                        }
                        double nelec = 0.0;
                        ifdftu >> nelec;

                        double value = 0.0;
                        for (int m0 = 0; m0 < nm; m0++)
                        {
                            for (int m1 = 0; m1 < nm; m1++)
                            {
                                ifdftu >> value;
                                occ.set(iat, L, is, m0, m1, value);
                            }
                        }
                    }
                }
            }
            else if (nspin == 4) // SOC
            {
                if (legacy_layout)
                {
                    double value = 0.0;
                    for (int m0 = 0; m0 < nm; m0++)
                    {
                        for (int ipol0 = 0; ipol0 < npol; ipol0++)
                        {
                            const int m0_all = m0 + nm * ipol0;

                            for (int m1 = 0; m1 < nm; m1++)
                            {
                                for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                {
                                    int m1_all = m1 + nm * ipol1;
                                    ifdftu >> value;
                                    occ.set(iat, L, 0, m0_all, m1_all, value);
                                }
                            }
                            ifdftu.ignore(150, '\n');
                        }
                    }
                }
                else
                {
                    // Four labeled blocks, in the same order they are
                    // written: up/up, up/down Re, up/down Im, down/down.
                    std::vector<double> uu(nm * nm);
                    std::vector<double> re(nm * nm);
                    std::vector<double> im(nm * nm);
                    std::vector<double> dd(nm * nm);

                    skip_compact_soc_label(ifdftu, word);
                    read_matrix_block(ifdftu, nm, uu.data());

                    skip_compact_soc_label(ifdftu, word);
                    read_matrix_block(ifdftu, nm, re.data());

                    // The imaginary block is needed by the PW (Pauli)
                    // path to reconstruct b2 = 2*Im(n_ud); the LCAO
                    // (real spin-basis) path discards it because its
                    // storage can not represent Im(n_ud).
                    skip_compact_soc_label(ifdftu, word);
                    read_matrix_block(ifdftu, nm, im.data());

                    skip_compact_soc_label(ifdftu, word);
                    read_matrix_block(ifdftu, nm, dd.data());

                    if (soc_layout == SOC_LAYOUT_PAULI)
                    {
                        // Reconstruct the 4 contiguous Pauli blocks
                        // [b0, b1, b2, b3] in the 2m x 2m flat buffer
                        // from the file's (n_uu, Re(n_ud), Im(n_ud),
                        // n_dd) representation. The writer prints
                        //   n_uu     = (b0 + b3)/2
                        //   Re(n_ud) = b1/2
                        //   Im(n_ud) = b2/2
                        //   n_dd     = (b0 - b3)/2
                        // so the inverse is
                        //   b0 = n_uu + n_dd
                        //   b1 = 2 * Re(n_ud)
                        //   b2 = 2 * Im(n_ud)
                        //   b3 = n_uu - n_dd.
                        ModuleBase::matrix& occ0 = occ.mat(iat, L, 0);
                        const int m2 = nm * nm;
                        for (int k = 0; k < m2; ++k)
                        {
                            occ0.c[0 * m2 + k] = uu[k] + dd[k];
                            occ0.c[1 * m2 + k] = 2.0 * re[k];
                            occ0.c[2 * m2 + k] = 2.0 * im[k];
                            occ0.c[3 * m2 + k] = uu[k] - dd[k];
                        }
                    }
                    else // SOC_LAYOUT_SPIN_BASIS_REAL
                    {
                        // Reconstruct the real 2m x 2m spin-basis matrix.
                        // Mirror the legacy writer: the Re block fills both
                        // (up, down) and (down, up) at the same (m0, m1).
                        // Im(n_ud) is dropped because the storage can not
                        // represent it.
                        for (int m0 = 0; m0 < nm; m0++)
                        {
                            for (int m1 = 0; m1 < nm; m1++)
                            {
                                const int k = m0 * nm + m1;
                                occ.set(iat, L, 0, m0, m1, uu[k]);
                                occ.set(iat, L, 0, m0, nm + m1, re[k]);
                                occ.set(iat, L, 0, nm + m0, m1, re[k]);
                                occ.set(iat, L, 0, nm + m0, nm + m1, dd[k]);
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

} // namespace DFTU_BASE
