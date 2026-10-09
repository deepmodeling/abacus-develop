#include "write_pdos_text.h"

#include "source_base/global_variable.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

void ModuleIO::write_pdos_text(
        const UnitCell& ucell,
        const ModuleBase::matrix* pdos,
        const int nspin,
        const int nlocal,
        const int npoints,
        const double& emin,
        const double& dos_edelta_ev,
        const std::string& out_dir,
        const std::string& basis,
        const int istep)
{
    ModuleBase::TITLE("ModuleIO", "write_pdos_text");

    // number of spin channels to write
    const int nspin0 = (nspin == 2) ? 2 : 1;

    for (int is = 0; is < nspin0; ++is)
    {
        std::stringstream ss;
        ss << out_dir << "pdoss" << is + 1;
        if (istep >= 0)
        {
            ss << "g" << istep + 1;
        }
        ss << "_" << basis << ".txt";

        std::ofstream ofs(ss.str().c_str());

        // iw2m stores m as 0..2l, mapping to physical m values:
        // 0->0, 1->+1, 2->-1, 3->+2, 4->-2, 5->+3, 6->-3
        // So columns are ordered: s, p(m=0,+1,-1), d(m=0,+1,-1,+2,-2), f(...)
        ofs << "# istep: " << istep + 1 << std::endl;
        ofs << "# npoints: " << npoints << std::endl;
        ofs << "# energy(eV)  atom  species  pdos(1/eV), columns: s(m=0) p(m=0,+1,-1) d(m=0,+1,-1,+2,-2) f(m=0,+1,-1,+2,-2,+3,-3)" << std::endl;

        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            const int ia = ucell.iat2ia[iat];
            const int it = ucell.iat2it[iat];
            const Atom* atom = &ucell.atoms[it];
            const int s0 = ucell.itiaiw2iwt(it, ia, 0);
            const int max_l = atom->nwl;

            // collect all (l, m) for this atom, ordered s,p,d,f
            // iw2m stores m as 0..2l (non-negative), so loop m = 0..2*L
            struct OrbInfo
            {
                int l;
                int m; // 0..2l
            };
            std::vector<OrbInfo> orb_list;

            for (int L = 0; L <= max_l; ++L)
            {
                for (int m = 0; m <= 2 * L; ++m)
                {
                    for (int j = 0; j < atom->nw; ++j)
                    {
                        if (atom->iw2l[j] == L && atom->iw2m[j] == m)
                        {
                            orb_list.push_back({L, m});
                            break;
                        }
                    }
                }
            }

            for (int n = 0; n < npoints; ++n)
            {
                const double en = emin + n * dos_edelta_ev;

                ofs << std::setw(12) << std::fixed << std::setprecision(6) << en
                    << "  " << std::setw(3) << iat + 1
                    << "  " << std::setw(4) << ucell.atoms[it].label;

                for (size_t io = 0; io < orb_list.size(); ++io)
                {
                    const int L = orb_list[io].l;
                    const int m = orb_list[io].m;
                    double pdos_val = 0.0;

                    // sum over all zeta for this (l, m)
                    for (int j = 0; j < atom->nw; ++j)
                    {
                        if (atom->iw2l[j] != L || atom->iw2m[j] != m)
                        {
                            continue;
                        }
                        const int w = ucell.itiaiw2iwt(it, ia, j);

                        if (nspin == 4)
                        {
                            const int w0 = w - s0;
                            pdos_val += pdos[0](s0 + 2 * w0, n) + pdos[0](s0 + 2 * w0 + 1, n);
                        }
                        else
                        {
                            pdos_val += pdos[is](w, n);
                        }
                    }

                    // zero out negligible values
                    if (std::abs(pdos_val) < 1e-6)
                    {
                        pdos_val = 0.0;
                    }

                    ofs << "  " << std::setw(8) << std::fixed << std::setprecision(6) << pdos_val;
                }
                ofs << std::endl;
            }
        }

        ofs.close();
    }
}
