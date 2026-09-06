#include "to_w90_lcao.h"

#include "source_lcao/module_ri/abfs_vector3_order.h"
#include "source_io/module_parameter/parameter.h"
#include "fr_overlap.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_hamilt/module_hcontainer/atom_pair.h"

#include <fstream>
#include <functional>

#ifdef __LCAO
toW90_LCAO::toW90_LCAO(const bool& out_wannier_mmn,
                                   const bool& out_wannier_amn,
                                   const bool& out_wannier_unk,
                                   const bool& out_wannier_eig,
                                   const bool& out_wannier_wvfn_formatted,
                                   const std::string& nnkpfile,
                                   const std::string& wannier_spin,
                                   const int& nspin,
                                   const int& nbands,
                                   const int& nqx,
                                   const double& dq,
                                   const int& npol,
                                   const LCAO_Orbitals& orb)
    : toW90(out_wannier_mmn,
                  out_wannier_amn,
                  out_wannier_unk,
                  out_wannier_eig,
                  out_wannier_wvfn_formatted,
                  nnkpfile,
                  wannier_spin,
                  nspin,
                  nbands,
                  nqx,
                  dq,
                  npol),
    orb_(orb)
{
}

toW90_LCAO::~toW90_LCAO()
{
}

void toW90_LCAO::calculate(const UnitCell& ucell,
                                 const Grid_Driver& gd,
                                 const ModuleBase::matrix& ekb,
                                 const K_Vectors& kv,
                                 const psi::Psi<std::complex<double>>& psi,
                                 const Parallel_Orbitals* pv)
{
    this->ParaV = pv;

    read_nnkp(ucell,kv);

    if (nspin_ == 2)
    {
        if (wannier_spin == "up")
        {
            start_k_index = 0;
        }
        else if (wannier_spin == "down")
        {
            start_k_index = num_kpts / 2;
        }
        else
        {
            ModuleBase::WARNING_QUIT("toW90::calculate", "Error wannier_spin set,is not \"up\" or \"down\" ");
        }
    }

    if (out_wannier_mmn || out_wannier_amn)
    {
        iw2it.resize(PARAM.globalv.nlocal);
        iw2ia.resize(PARAM.globalv.nlocal);
        iw2iL.resize(PARAM.globalv.nlocal);
        iw2iN.resize(PARAM.globalv.nlocal);
        iw2im.resize(PARAM.globalv.nlocal);
        iw2iorb.resize(PARAM.globalv.nlocal);

        std::map<size_t, std::map<size_t, std::map<size_t, size_t>>> temp_orb_index;
        int count = 0;
        for (int it = 0; it < ucell.ntype; it++)
        {
            for (int iL = 0; iL < ucell.atoms[it].nwl + 1; iL++)
            {
                for (int iN = 0; iN < ucell.atoms[it].l_nchi[iL]; iN++)
                {
                    temp_orb_index[it][iL][iN] = count;
                    count++;
                }
            }
        }

        int iw = 0;
        for (int it = 0; it < ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                for (int iL = 0; iL < ucell.atoms[it].nwl + 1; iL++)
                {
                    for (int iN = 0; iN < ucell.atoms[it].l_nchi[iL]; iN++)
                    {
                        for (int im = 0; im < (2 * iL + 1); im++)
                        {
                            iw2it[iw] = it;
                            iw2ia[iw] = ia;
                            iw2iL[iw] = iL;
                            iw2iN[iw] = iN;
                            iw2im[iw] = im;
                            iw2iorb[iw] = temp_orb_index[it][iL][iN];

                            iw++;
                        }
                    }
                }
            }
        }

        initialize_orb_table(ucell);
        produce_basis_orb();
        set_R_coor(ucell, gd);
        count_delta_k(ucell,kv);
    }

    if (out_wannier_eig)
    {
        out_eig(ekb);
    }

    if (out_wannier_mmn)
    {
        int dk_size = delta_k_all.size();
        this->FR.resize(dk_size);
        std::vector<std::function<std::complex<double>(ModuleBase::Vector3<double>)>> fr_ptr(dk_size);
        for (int i = 0; i < dk_size; i++)
        {
            ModuleBase::Vector3<double> delta_k(delta_k_all[i].x, delta_k_all[i].y, delta_k_all[i].z);

            fr_ptr[i] = [delta_k](ModuleBase::Vector3<double> r) -> std::complex<double> {
                double phase = delta_k * r;
                std::complex<double> exp_idkr = std::exp(-1.0 * ModuleBase::IMAG_UNIT * phase);
                return exp_idkr;
            };

            FR[i].set_parameters(fr_ptr[i], &ucell, &orb_, &gd, ParaV, 140, 110);
            FR[i].calculate_FR();
        }

        cal_Mmn(ucell,kv, psi);
    }

    if (out_wannier_amn)
    {
        cal_Amn(ucell,kv, psi);
    }

    if (out_wannier_unk)
    {
        out_unk(psi);
    }
}

#endif
