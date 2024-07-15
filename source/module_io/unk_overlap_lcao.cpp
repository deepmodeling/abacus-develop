#include "unk_overlap_lcao.h"

#include "ctime"
#include "module_base/scalapack_connector.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

unkOverlap_lcao::unkOverlap_lcao()
{
    allocate_flag = false;
}

unkOverlap_lcao::~unkOverlap_lcao()
{
    if (allocate_flag)
    {
        for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
        {
            delete [] cal_tag[iw];
        }
        delete [] cal_tag;
    }

    // GlobalV::ofs_running << "this is ~unkOverlap_lcao()" << std::endl;
}

void unkOverlap_lcao::init(const Grid_Technique& gt, const int nkstot)
{
    // std::cout << "unkOverlap_lcao::init start" << std::endl;

    int Lmax_used, Lmax;
    int exx_lmax = 0;
#ifdef __EXX
    exx_lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif

    const int ntype = GlobalC::ORB.get_ntype();
    int lmax_orb = -1, lmax_beta = -1;
    for (int it = 0; it < ntype; it++)
    {
        lmax_orb = std::max(lmax_orb, GlobalC::ORB.Phi[it].getLmax());
        lmax_beta = std::max(lmax_beta, GlobalC::ucell.infoNL.Beta[it].getLmax());
    }
    const double dr = GlobalC::ORB.get_dR();
    const double dk = GlobalC::ORB.get_dk();
    const int kmesh = GlobalC::ORB.get_kmesh() * 4 + 1;
    int Rmesh = static_cast<int>(GlobalC::ORB.get_Rmax() / dr) + 4;
    Rmesh += 1 - Rmesh % 2;

    Center2_Orb::init_Table_Spherical_Bessel(2,
                                             3,
                                             Lmax_used,
                                             Lmax,
                                             exx_lmax,
                                             lmax_orb,
                                             lmax_beta,
                                             dr,
                                             dk,
                                             kmesh,
                                             Rmesh,
                                             psb_);

    ModuleBase::Ylm::set_coefficients();

    MGT.init_Gaunt_CH(Lmax);
    MGT.init_Gaunt(Lmax);

    const int T = 0;                                                    // any selected element type
    orb_r.set_orbital_info(GlobalC::ORB.Phi[T].PhiLN(0, 0).getLabel(),  // atom label
                           T,                                           // atom type
                           1,                                           // angular momentum L
                           1,                                           // number of orbitals of this L , just N
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getNr(),     // number of radial mesh
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getRab(),    // the mesh interval in radial mesh
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getRadial(), // radial mesh value(a.u.)
                           Numerical_Orbital_Lm::Psi_Type::Psi,
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getRadial(), // radial wave function
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getNk(),
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getDk(),
                           GlobalC::ORB.Phi[T].PhiLN(0, 0).getDruniform(),
                           false,
                           true,
                           GlobalV::CAL_FORCE);

    // array initialization
    allocate_flag = true;
    this->kpoints_number = nkstot;
    if (allocate_flag)
    {
        cal_tag = new int*[GlobalV::NLOCAL];
        for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
        {
            cal_tag[iw] = new int[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(cal_tag[iw], GlobalV::NLOCAL);
        }
    }

#ifdef __MPI
    // parallel scheme
    int nproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    const int total_term = GlobalV::NLOCAL * GlobalV::NLOCAL;
    const int remain = total_term % nproc;
    int local_term = total_term / nproc;
    if (myrank < remain)
    {
        local_term++;
    }
    int start;
    for (int rank = 0; rank < nproc; rank++)
    {
        if (rank == myrank)
        {
            if (myrank < remain)
            {
                start = myrank * local_term;
            }
            else
            {
                start = myrank * local_term + remain;
            }
        }
    }
#else
    int start = 0;
    int local_term = GlobalV::NLOCAL * GlobalV::NLOCAL;
#endif
    int count = -1;
    for (int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
    {
        for (int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
        {
            count++;
            if (count >= start && count < (start + local_term))
            {
                cal_tag[iw1][iw2] = 1;
            }
        }
    }

    for (int TA = 0; TA < GlobalC::ucell.ntype; TA++)
    {
        for (int TB = 0; TB < GlobalC::ucell.ntype; TB++)
        {
            for (int LA = 0; LA <= GlobalC::ORB.Phi[TA].getLmax(); LA++)
            {
                for (int NA = 0; NA < GlobalC::ORB.Phi[TA].getNchi(LA); ++NA)
                {
                    for (int LB = 0; LB <= GlobalC::ORB.Phi[TB].getLmax(); ++LB)
                    {
                        for (int NB = 0; NB < GlobalC::ORB.Phi[TB].getNchi(LB); ++NB)
                        {
                            center2_orb11[TA][TB][LA][NA][LB].insert(
                                std::make_pair(NB,
                                               Center2_Orb::Orb11(GlobalC::ORB.Phi[TA].PhiLN(LA, NA),
                                                                  GlobalC::ORB.Phi[TB].PhiLN(LB, NB),
                                                                  psb_,
                                                                  MGT)));
                        }
                    }
                }
            }
        }
    }

    for (int TA = 0; TA < GlobalC::ucell.ntype; TA++)
    {
        for (int TB = 0; TB < GlobalC::ucell.ntype; TB++)
        {
            for (int LA = 0; LA <= GlobalC::ORB.Phi[TA].getLmax(); LA++)
            {
                for (int NA = 0; NA < GlobalC::ORB.Phi[TA].getNchi(LA); ++NA)
                {
                    for (int LB = 0; LB <= GlobalC::ORB.Phi[TB].getLmax(); ++LB)
                    {
                        for (int NB = 0; NB < GlobalC::ORB.Phi[TB].getNchi(LB); ++NB)
                        {
                            center2_orb21_r[TA][TB][LA][NA][LB].insert(
                                std::make_pair(NB,
                                               Center2_Orb::Orb21(GlobalC::ORB.Phi[TA].PhiLN(LA, NA),
                                                                  orb_r,
                                                                  GlobalC::ORB.Phi[TB].PhiLN(LB, NB),
                                                                  psb_,
                                                                  MGT)));
                        }
                    }
                }
            }
        }
    }

    for (auto& co1: center2_orb11)
        for (auto& co2: co1.second)
            for (auto& co3: co2.second)
                for (auto& co4: co3.second)
                    for (auto& co5: co4.second)
                        for (auto& co6: co5.second)
                            co6.second.init_radial_table();

    for (auto& co1: center2_orb21_r)
        for (auto& co2: co1.second)
            for (auto& co3: co2.second)
                for (auto& co4: co3.second)
                    for (auto& co5: co4.second)
                        for (auto& co6: co5.second)
                            co6.second.init_radial_table();

    // std::cout << "unkOverlap_lcao::init end" << std::endl;
    return;
}

int unkOverlap_lcao::iw2it(int iw)
{
    int ic, type;
    ic = 0;
    type = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            for (int L = 0; L < GlobalC::ucell.atoms[it].nwl + 1; L++)
            {
                for (int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for (int i = 0; i < (2 * L + 1); i++)
                    {
                        if (ic == iw)
                        {
                            type = it;
                        }
                        ic++;
                    }
                }
            }
        }
    }
    return type;
}

int unkOverlap_lcao::iw2ia(int iw)
{
    int ic, na;
    ic = 0;
    na = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            for (int L = 0; L < GlobalC::ucell.atoms[it].nwl + 1; L++)
            {
                for (int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for (int i = 0; i < (2 * L + 1); i++)
                    {
                        if (ic == iw)
                        {
                            na = ia;
                        }
                        ic++;
                    }
                }
            }
        }
    }
    return na;
}

int unkOverlap_lcao::iw2iL(int iw)
{
    int ic, iL;
    ic = 0;
    iL = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            for (int L = 0; L < GlobalC::ucell.atoms[it].nwl + 1; L++)
            {
                for (int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for (int i = 0; i < (2 * L + 1); i++)
                    {
                        if (ic == iw)
                        {
                            iL = L;
                        }
                        ic++;
                    }
                }
            }
        }
    }
    return iL;
}

int unkOverlap_lcao::iw2iN(int iw)
{
    int ic, iN;
    ic = 0;
    iN = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            for (int L = 0; L < GlobalC::ucell.atoms[it].nwl + 1; L++)
            {
                for (int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for (int i = 0; i < (2 * L + 1); i++)
                    {
                        if (ic == iw)
                        {
                            iN = N;
                        }
                        ic++;
                    }
                }
            }
        }
    }
    return iN;
}

int unkOverlap_lcao::iw2im(int iw)
{
    int ic, im;
    ic = 0;
    im = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            for (int L = 0; L < GlobalC::ucell.atoms[it].nwl + 1; L++)
            {
                for (int N = 0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for (int i = 0; i < (2 * L + 1); i++)
                    {
                        if (ic == iw)
                        {
                            im = i;
                        }
                        ic++;
                    }
                }
            }
        }
    }
    return im;
}

// search for the nearest neighbor atoms
void unkOverlap_lcao::cal_R_number()
{
    // The number of overlaps between atomic orbitals 1 and atomic orbitals 2,
    // or the number of R, is empty when there is no overlap
    orb1_orb2_R.resize(GlobalV::NLOCAL);
    for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
    {
        orb1_orb2_R[iw].resize(GlobalV::NLOCAL);
    }

    ModuleBase::Vector3<double> tau1, tau2, dtau;
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);

            for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];
                const double R_direct_x = (double)GlobalC::GridD.getBox(ad).x;
                const double R_direct_y = (double)GlobalC::GridD.getBox(ad).y;
                const double R_direct_z = (double)GlobalC::GridD.getBox(ad).z;

                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
                if (distance < rcut - 1.0e-15)
                {
                    // translate: the unit of R_car is GlobalC::ucell.lat0
                    ModuleBase::Vector3<double> R_car = R_direct_x * GlobalC::ucell.a1 + R_direct_y * GlobalC::ucell.a2
                                                        + R_direct_z * GlobalC::ucell.a3;

                    for (int iw1 = 0; iw1 < atom1->nw; iw1++)
                    {
                        int orb_index_in_NLOCAL_1 = GlobalC::ucell.itiaiw2iwt(T1, I1, iw1);
                        for (int iw2 = 0; iw2 < atom2->nw; iw2++)
                        {
                            int orb_index_in_NLOCAL_2 = GlobalC::ucell.itiaiw2iwt(T2, I2, iw2);
                            orb1_orb2_R[orb_index_in_NLOCAL_1][orb_index_in_NLOCAL_2].push_back(R_car);
                        } // end iw2

                    } // end iw1
                }

            } //  end ad

        } // end I1

    } // end T1

    return;
}

void unkOverlap_lcao::cal_orb_overlap()
{
    // std::cout << "the cal_orb_overlap is start" << std::endl;
    psi_psi.resize(GlobalV::NLOCAL);
    psi_r_psi.resize(GlobalV::NLOCAL);
    for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
    {
        psi_psi[iw].resize(GlobalV::NLOCAL);
        psi_r_psi[iw].resize(GlobalV::NLOCAL);
    }

    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);

    for (int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
    {
        for (int iw2 = 0; iw2 < GlobalV::NLOCAL; iw2++)
        {
            // if ( !pv.in_this_processor(iw1,iw2) ) continue;

            // iw1 and iw2 never have overlap
            if (orb1_orb2_R[iw1][iw2].empty())
                continue;

            int atomType1 = iw2it(iw1);
            int ia1 = iw2ia(iw1);
            int N1 = iw2iN(iw1);
            int L1 = iw2iL(iw1);
            int m1 = iw2im(iw1);
            int atomType2 = iw2it(iw2);
            int ia2 = iw2ia(iw2);
            int N2 = iw2iN(iw2);
            int L2 = iw2iL(iw2);
            int m2 = iw2im(iw2);

            for (int iR = 0; iR < orb1_orb2_R[iw1][iw2].size(); iR++)
            {
                ModuleBase::Vector3<double> r_distance
                    = (GlobalC::ucell.atoms[atomType2].tau[ia2] - GlobalC::ucell.atoms[atomType1].tau[ia1]
                       + orb1_orb2_R[iw1][iw2][iR])
                      * GlobalC::ucell.lat0;
                psi_psi[iw1][iw2].push_back(
                    center2_orb11[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap(origin_point,
                                                                                       r_distance,
                                                                                       m1,
                                                                                       m2));

                double overlap_x = -1 * sqrt(ModuleBase::FOUR_PI / 3.0)
                                   * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap(origin_point,
                                                                                                          r_distance,
                                                                                                          m1,
                                                                                                          1,
                                                                                                          m2); // m = 1
                double overlap_y = -1 * sqrt(ModuleBase::FOUR_PI / 3.0)
                                   * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap(origin_point,
                                                                                                          r_distance,
                                                                                                          m1,
                                                                                                          2,
                                                                                                          m2); // m = -1
                double overlap_z = sqrt(ModuleBase::FOUR_PI / 3.0)
                                   * center2_orb21_r[atomType1][atomType2][L1][N1][L2].at(N2).cal_overlap(origin_point,
                                                                                                          r_distance,
                                                                                                          m1,
                                                                                                          0,
                                                                                                          m2); // m =0
                ModuleBase::Vector3<double> overlap(overlap_x, overlap_y, overlap_z);

                psi_r_psi[iw1][iw2].push_back(overlap);
            }
        }
    }

    // std::cout << "the cal_orb_overlap is end" << std::endl;
    return;
}

void unkOverlap_lcao::prepare_midmatrix_pblas(const int ik_L,
                                              const int ik_R,
                                              const ModuleBase::Vector3<double> dk,
                                              std::complex<double>*& midmatrix,
                                              const Parallel_Orbitals& pv,
                                              const K_Vectors& kv)
{
    // ModuleBase::Vector3<double> dk = kv.kvec_c[ik_R] - kv.kvec_c[ik_L];
    midmatrix = new std::complex<double>[pv.nloc];
    ModuleBase::GlobalFunc::ZEROS(midmatrix, pv.nloc);
    for (int iw_row = 0; iw_row < GlobalV::NLOCAL; iw_row++) // global
    {
        for (int iw_col = 0; iw_col < GlobalV::NLOCAL; iw_col++) // global
        {
            int ir = pv.global2local_row(iw_row); // local
            int ic = pv.global2local_col(iw_col); // local

            if (ir >= 0 && ic >= 0)
            {
                int index = ic * pv.nrow + ir;
                ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[iw2it(iw_row)].tau[iw2ia(iw_row)];
                for (int iR = 0; iR < orb1_orb2_R[iw_row][iw_col].size(); iR++)
                {
                    double kRn = (kv.kvec_c[ik_R] * orb1_orb2_R[iw_row][iw_col][iR] - dk * tau1) * ModuleBase::TWO_PI;
                    std::complex<double> kRn_phase(cos(kRn), sin(kRn));
                    std::complex<double> orb_overlap(psi_psi[iw_row][iw_col][iR],
                                                     (-dk * GlobalC::ucell.tpiba * psi_r_psi[iw_row][iw_col][iR]));
                    midmatrix[index] = midmatrix[index] + kRn_phase * orb_overlap;
                }
            }
        }
    }
}

std::complex<double> unkOverlap_lcao::det_berryphase(const int ik_L,
                                                     const int ik_R,
                                                     const ModuleBase::Vector3<double> dk,
                                                     const int occ_bands,
                                                     const Parallel_Orbitals& para_orb,
                                                     const psi::Psi<std::complex<double>>* psi_in,
                                                     const K_Vectors& kv)
{
    const std::complex<double> minus = std::complex<double>(-1.0, 0.0);
    std::complex<double> det = std::complex<double>(1.0, 0.0);
    std::complex<double>* midmatrix = nullptr;
    std::complex<double>* C_matrix = new std::complex<double>[para_orb.nloc];
    std::complex<double>* out_matrix = new std::complex<double>[para_orb.nloc];
    ModuleBase::GlobalFunc::ZEROS(C_matrix, para_orb.nloc);
    ModuleBase::GlobalFunc::ZEROS(out_matrix, para_orb.nloc);

    this->prepare_midmatrix_pblas(ik_L, ik_R, dk, midmatrix, para_orb, kv);

    char transa = 'C';
    char transb = 'N';
    int occBands = occ_bands;
    int nlocal = GlobalV::NLOCAL;
    std::complex<double> alpha = {1.0, 0.0}, beta = {0.0, 0.0};
    int one = 1;
#ifdef __MPI
    pzgemm_(&transa,
            &transb,
            &occBands,
            &nlocal,
            &nlocal,
            &alpha,
            &psi_in[0](ik_L, 0, 0),
            &one,
            &one,
            para_orb.desc,
            midmatrix,
            &one,
            &one,
            para_orb.desc,
            &beta,
            C_matrix,
            &one,
            &one,
            para_orb.desc);

    pzgemm_(&transb,
            &transb,
            &occBands,
            &occBands,
            &nlocal,
            &alpha,
            C_matrix,
            &one,
            &one,
            para_orb.desc,
            &psi_in[0](ik_R, 0, 0),
            &one,
            &one,
            para_orb.desc,
            &beta,
            out_matrix,
            &one,
            &one,
            para_orb.desc);

    int* ipiv = new int[para_orb.nrow];
    int info;
    pzgetrf_(&occBands, &occBands, out_matrix, &one, &one, para_orb.desc, ipiv, &info);

    for (int i = 0; i < occBands; i++) // global
    {
        int ir = para_orb.global2local_row(i); // local
        int ic = para_orb.global2local_col(i); // local
        if (ir >= 0 && ic >= 0)
        {
            int index = ic * para_orb.nrow + ir;
            if (ipiv[ir] != (i + 1))
            {
                det = minus * det * out_matrix[index];
            }
            else
            {
                det = det * out_matrix[index];
            }
        }
    }
    delete[] ipiv;
#endif
    delete[] midmatrix;
    delete[] C_matrix;
    delete[] out_matrix;

#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the GlobalV::KPAR = 1.
    std::complex<double> result;
    MPI_Allreduce(&det, &result, 1, MPI_DOUBLE_COMPLEX, MPI_PROD, DIAG_WORLD);
    return result;
#endif

    return det;
}