#include "td_current.h"
#ifdef __LCAO
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_lcao/module_tddft/snap_psibeta_half_tddft.h"
#ifdef _OPENMP
#include <unordered_set>
#include <omp.h>
#endif

TD_current::TD_current(const UnitCell* ucell_in,
                      Grid_Driver* GridD_in,
                      const Parallel_Orbitals* paraV,
                      const LCAO_Orbitals& orb,
                      const TwoCenterIntegrator* intor)
    : ucell(ucell_in), paraV(paraV) , orb_(orb), Grid(GridD_in), intor_(intor)
{   
    // for length gague, the A(t) = 0 for all the time.
    this->cart_At = ModuleBase::Vector3<double>(0,0,0);
    this->initialize_vcomm_r(GridD_in, paraV);
    this->initialize_grad_term(GridD_in, paraV);
}
TD_current::~TD_current()
{
    for (int dir=0;dir<3;dir++)
    {
        delete this->current_term[dir];
    }
}
//allocate space for current_term
void TD_current::initialize_vcomm_r(Grid_Driver* GridD, const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("TD_current", "initialize_vcomm_r");
    ModuleBase::timer::tick("TD_current", "initialize_vcomm_r");
    for (int dir=0;dir<3;dir++)
    {
        if (this->current_term[dir] == nullptr)
        this->current_term[dir] = new hamilt::HContainer<std::complex<double>>(paraV);
    }

    this->adjs_vcommr.clear();
    this->adjs_vcommr.reserve(this->ucell->nat);
    for (int iat0 = 0; iat0 < ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell, tau0, T0, I0, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad1];
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0
                < orb_.Phi[T1].getRcut() + this->ucell->infoNL.Beta[T0].get_rcut_max())
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        this->adjs_vcommr.push_back(adjs);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
            {
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int iat2 = ucell->itia2iat(T2, I2);
                ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                if (paraV->get_col_size(iat2) <= 0 || paraV->get_row_size(iat1) <= 0)
                {
                    continue;
                }
                hamilt::AtomPair<std::complex<double>> tmp(iat1,
                                         iat2,
                                         R_index2.x - R_index1.x,
                                         R_index2.y - R_index1.y,
                                         R_index2.z - R_index1.z,
                                         paraV);
                for (int dir=0;dir<3;dir++)
                {
                    this->current_term[dir]->insert_pair(tmp);
                }
            }
        }
    }
    // allocate the memory of BaseMatrix in cal_vcomm_r_IJR, and set the new values to zero
    for (int dir=0;dir<3;dir++)
    {
        this->current_term[dir]->allocate(nullptr, true);
    }
    ModuleBase::timer::tick("TD_current", "initialize_vcomm_r");
}
void TD_current::initialize_grad_term(Grid_Driver* GridD, const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("TD_current", "initialize_grad_term");
    ModuleBase::timer::tick("TD_current", "initialize_grad_term");

    for (int dir=0;dir<3;dir++)
    {
        if (this->current_term[dir] == nullptr)
        this->current_term[dir] = new hamilt::HContainer<std::complex<double>>(paraV);
    }
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell, tau1, T1, I1, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T2 = adjs.ntype[ad1];
            const int I2 = adjs.natom[ad1];
            const int iat2 = ucell->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad1];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat1, iat2, R_index2).norm() * this->ucell->lat0
                < orb_.Phi[T1].getRcut() + orb_.Phi[T2].getRcut())
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        this->adjs_grad.push_back(adjs);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell->itia2iat(T2, I2);
            ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            hamilt::AtomPair<std::complex<double>> tmp(iat1, iat2, R_index.x, R_index.y, R_index.z, paraV);
            for (int dir=0;dir<3;dir++)
            {
                this->current_term[dir]->insert_pair(tmp);
            }
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    for (int dir=0;dir<3;dir++)
    {
        this->current_term[dir]->allocate(nullptr, true);
    }

    ModuleBase::timer::tick("TD_current", "initialize_grad_term");
}

void TD_current::calculate_vcomm_r()
{
    ModuleBase::TITLE("TD_current", "calculate_vcomm_r");
    ModuleBase::timer::tick("TD_current", "calculate_vcomm_r");

    const Parallel_Orbitals* paraV = this->current_term[0]->get_atom_pair(0).get_paraV();
    const int npol = this->ucell->get_npol();

    // 1. calculate <psi|beta> for each pair of atoms
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        AdjacentAtomInfo& adjs = this->adjs_vcommr[iat0];

        std::vector<std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>> nlm_tot;
        nlm_tot.resize(adjs.adj_num + 1);
        for (int i = 0; i < adjs.adj_num + 1; i++)
        {
            nlm_tot[i].resize(4);
        }

        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic)
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T1 = adjs.ntype[ad];
                const int I1 = adjs.natom[ad];
                const int iat1 = ucell->itia2iat(T1, I1);
                const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
                const Atom* atom1 = &ucell->atoms[T1];
                auto all_indexes = paraV->get_indexes_row(iat1);
                auto col_indexes = paraV->get_indexes_col(iat1);
                // insert col_indexes into all_indexes to get universal set with no repeat elements
                all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
                std::sort(all_indexes.begin(), all_indexes.end());
                all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());
                for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
                {
                    const int iw1 = all_indexes[iw1l] / npol;
                    std::vector<std::vector<std::complex<double>>> nlm;
                    // nlm is a vector of vectors, but size of outer vector is only 1 when out_current is false
                    // and size of outer vector is 4 when out_current is true (3 for <psi|r_i * exp(-iAr)|beta>, 1 for
                    // <psi|exp(-iAr)|beta>) inner loop : all projectors (L0,M0)

                    // snap_psibeta_half_tddft() are used to calculate <psi|exp(-iAr)|beta>
                    // and <psi|rexp(-iAr)|beta> as well if current are needed
                    
                    module_tddft::snap_psibeta_half_tddft(orb_,
                                                            this->ucell->infoNL,
                                                            nlm,
                                                            tau1 * this->ucell->lat0,
                                                            T1,
                                                            atom1->iw2l[iw1],
                                                            atom1->iw2m[iw1],
                                                            atom1->iw2n[iw1],
                                                            tau0 * this->ucell->lat0,
                                                            T0,
                                                            this->cart_At,
                                                            true);
                    for (int dir = 0; dir < 4; dir++)
                    {
                        nlm_tot[ad][dir].insert({all_indexes[iw1l], nlm[dir]});
                    }
                }
            }

    #ifdef _OPENMP
            // record the iat number of the adjacent atoms
            std::set<int> ad_atom_set;
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T1 = adjs.ntype[ad];
                const int I1 = adjs.natom[ad];
                const int iat1 = ucell->itia2iat(T1, I1);
                ad_atom_set.insert(iat1);
            }

            // split the ad_atom_set into num_threads parts
            const int num_threads = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
            std::set<int> ad_atom_set_thread;
            int i = 0;
            for(const auto iat1 : ad_atom_set)
            {
                if (i % num_threads == thread_id)
                {
                    ad_atom_set_thread.insert(iat1);
                }
                i++;
            }
    #endif

    // 2. calculate <psi_I|beta>D<beta|psi_{J,R}> for each pair of <IJR> atoms
            for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
            {
                const int T1 = adjs.ntype[ad1];
                const int I1 = adjs.natom[ad1];
                const int iat1 = ucell->itia2iat(T1, I1);
    #ifdef _OPENMP
                if (ad_atom_set_thread.find(iat1) == ad_atom_set_thread.end())
                    {
                        continue;
                    }
    #endif
                ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
                for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
                {
                    const int T2 = adjs.ntype[ad2];
                    const int I2 = adjs.natom[ad2];
                    const int iat2 = ucell->itia2iat(T2, I2);
                    ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                    ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                                    R_index2[1] - R_index1[1],
                                                    R_index2[2] - R_index1[2]);
                    std::complex<double>* tmp_c[3] = {nullptr, nullptr, nullptr};
                    for (int i = 0; i < 3; i++)
                    {
                        tmp_c[i] = this->current_term[i]->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2])->get_pointer();
                    }
                    // if not found , skip this pair of atoms
                    if (tmp_c[0] != nullptr)
                    {
                        this->cal_vcomm_r_IJR(iat1,
                                            iat2,
                                            T0,
                                            paraV,
                                            nlm_tot[ad1],
                                            nlm_tot[ad2],
                                            tmp_c);
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("TD_current", "calculate_vcomm_r");
}

// cal_HR_IJR()
void TD_current::cal_vcomm_r_IJR(
    const int& iat1,
    const int& iat2,
    const int& T0,
    const Parallel_Orbitals* paraV,
    const std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>& nlm1_all,
    const std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>& nlm2_all,
    std::complex<double>** current_mat_p)
{
    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = col_indexes.size() * is + is2;
        }
    }
    // calculate the local matrix
    const std::complex<double>* tmp_d = nullptr;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        // const std::vector<std::complex<double>>* nlm1 = &(nlm1_all[0].find(row_indexes[iw1l])->second);
        std::vector<const std::vector<std::complex<double>>*> nlm1;
        for (int dir = 0; dir < 4; dir++)
        {
            nlm1.push_back(&(nlm1_all[dir].find(row_indexes[iw1l])->second));
        }

        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            std::vector<const std::vector<std::complex<double>>*> nlm2;
            for (int dir = 0; dir < 4; dir++)
            {
                nlm2.push_back(&(nlm2_all[dir].find(col_indexes[iw2l])->second));
            }
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            for (int is = 0; is < npol * npol; ++is)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    std::complex<double> nlm_r_tmp = std::complex<double>{0, 0};
                    std::complex<double> imag_unit = std::complex<double>{0, 1};
                    for (int no = 0; no < this->ucell->atoms[T0].ncpp.non_zero_count_soc[is]; no++)
                    {
                        const int p1 = this->ucell->atoms[T0].ncpp.index1_soc[is][no];
                        const int p2 = this->ucell->atoms[T0].ncpp.index2_soc[is][no];
                        this->ucell->atoms[T0].ncpp.get_d(is, p1, p2, tmp_d);
                        //<psi|rexp(-iAr)|beta><beta|exp(iAr)|psi>-<psi|exp(-iAr)|beta><beta|rexp(iAr)|psi>
                        // multiply d in the end
                        nlm_r_tmp += (nlm1[dir + 1]->at(p1) * std::conj(nlm2[0]->at(p2))
                                        - nlm1[0]->at(p1) * std::conj(nlm2[dir + 1]->at(p2)))
                                        * (*tmp_d);
                    }
                    // -i[r,Vnl], 2.0 due to the unit transformation
                    current_mat_p[dir][step_trace[is]] -= imag_unit * nlm_r_tmp / 2.0;
                }
            }
            for (int dir = 0; dir < 3; dir++)
            {
                current_mat_p[dir] += npol;
            }
        }
        for (int dir = 0; dir < 3; dir++)
        {
            current_mat_p[dir] += (npol - 1) * col_indexes.size();
        }
    }
}

void TD_current::calculate_grad_term()
{
    ModuleBase::TITLE("TD_current", "calculate_grad_term");
    if(this->current_term[0]==nullptr || this->current_term[0]->size_atom_pairs()<=0)
    {
        ModuleBase::WARNING_QUIT("TD_current::calculate_grad_term", "grad_term is nullptr or empty");
    }
    ModuleBase::timer::tick("TD_current", "calculate_grad_term");

    const Parallel_Orbitals* paraV = this->current_term[0]->get_atom_pair(0).get_paraV();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iat1 = 0; iat1 < this->ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo& adjs = this->adjs_grad[iat1];
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = ucell->itia2iat(T2, I2);
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad];
            ModuleBase::Vector3<double> dtau = this->ucell->cal_dtau(iat1, iat2, R_index2);

            std::complex<double>* tmp_c[3] = {nullptr, nullptr, nullptr};
            for (int i = 0; i < 3; i++)
            {
                tmp_c[i] = this->current_term[i]->find_matrix(iat1, iat2, R_index2)->get_pointer();
            }
            if (tmp_c[0] != nullptr)
            {
                this->cal_grad_IJR(iat1, iat2, paraV, dtau, tmp_c);
            }
            else
            {
                ModuleBase::WARNING_QUIT("TD_current::calculate_grad_term", "R_index not found in HR");
            }
        }
    }
    ModuleBase::timer::tick("TD_current", "calculate_grad_term");
}

void TD_current::cal_grad_IJR(const int& iat1,
                              const int& iat2,
                              const Parallel_Orbitals* paraV,
                              const ModuleBase::Vector3<double>& dtau,
                              std::complex<double>** current_mat_p)
{
    // ---------------------------------------------
    // get info of orbitals of atom1 and atom2 from ucell
    // ---------------------------------------------
    int T1, I1;
    this->ucell->iat2iait(iat1, &I1, &T1);
    int T2, I2;
    this->ucell->iat2iait(iat2, &I2, &T2);
    Atom& atom1 = this->ucell->atoms[T1];
    Atom& atom2 = this->ucell->atoms[T2];

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();

    const int* iw2l1 = atom1.iw2l.data();
    const int* iw2n1 = atom1.iw2n.data();
    const int* iw2m1 = atom1.iw2m.data();
    const int* iw2l2 = atom2.iw2l.data();
    const int* iw2n2 = atom2.iw2n.data();
    const int* iw2m2 = atom2.iw2m.data();
    // ---------------------------------------------
    // get tau1 (in cell <0,0,0>) and tau2 (in cell R)
    // in principle, only dtau is needed in this function
    // snap_psipsi should be refactored to use dtau directly
    // ---------------------------------------------
    const ModuleBase::Vector3<double>& tau1 = this->ucell->get_tau(iat1);
    const ModuleBase::Vector3<double> tau2 = tau1 + dtau;
    // ---------------------------------------------
    // calculate the Ekinetic matrix for each pair of orbitals
    // ---------------------------------------------
    double grad[3] = {0, 0, 0};
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    const int step_trace = col_indexes.size() + 1;
    for(int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const int iw1 = row_indexes[iw1l] / npol;
        const int L1 = iw2l1[iw1];
        const int N1 = iw2n1[iw1];
        const int m1 = iw2m1[iw1];
        // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
        int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const int iw2 = col_indexes[iw2l] / npol;
            const int L2 = iw2l2[iw2];
            const int N2 = iw2n2[iw2];
            const int m2 = iw2m2[iw2];
            // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
            int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;

            // calculate <psi|∇R|psi>, which equals to -<psi|∇r|psi>.
            intor_->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * this->ucell->lat0, nullptr, grad);
            ModuleBase::Vector3<double> grad_overlap(grad[0], grad[1], grad[2]);

            for (int dir = 0; dir < 3; dir++)
            {
                for (int ipol = 0; ipol < npol; ipol++)
                {
                    // part of Momentum operator, -i∇r,used to calculate the current
                    // here is actually i∇R
                    current_mat_p[dir][ipol * step_trace] += std::complex<double>(0, grad_overlap[dir]);
                }
                current_mat_p[dir] += npol;
            }
        }
        for (int dir = 0; dir < 3; dir++)
        {
            current_mat_p[dir] += (npol - 1) * col_indexes.size();
        }   
    }
}

#endif // __LCAO
