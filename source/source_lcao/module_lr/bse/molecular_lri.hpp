//=======================
// AUTHOR : Ziqing Guan
// DATE :   2026-03-22
//=======================
#include "molecular_lri.h"
#include <RI/distribute/Distribute_Equally.h>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif
namespace BSE
{

template <typename T>
void MolecularLRI<T>::init(TLRI<T>& Cs_in, TLRI<T>& Vs_in, TLRI<T>& Ws_in, const Exx_Info_RI& info_ri)
{
    ModuleBase::TITLE("MolecularLRI", "init");
    ModuleBase::timer::start("MolecularLRI", "init");

    // 1. build kindex_map: index → fractional coordinate for all k-points
    std::vector<Tk> kindex_map_in(this->nk);
    for (int ik = 0; ik < this->nk; ++ik)
        kindex_map_in[ik] = RI_Util::Vector3_to_array3(this->kv.kvec_d.at(ik));
    this->LR_lri.init(kindex_map_in, this->nocc, this->nvirt);
    
    // 2. distribute atom and k-point tasks among MPI processes
    int nproc = GlobalV::NPROC;
    int task_sizes = this->ucell.nat * this->ucell.nat * this->nk * this->nk;
    std::cout << "Total Molecular LRI tasks: " << task_sizes << ", number of MPI processes: " << nproc << std::endl;
    this->LR_lri.set_parallel(MPI_COMM_WORLD, (std::size_t)this->ucell.nat, (std::size_t)this->nk, this->kRlist.period);
    for (int k1 : this->LR_lri.k1_indices)
    {
        this->is_local_k1[k1] = true;
    }

    int proc_ntasks = this->LR_lri.list_I.size() * this->LR_lri.list_J.size() * this->LR_lri.k1_indices.size() * this->LR_lri.k2_indices.size();
    GlobalV::ofs_running << "Molecular LRI init: Process " << GlobalV::MY_RANK
        << " handles " << proc_ntasks << " tasks." << std::endl;
    print_a(GlobalV::ofs_running, this->LR_lri.list_I, "list_I");
    print_a(GlobalV::ofs_running, this->LR_lri.list_J, "list_J");
    print_k(GlobalV::ofs_running, this->LR_lri.kindex_map, this->LR_lri.k1_indices, "list_k1");
    print_k(GlobalV::ofs_running, this->LR_lri.kindex_map, this->LR_lri.k2_indices, "list_k2");
    // LRI_CV_Tools::write_Cs_ao(Cs_in, PARAM.globalv.global_out_dir + "Cs_in_test_" + std::to_string(GlobalV::MY_RANK));
    // LRI_CV_Tools::write_Vs_abf(Vs_in, PARAM.globalv.global_out_dir + "Vs_in_test_" + std::to_string(GlobalV::MY_RANK));
    // LRI_CV_Tools::write_Vs_abf(Ws_in, PARAM.globalv.global_out_dir + "Ws_in_test_" + std::to_string(GlobalV::MY_RANK));

    // build q→kpair mapping (for W term q-approximation)
    this->build_q_to_kpair_map(PARAM.inp.bse_q_approx_mode, PARAM.inp.bse_q_approx_threshold);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "distribute_atom_and_k");

    // 3. move R of tensors to nearest image Bvk cell for valid Fourier interpolation
    double dist;
    std::set<int> all_atoms;
    for (int i = 0; i < this->ucell.nat; ++i)
    {
        all_atoms.insert(i);
        for (int j = 0; j < this->ucell.nat; ++j)
        {
            for (const TC R_original : this->kRlist.Rlist)
            {
				const TC R = cell_nearest.cell_nearest_direction(i, j, R_original, dist);
                if (R != R_original)
                {
                    BSE_Util::move_R_tensor(Cs_in, i, j, R_original, R);
                    BSE_Util::move_R_tensor(Vs_in, i, j, R_original, R);
                    BSE_Util::move_R_tensor(Ws_in, i, j, R_original, R);
                }
            }
        }
    }
    std::set<TA> set_I(this->LR_lri.list_I.begin(), this->LR_lri.list_I.end());
    std::set<TA> set_J(this->LR_lri.list_J.begin(), this->LR_lri.list_J.end());
    std::set<TA> set_IJ(this->LR_lri.list_IJ.begin(), this->LR_lri.list_IJ.end());

    // 4. set tensors, in these functions MPI distribution will be performed
    ModuleBase::TITLE("MolecularLRI", "before_set_Cs");
    this->LR_lri.set_Cs(Cs_in, info_ri.C_threshold, set_IJ, all_atoms);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "set_Cs");
    ModuleBase::TITLE("MolecularLRI", "before_set_Vs");
    this->LR_lri.set_Vs(Vs_in, info_ri.V_threshold, set_I, set_J);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "set_Vs");
    ModuleBase::TITLE("MolecularLRI", "before_set_Ws");
    this->LR_lri.set_Ws(Ws_in, info_ri.V_threshold, set_I, set_J);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "set_Ws");
    if (PARAM.inp.out_ri_cv)
    {        
        TLRI<T>& Cs_LRI = this->LR_lri.lrik.data_pool.at("Cs_").Ds_ab;// see LRI::set_tensor_map2
        LRI_CV_Tools::write_Cs_ao(Cs_LRI, PARAM.globalv.global_out_dir + "Cs_lrik_test_" + std::to_string(GlobalV::MY_RANK));
        TLRI<T>& Vs_LRI = this->LR_lri.lrik.data_pool.at("Vs_").Ds_ab;
        LRI_CV_Tools::write_Vs_abf(Vs_LRI, PARAM.globalv.global_out_dir + "Vs_lrik_test_" + std::to_string(GlobalV::MY_RANK));
        TLRI<T>& Ws_LRI = this->LR_lri.lrik.data_pool.at("Ws_").Ds_ab;
        LRI_CV_Tools::write_Vs_abf(Ws_LRI, PARAM.globalv.global_out_dir + "Ws_lrik_test_" + std::to_string(GlobalV::MY_RANK));
    }
    // 5. prepare mo-type tensors  
    this->LR_lri.map_psi = this->transform_psi_k(this->psi_ks, this->LR_lri.k_indices);

    ModuleBase::TITLE("MolecularLRI", "cal_Csk_ao_mo");
    ModuleBase::timer::start("MolecularLRI", "cal_Csk_ao_mo");
    this->LR_lri.cal_Csk_ao_mo("Cs_", GlobalV::ofs_running);
    ModuleBase::timer::end("MolecularLRI", "cal_Csk_ao_mo");
    
    ModuleBase::TITLE("MolecularLRI", "before_free_Cs");
    this->LR_lri.free_Cs(); // free Cs_ao to save memory
    ModuleBase::TITLE("MolecularLRI", "after_free_Cs");
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "MolecularLRI_init");
    ModuleBase::timer::end("MolecularLRI", "init");
}

/// @brief build q_list and q2kpair mapping according to approximation mode
template <typename T>
void MolecularLRI<T>::build_q_to_kpair_map(int mode, double threshold)
{
    ModuleBase::TITLE("MolecularLRI", "build_q_to_kpair_map");
    using namespace RI::Array_Operator;
    const Tk k_unit{1.0, 1.0, 1.0};

    // Fuzzy comparator for q_set
    struct kComparator {
        bool operator()(const Tk& lhs, const Tk& rhs) const {
            constexpr double eps = 1e-6;
            for (int i = 0; i < 3; ++i)
                if (std::abs(lhs[i] - rhs[i]) > eps) return lhs[i] < rhs[i];
            return false;
        }
    };

    auto frac_to_cart = [&](const Tk& k_frac) -> Tk {
        ModuleBase::Vector3<double> v_frac(k_frac[0], k_frac[1], k_frac[2]);
        ModuleBase::Vector3<double> v_cart = v_frac * ucell.G;
        return Tk{v_cart.x, v_cart.y, v_cart.z};
    };
    auto cart_to_frac = [&](const Tk& k_cart) -> Tk {
        ModuleBase::Vector3<double> v_cart(k_cart[0], k_cart[1], k_cart[2]);
        ModuleBase::Vector3<double> v_frac = v_cart * ucell.latvec.Transpose();
        return Tk{v_frac.x, v_frac.y, v_frac.z};
    };
    auto cart_dist = [&](const Tk& q1, const Tk& q2) -> double {
        Tk dq_frac = cart_to_frac({q1[0]-q2[0], q1[1]-q2[1], q1[2]-q2[2]});
        for (int i = 0; i < 3; ++i) dq_frac[i] -= std::round(dq_frac[i]);
        Tk dq_cart = frac_to_cart(dq_frac);
        return std::sqrt(dq_cart[0]*dq_cart[0] + dq_cart[1]*dq_cart[1] + dq_cart[2]*dq_cart[2]);
    };

    if (mode == 0)
    {
        std::map<Tk, std::vector<std::pair<int, int>>, kComparator> q2kpair_tmp;
        for (int k1 : this->LR_lri.k1_indices)
            for (int k2 : this->LR_lri.k2_indices)
            {
                Tk q = (this->LR_lri.kindex_map[k2] - this->LR_lri.kindex_map[k1]) % k_unit;
                q2kpair_tmp[q].emplace_back(k1, k2);
            }
        this->LR_lri.q2kpair.insert(std::make_move_iterator(q2kpair_tmp.begin()), std::make_move_iterator(q2kpair_tmp.end()));
    }
    else
    {
        std::set<Tk, kComparator> q_coarse_set;
        const K_Vectors& kv_coarse = this->kRlist.klist_coarse;
        int nk_coarse = kv_coarse.get_nkstot_full();
        for (int ik1 = 0; ik1 < nk_coarse; ++ik1)
        {
            Tk ck1 = RI_Util::Vector3_to_array3(kv_coarse.kvec_d.at(ik1));
            for (int ik2 = 0; ik2 < nk_coarse; ++ik2)
            {
                Tk ck2 = RI_Util::Vector3_to_array3(kv_coarse.kvec_d.at(ik2));
                q_coarse_set.insert((ck2 - ck1) % k_unit);
            }
        }
        std::vector<Tk> q_coarse_list(q_coarse_set.begin(), q_coarse_set.end());
        std::vector<Tk> q_coarse_cart;
        for (const Tk& q : q_coarse_list) q_coarse_cart.push_back(frac_to_cart(q));

        for (int k1 : this->LR_lri.k1_indices)
        {
            Tk k1_frac = this->LR_lri.kindex_map[k1];
            for (int k2 : this->LR_lri.k2_indices)
            {
                Tk k2_frac = this->LR_lri.kindex_map[k2];
                Tk q_frac = (k2_frac - k1_frac) % k_unit;
                Tk q_cart = frac_to_cart(q_frac);
                double q_norm = std::sqrt(q_cart[0]*q_cart[0] + q_cart[1]*q_cart[1] + q_cart[2]*q_cart[2]);

                Tk q_target;
                if (mode == 2 && q_norm < threshold)
                {
                    q_target = q_frac;
                }
                else
                {
                    double min_dist = std::numeric_limits<double>::max();
                    int nearest = -1;
                    for (std::size_t iq = 0; iq < q_coarse_cart.size(); ++iq)
                    {
                        double d = cart_dist(q_cart, q_coarse_cart[iq]);
                        if (d < min_dist) { min_dist = d; nearest = iq; }
                    }
                    q_target = q_coarse_list[nearest];
                }
                this->LR_lri.q2kpair[q_target].emplace_back(k1, k2);
            }
        }
    }

    for (auto& pair : this->LR_lri.q2kpair)
    {
        this->LR_lri.q_list.push_back(pair.first);
    }
    // Print q_list with both direct and Cartesian coordinates
    GlobalV::ofs_running << "q_list: size = " << this->LR_lri.q_list.size() << std::endl;
    GlobalV::ofs_running << "detailed q_list see 'qlist_{rank}.dat'" << std::endl;
    std::ofstream ofs_qlist(PARAM.globalv.global_out_dir + "qlist_" + std::to_string(GlobalV::MY_RANK) + ".dat");
    ofs_qlist << "q_list: size = " << this->LR_lri.q_list.size() << std::endl;
    ofs_qlist << "(direct coords)       | (Cartesian, Bohr^-1)" << std::endl;
    ofs_qlist << std::fixed << std::setprecision(4);
    for (const Tk& q : this->LR_lri.q_list)
    {
        Tk qc = frac_to_cart(q);
        ofs_qlist << "(" << std::setw(6) << q[0] << "," << std::setw(6) << q[1]
            << "," << std::setw(6) << q[2] << ") | (" << std::setw(7) << qc[0]
            << "," << std::setw(7) << qc[1] << "," << std::setw(7) << qc[2] << ")" << std::endl;
    }
    ofs_qlist << std::defaultfloat;
}

/// @brief transform psi to <k, <iat, tensor{nmo, iat.nw}>>, mo is not sliced
template <typename T>
std::map<int, std::map<TA, RI::Tensor<T>>>
MolecularLRI<T>::transform_psi_k(const psi::Psi<T>& psi_ks, const std::vector<int>& k_list)
{
    ModuleBase::TITLE("MolecularLRI", "transform_psi_k");
    ModuleBase::timer::start("MolecularLRI", "transform_psi_k");
    std::map<int, std::map<TA, RI::Tensor<T>>> psi_map;
    const std::size_t nmo = psi_ks.get_nbands();
    for (const int ik : k_list) // initialize
    {
        auto& psi_map_k = psi_map[ik];
        for (int iat = 0; iat < this->ucell.nat; ++iat)
        {
            psi_map_k[iat];
        }
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (const auto& ik : k_list)
    {
        auto& psi_map_k = psi_map.at(ik);
        for (int iat = 0; iat < this->ucell.nat; ++iat)
        {
            const int it = this->ucell.iat2it[iat];
            const std::size_t nw = this->ucell.atoms[it].nw;
            RI::Tensor<T> t({nmo, nw});
            for (int im = 0; im < nmo; ++im)
            {
                for (int iw = 0; iw < nw; ++iw)
                {
                    t(im, iw) = psi_ks(ik, im, this->ucell.get_iat2iwt()[iat]+iw);
                }
            }
            psi_map_k.at(iat) = std::move(t);
        }
    }
    ModuleBase::timer::end("MolecularLRI", "transform_psi_k");
    return psi_map;
}

}// namespace BSE