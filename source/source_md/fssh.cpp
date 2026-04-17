#include "fssh.h"

#include "source_esolver/esolver_ks.h"
#include "source_lcao/module_lr/esolver_lrtd_lcao.h"
#include "source_lcao/module_lr/utils/lr_util.hpp"
#include "source_base/parallel_2d.h"
#include "source_base/parallel_global.h"
#include "source_psi/psi.h"
#include "source_cell/print_cell.h"
#include <mpi.h>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <sstream>

// [FSSH修改说明] 增加宏定义以控制是否使用有限差分(FD)受力
// 1: 调用外部 abacus-fd 工具获取当前活性态的数值受力 (物理上更合理但极慢)
// 0: 使用 ESolver 返回的受力 (目前 ABACUS 仅支持基态受力)
#ifndef FSSH_USE_FD_FORCE
#define FSSH_USE_FD_FORCE 1
#endif

// ====================================================================
// Accessor: 用于访问ESolver的protected成员
// [FSSH修改说明] 大幅简化了原有的Accessor类:
// - 原代码: EsolverKsLcaoAccessor + EsolverLrAccessor 共含11个static方法，
//   其中包括get_X_ref, get_paraX_ref, get_paraC_ref, get_eig_ks_ref等
//   专门用于手动析构LR独占内存的方法(配合safe_destroy使用)
// - 修改后: 仅保留提取KS波函数和LR结果所需的最小接口
// - 删除原因: 使用借用构造函数后，LR可正常析构，无需手动定向销毁
// ====================================================================
namespace FsshIntegration {

    template <typename T, typename TR>
    class EsolverKsLcaoAccessor : public ModuleESolver::ESolver_KS_LCAO<T, TR> {
    public:
        static auto get_psi(ModuleESolver::ESolver_KS_LCAO<T, TR>* ks) -> decltype(static_cast<EsolverKsLcaoAccessor*>(ks)->psi)& { return static_cast<EsolverKsLcaoAccessor*>(ks)->psi; }
        static auto get_pelec(ModuleESolver::ESolver_KS_LCAO<T, TR>* ks) -> decltype(static_cast<EsolverKsLcaoAccessor*>(ks)->pelec)& { return static_cast<EsolverKsLcaoAccessor*>(ks)->pelec; }
    };

    template <typename T, typename TR>
    class EsolverLrAccessor : public LR::ESolver_LR<T, TR> {
    public:
        static const std::vector<ct::Tensor>* get_X(LR::ESolver_LR<T, TR>* esolver) { return &(static_cast<const EsolverLrAccessor*>(esolver)->X); }
        static int get_nstates(LR::ESolver_LR<T, TR>* esolver) { return static_cast<const EsolverLrAccessor*>(esolver)->nstates; }
        static auto get_pelec(LR::ESolver_LR<T, TR>* lr) -> decltype(static_cast<EsolverLrAccessor*>(lr)->pelec)& { return static_cast<EsolverLrAccessor*>(lr)->pelec; }
        // [FSSH修改说明] 新增: 获取 LR 模块实际使用的 nocc 和 nvirt
        // LR 模块的 nocc/nvirt 可能与 KS 的不同 (由 INPUT 参数 nocc/nvirt/nbands 控制)
        // X 系数索引为 X[i_lr * nvirt_lr + a_lr], 必须使用 LR 的维度而非 KS 的
        static int get_nocc(LR::ESolver_LR<T, TR>* esolver) { return static_cast<const EsolverLrAccessor*>(esolver)->nocc[0]; }
        static int get_nvirt(LR::ESolver_LR<T, TR>* esolver) { return static_cast<const EsolverLrAccessor*>(esolver)->nvirt[0]; }
        // [FSSH修改说明] 新增: 获取 X 系数的 Parallel_2D 分布信息
        // X 张量在 MPI 进程间通过 2D block-cyclic 分布 (paraX_)
        // 需要此信息将分布式 X 收集为全局 X 系数
        static const std::vector<Parallel_2D>& get_paraX(LR::ESolver_LR<T, TR>* esolver) { return static_cast<const EsolverLrAccessor*>(esolver)->paraX_; }
    };
} // namespace FsshIntegration

FsshMD::FsshMD(const Parameter& param_in, UnitCell& unit_in) 
    : MD_base(param_in, unit_in), fssh_initialized_(false) {
}

FsshMD::~FsshMD() {}

// Standard Verlet first-half integration
void FsshMD::first_half(std::ofstream& ofs) {
    update_vel(this->force);
    update_pos();
}

// Standard Verlet second-half integration
void FsshMD::second_half() {
    update_vel(this->force);
}

// Map flat MD_base::vel to UnitCell atoms' velocity
void FsshMD::sync_vel_to_ucell() {
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it) {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia) {
            ucell.atoms[it].vel[ia] = this->vel[iat];
            iat++;
        }
    }
}

// Map UnitCell atoms' velocity back to MD_base::vel after FSSH scaling
void FsshMD::sync_vel_from_ucell() {
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it) {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia) {
            this->vel[iat] = ucell.atoms[it].vel[ia];
            iat++;
        }
    }
}

void FsshMD::execute_hopping(ModuleESolver::ESolver* p_esolver, const Parameter& param_in) {
    if (!param_in.inp.cal_syns) return;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // ====================================================================
    // [FSSH修改说明] 重新设计 use_tddft 的判断逻辑
    // 原代码:
    //   bool use_tddft = (param_in.inp.esolver_type == "ks-lr" || param_in.inp.esolver_type == "lr");
    //   if (param_in.inp.esolver_type != "ksdft") { exit(1); }
    // 问题: use_tddft 在 esolver_type=="ks-lr" 时为 true，但紧接着因 != "ksdft" 而 exit(1)
    // 修改后:
    //   1. 仅当 esolver_type == "ks-lr" 时 use_tddft = true
    //   2. 同时支持 esolver_type == "ksdft" (纯KS模式) 和 "ks-lr" (KS+LR模式)
    // ====================================================================
    bool use_tddft = (param_in.inp.esolver_type == "ks-lr");

    if (param_in.inp.esolver_type != "ksdft" && param_in.inp.esolver_type != "ks-lr") {
        if (my_rank == 0) {
            GlobalV::ofs_warning << "[FSSH ERROR] Only esolver_type = ksdft or ks-lr is supported for FSSH!" << std::endl;
        }
        exit(1);
    }

    // ====================================================================
    // [FSSH修改说明] 读取新增 FSSH 输入参数
    // 原代码: fssh_nstate 和 fssh_init_state 通过硬编码或 out_dmr hack 获取
    // 修改后: 从 INPUT 参数 param_in.mdp.fssh_nstate / fssh_init_state 读取
    // ====================================================================
    const int fssh_nstate = param_in.mdp.fssh_nstate;
    const int fssh_init_state = param_in.mdp.fssh_init_state;

    int nbasis = 0;
    int nstates_ks = 0;
    int nocc = 0;  // [FSSH修改说明] 提升 nocc 作用域到函数顶层，供后续 init 调用使用
    ModuleBase::ComplexMatrix coef_new;
    std::vector<double> ks_bands;
    std::vector<double> tddft_energies;
    std::vector<CasidaWavefunction> casida_wfcs;
    bool extract_success = false;

    if (auto ks_gamma = dynamic_cast<ModuleESolver::ESolver_KS_LCAO<double, double>*>(p_esolver)) {

        auto psi_ptr = FsshIntegration::EsolverKsLcaoAccessor<double, double>::get_psi(ks_gamma);
        auto pelec_ptr = FsshIntegration::EsolverKsLcaoAccessor<double, double>::get_pelec(ks_gamma);

        // [FSSH修改说明] 使用全局维度代替局部维度
        // 原代码:
        //   nbasis = psi_ptr[0].get_nbasis();     // 返回 LOCAL nrow
        //   nstates_ks = psi_ptr[0].get_nbands(); // 返回 LOCAL ncol_bands
        // 问题: psi 通过 Parallel_2D 在 MPI 进程间 2D block-cyclic 分布，
        //   get_nbasis()/get_nbands() 返回当前进程的局部维度而非全局维度
        //   导致 coef_new 矩阵只有局部数据，TDA 伪波函数构造缺失大部分虚轨道
        // 修改后: 从 param_in 获取全局维度
        nbasis = param_in.globalv.nlocal;       // 全局基函数数 (nlocal)
        nstates_ks = param_in.inp.nbands;       // 全局 band 数 (nbands)

        // [FSSH修改说明] 收集分布式 psi 到全局 coef_new 矩阵
        // 原代码: 直接从 psi 局部数据构建 coef_new(nstates_ks_local, nbasis_local)
        // 问题: psi 在 MPI 进程间按 2D block-cyclic 分布:
        //   - 行维度: nlocal 个基函数分布到 nrow 个局部行
        //   - 列维度: nbands 个 band 分布到 ncol_bands 个局部列
        //   直接读局部数据导致 coef_new 只有 ~1/nprocs 的数据
        // 修改后: 使用 Parallel_Orbitals 的 local2global 映射 + MPI_Allreduce
        //   将每个进程的局部 psi 数据放到全局矩阵的正确位置，然后跨进程求和
        const auto& pv = ks_gamma->get_pv();
        const int nrow_local = psi_ptr[0].get_nbasis();   // 局部基函数数
        const int ncol_local = psi_ptr[0].get_nbands();   // 局部 band 数

        coef_new.create(nstates_ks, nbasis);
        {
            // 使用 double 临时数组进行 MPI 收集
            std::vector<double> full_psi(static_cast<size_t>(nstates_ks) * nbasis, 0.0);

            // 将局部 psi 数据放到全局矩阵的正确位置
            for (int j = 0; j < ncol_local; ++j) {          // 局部 band 索引
                int g_band = pv.local2global_col(j);         // 全局 band 索引
                if (g_band >= nstates_ks) continue;
                for (int i = 0; i < nrow_local; ++i) {      // 局部基函数索引
                    int g_basis = pv.local2global_row(i);    // 全局基函数索引
                    full_psi[g_band * nbasis + g_basis] = psi_ptr[0].operator()(0, j, i);
                }
            }

            // 跨所有进程求和 (每个进程贡献不同的元素，其余为 0)
#ifdef __MPI
            MPI_Allreduce(MPI_IN_PLACE, full_psi.data(),
                          nstates_ks * nbasis, MPI_DOUBLE, MPI_SUM, pv.comm());
#endif

            // 转换为 ComplexMatrix
            for (int ib = 0; ib < nstates_ks; ++ib) {
                for (int mu = 0; mu < nbasis; ++mu) {
                    coef_new(ib, mu) = std::complex<double>(full_psi[ib * nbasis + mu], 0.0);
                }
            }
        }

        // ekb/wg 存储全局数据 (所有 band 的能量和占据数)，可直接用全局索引访问
        ks_bands.resize(nstates_ks);
        for (int is = 0; is < nstates_ks; ++is) { ks_bands[is] = pelec_ptr->ekb(0, is); }

        // [FSSH修改说明] 从占据数(occupation weights)计算 nocc
        // ekb/wg 使用全局 band 索引，遍历所有 band
        nocc = 0;
        for (int ib = 0; ib < nstates_ks; ++ib) {
            if (pelec_ptr->wg(0, ib) > 0.5) nocc++;
        }

        // ====================================================================
        // [FSSH修改说明] 条件性调用 LR-TDDFT
        // 原代码: 无条件调用 LR::ESolver_LR 并始终 use_tddft_in_engine = true
        // 修改后:
        //   - 当 use_tddft == true (esolver_type == "ks-lr"):
        //     构造 LR 求解器，运行 Casida 方程获取激发能和 TDA 伪波函数
        //     两态能量差使用 TDDFT 计算结果
        //     非绝热耦合使用 TDA 伪波函数的 NAC 矩阵元
        //   - 当 use_tddft == false (esolver_type == "ksdft"):
        //     不调用 LR-TDDFT，使用 KSDFT 轨道能量差
        //     直接读取 KSDFT 的非绝热耦合矩阵元（通过 CSR 文件）
        // ====================================================================
        if (use_tddft) {
            LR::ESolver_LR<double, double> lr_solver(*ks_gamma, param_in.inp, this->ucell);

            // 运行 LR-TDDFT (Casida 方程)
            lr_solver.runner(this->ucell, param_in.mdp.md_nstep);

            // 提取 Casida 数据
            const std::vector<ct::Tensor>* X_tensor = FsshIntegration::EsolverLrAccessor<double, double>::get_X(&lr_solver);
            int num_excitations = FsshIntegration::EsolverLrAccessor<double, double>::get_nstates(&lr_solver);
            auto lr_pelec = FsshIntegration::EsolverLrAccessor<double, double>::get_pelec(&lr_solver);

            // [FSSH修改说明] 获取 LR 模块实际使用的 nocc/nvirt 维度
            // LR 模块的 nocc_lr/nvirt_lr 可能小于 KS 的 nocc/nvirt (由 INPUT 参数控制)
            // X 系数布局为 X[i_lr * nvirt_lr + a_lr], 必须使用 LR 维度正确索引
            int nocc_lr = FsshIntegration::EsolverLrAccessor<double, double>::get_nocc(&lr_solver);
            int nvirt_lr = FsshIntegration::EsolverLrAccessor<double, double>::get_nvirt(&lr_solver);

            // 态 0 为基态 (能量 = 0 作为参考)
            tddft_energies.push_back(0.0);
            casida_wfcs.push_back(CasidaWavefunction(0.0, {}, {}, nocc_lr, nvirt_lr));

            if (X_tensor && !X_tensor->empty() && lr_pelec) {
                const ct::Tensor& X_mat = (*X_tensor)[0];
                const double* x_data = X_mat.template data<double>();
                int nloc_per_band = X_mat.NumElements() / num_excitations;

                // [FSSH修改说明] 获取 paraX_ 以收集分布式 X 系数
                // 原代码: 直接复制局部 X 数据到 CasidaWavefunction::X_coeffs
                // 问题: X 张量通过 Parallel_2D (paraX_) 在 MPI 进程间 2D block-cyclic 分布，
                //   每个进程只有 nloc_per_band 个元素 (= nrow_local * ncol_local)，
                //   而全局大小为 nocc_lr * nvirt_lr。直接复制只获得局部数据，
                //   导致 ||X||² ≈ 1/nprocs (应为 ~1)，NAC 矩阵几乎为零
                // 修改后: 使用 LR_Util::gather_2d_to_full() 将分布式 X 收集到全局数组
                //   paraX_[0] 行维度 = nvirt_lr (虚轨道), 列维度 = nocc_lr (占据轨道)
                //   col_first=false → 输出布局: fullX[i_occ * nvirt_lr + a_virt]
                //   与 fssh_driver.cpp 中的索引 idx_ia = i_lr * nvirt_lr + a_lr 一致
                const auto& paraX = FsshIntegration::EsolverLrAccessor<double, double>::get_paraX(&lr_solver);
                const Parallel_2D& px = paraX[0];  // spin 0
                int full_size = nocc_lr * nvirt_lr;  // 全局 X 系数大小

                for (int i = 0; i < num_excitations; ++i) {
                    CasidaWavefunction wfc;
                    wfc.omega = lr_pelec->ekb(0, i);
                    wfc.nocc_lr = nocc_lr;
                    wfc.nvirt_lr = nvirt_lr;
                    tddft_energies.push_back(wfc.omega);

                    // 收集分布式 X 到全局数组
                    wfc.X_coeffs.resize(full_size, 0.0);
                    const double* band_data = x_data + i * nloc_per_band;
#ifdef __MPI
                    LR_Util::gather_2d_to_full(px, band_data, wfc.X_coeffs.data(),
                                                false, nvirt_lr, nocc_lr);
#else
                    // 串行: 直接复制 (nloc_per_band == full_size)
                    for (int j = 0; j < nloc_per_band; ++j) {
                        wfc.X_coeffs[j] = band_data[j];
                    }
#endif
                    casida_wfcs.push_back(wfc);
                }
            }

            if (my_rank == 0) {
                std::cout << "[FSSH INFO] LR-TDDFT completed. "
                          << tddft_energies.size() << " states available (incl. ground state)." << std::endl;

                // 调试: 检查 X 系数维度和量级
                int nvirt_check = nstates_ks - nocc;
                int expected_nloc = nocc * nvirt_check;
                std::ofstream dbg("fssh_nac_diagnostic.log", std::ios::app);
                dbg << std::scientific;
                dbg << "=== X coefficient debug ===\n";
                dbg << "  nocc_ks=" << nocc << " nvirt_ks=" << nvirt_check
                    << " nstates_ks=" << nstates_ks << " nbasis=" << nbasis << "\n";
                if (!casida_wfcs.empty() && casida_wfcs[0].nocc_lr > 0) {
                    dbg << "  nocc_lr=" << casida_wfcs[0].nocc_lr
                        << " nvirt_lr=" << casida_wfcs[0].nvirt_lr
                        << " expected_nloc(lr)=" << casida_wfcs[0].nocc_lr * casida_wfcs[0].nvirt_lr << "\n";
                }
                if (X_tensor && !X_tensor->empty()) {
                    const ct::Tensor& X_mat = (*X_tensor)[0];
                    int nlpb = X_mat.NumElements() / num_excitations;
                    int full_sz = nocc_lr * nvirt_lr;
                    dbg << "  X_mat.NumElements()=" << X_mat.NumElements()
                        << " num_excitations=" << num_excitations
                        << " nloc_per_band(local)=" << nlpb
                        << " full_size(gathered)=" << full_sz
                        << " expected(nocc*nvirt)=" << expected_nloc << "\n";
                    dbg << "  GATHERED: " << (casida_wfcs.size() > 1 && static_cast<int>(casida_wfcs[1].X_coeffs.size()) == full_sz ? "YES" : "NO") << "\n";
                }
                for (size_t I = 1; I < casida_wfcs.size() && I <= 2; ++I) {
                    const auto& X = casida_wfcs[I].X_coeffs;
                    dbg << "  State " << I << ": omega=" << casida_wfcs[I].omega
                        << " X_coeffs.size()=" << X.size() << "\n";
                    dbg << "    X values: ";
                    double max_x = 0.0, sum_x2 = 0.0;
                    for (size_t j = 0; j < X.size(); ++j) {
                        if (std::abs(X[j]) > max_x) max_x = std::abs(X[j]);
                        sum_x2 += X[j] * X[j];
                        if (j < 20) dbg << X[j] << " ";
                    }
                    dbg << "\n    max|X|=" << max_x << " ||X||²=" << sum_x2 << "\n";
                }
                dbg << "=== end X debug ===\n";
                dbg.close();
            }
        }
        // 当 use_tddft == false 时，tddft_energies 和 casida_wfcs 保持为空
        // fssh_driver 将根据 use_tddft 标志使用 KS 能量和 CSR 文件中的 NAC

        extract_success = true;
    }

    // ====================================================================
    // 交付 FSSH 动力学引擎推演
    // ====================================================================
    if (!fssh_initialized_) {
        double md_dt_au = param_in.mdp.md_dt * 41.341;

        // [FSSH修改说明] 使用新增输入参数初始化 FSSH 引擎
        // 原代码:
        //   int start_state = 1;
        //   if (param_in.inp.out_dmr.size() > 1) {
        //       start_state = param_in.inp.out_dmr[1];
        //   }
        //   fssh_engine_.init(nbasis, nstates_ks, md_dt_au, start_state);
        // 问题: 使用 out_dmr 参数 hack 来传递 start_state，语义不清且默认值硬编码为 1
        // 修改后: 使用 INPUT 文件中的 fssh_nstate 和 fssh_init_state 参数
        // [FSSH修改说明] init 增加 nocc 和 nstates_ks 参数
        // 原代码: fssh_engine_.init(nbasis, fssh_nstate, md_dt_au, fssh_init_state);
        // 修改后: 传递 nocc 和 nstates_ks 供 TDA NAC 使用
        fssh_engine_.init(nbasis, fssh_nstate, md_dt_au, fssh_init_state, nocc, nstates_ks);

        if (my_rank == 0) {
            std::cout << "[FSSH INFO] Initialized: nstate=" << fssh_nstate
                      << ", init_state=" << fssh_init_state
                      << ", nocc=" << nocc << ", nks_total=" << nstates_ks
                      << ", use_tddft=" << (use_tddft ? "True" : "False") << std::endl;
        }

        coef_old_ = coef_new;

        // [FSSH修改说明] 在初始化阶段缓存首步 Casida 数据
        // 原代码: (无此逻辑)
        // 修改原因: TDA NAC 公式需要前后两步的 X 系数 (X^I_{ia}(t) 和 X^I_{ia}(t+Δt))
        //   在 step 0 初始化时缓存 casida_wfcs，使 step 1 首次调用 run_step_advanced 时
        //   即可使用完整的 TDA NAC 计算 (casida_old=step0, casida_new=step1)
        if (use_tddft && !casida_wfcs.empty()) {
            fssh_engine_.set_casida_cache(casida_wfcs);
        }

        fssh_initialized_ = true;
    } else {
        std::string csr_file = param_in.globalv.global_out_dir + "syns_nao.csr";

        this->sync_vel_to_ucell();

        // [FSSH修改说明] 将 use_tddft 标志传给引擎
        // 原代码: bool use_tddft_in_engine = true;  // 硬编码为 true
        // 修改后: 直接传递根据 esolver_type 判断的 use_tddft
        fssh_engine_.run_step_advanced(
            coef_old_, coef_new, csr_file, ks_bands, this->ucell, use_tddft, tddft_energies, casida_wfcs
        );

        this->sync_vel_from_ucell();
        coef_old_ = coef_new;

#if FSSH_USE_FD_FORCE
        // [劫持逻辑] 调用有限差分工具获取当前活性态的真实受力并覆盖 Analytical Force
        this->update_force_active_state_fd(param_in);
#endif
    }
}

// ====================================================================
// [劫持逻辑实现] FsshMD::update_force_active_state_fd
// ====================================================================
void FsshMD::update_force_active_state_fd(const Parameter& param_in) {
    const int active_state = fssh_engine_.get_current_state();
    const int natom = ucell.nat;
    std::vector<double> fd_forces(natom * 3, 0.0);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (active_state == 0) {
        if (my_rank == 0) {
            std::cout << "[FSSH FD-HIJACK] Active state: 0 (Ground State). Using ESolver internal forces." << std::endl;
        }
        // Use the forces already present in this->force (MD_base::force)
        // No need for system call or file parsing.
        for (int i = 0; i < natom; ++i) {
            fd_forces[i * 3 + 0] = this->force[i].x;
            fd_forces[i * 3 + 1] = this->force[i].y;
            fd_forces[i * 3 + 2] = this->force[i].z;
        }
    } else {
        if (my_rank == 0) {
            std::cout << "[FSSH FD-HIJACK] Active state: " << active_state 
                      << ". Calling abacus-fd for numerical forces..." << std::endl;

            // 1. 落盘当前坐标到 STRU 文件
            unitcell::print_stru_file(this->ucell, this->ucell.atoms, this->ucell.latvec, "STRU", 
                                      param_in.inp.nspin, false, false, false, true, false, 0);

            // 2. 构造并执行系统调用
            // 简化 Python 脚本调用逻辑，使用绝对路径确保成功
            std::string py_script = std::string("import sys, os\n")
                + "from abacus_fd import main\n"
                + "abacus_bin = os.path.abspath('../../abacus-develop/build/abacus_2p')\n"
                + "sys.argv = ['abacus-fd', '" + (param_in.inp.esolver_type == "ks-lr" ? "kslr-all" : "gs-all") + "', '-a', abacus_bin, '-d', '.']\n"
                + "try:\n"
                + "    main()\n"
                + "except Exception as e:\n"
                + "    print(f'Python execution error: {e}')\n"
                + "    sys.exit(1)\n";

            std::ofstream ofs_py("run_fd.py");
            ofs_py << py_script;
            ofs_py.close();

            std::string cmd = "python3 run_fd.py > fd_force.log 2>&1";


            int ret = std::system(cmd.c_str());
            if (ret != 0) {
                std::cerr << "[FSSH ERROR] abacus-fd execution failed!" << std::endl;
            }

            // 3. 解析输出文件
            std::string force_file = (param_in.inp.esolver_type == "ks-lr") ? "excited_forces.txt" : "ground_forces.txt";
            std::ifstream ifs(force_file);
            if (ifs.is_open()) {
                std::string line;
                while (std::getline(ifs, line)) {
                    if (line.empty() || line[0] == '#') continue;
                    std::stringstream ss(line);
                    std::string type;
                    int s_idx, a_idx;
                    double fx, fy, fz;

                    if (param_in.inp.esolver_type == "ks-lr") {
                        ss >> type >> s_idx >> a_idx >> fx >> fy >> fz;
                        if (s_idx == active_state) {
                            if (a_idx < natom) {
                                fd_forces[a_idx * 3 + 0] = fx;
                                fd_forces[a_idx * 3 + 1] = fy;
                                fd_forces[a_idx * 3 + 2] = fz;
                            }
                        }
                    } else {
                        ss >> a_idx >> fx >> fy >> fz;
                        if (a_idx < natom) {
                            fd_forces[a_idx * 3 + 0] = fx;
                            fd_forces[a_idx * 3 + 1] = fy;
                            fd_forces[a_idx * 3 + 2] = fz;
                        }
                    }
                }
                ifs.close();
            } else {
                std::cerr << "[FSSH ERROR] Cannot open force file: " << force_file << std::endl;
            }
        }

        // 4. MPI 广播受力数据到所有进程
#ifdef __MPI
        MPI_Bcast(fd_forces.data(), natom * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    }

    // 5. 打印最终采用的受力信息 (到 fssh_md.log)
    if (my_rank == 0) {
        std::cout << " --- FSSH Nuclear Forces (eV/Angstrom) ---" << std::endl;
        std::cout << std::fixed << std::setprecision(8);
        for (int i = 0; i < natom; ++i) {
            std::cout << " Atom " << std::setw(3) << i << ": " 
                      << std::setw(14) << fd_forces[i*3+0] << " " 
                      << std::setw(14) << fd_forces[i*3+1] << " " 
                      << std::setw(14) << fd_forces[i*3+2] << std::endl;
        }
        std::cout << " ------------------------------------------" << std::endl;
    }

    // 6. 覆写内存中的力数组
    if (this->force != nullptr) {
        for (int i = 0; i < natom; ++i) {
            this->force[i].x = fd_forces[i * 3 + 0];
            this->force[i].y = fd_forces[i * 3 + 1];
            this->force[i].z = fd_forces[i * 3 + 2];
        }
    }

    // 同步到 UnitCell 对象
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it) {
        // 确保 Atom 对象的 force 向量大小正确
        if (ucell.atoms[it].force.size() != ucell.atoms[it].na) {
            ucell.atoms[it].force.resize(ucell.atoms[it].na);
        }
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia) {
            ucell.atoms[it].force[ia].x = fd_forces[iat * 3 + 0];
            ucell.atoms[it].force[ia].y = fd_forces[iat * 3 + 1];
            ucell.atoms[it].force[ia].z = fd_forces[iat * 3 + 2];
            iat++;
        }
    }
}