#include "dspin_lcao.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_base/blas_connector.h"
#include "module_base/timer.h"
#include "module_base/memory.h"
#include "module_base/tool_title.h"
#include "module_base/parallel_reduce.h"
#include "module_parameter/parameter.h"

template <typename TK, typename TR>
hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, TR>>::DeltaSpin(
    HS_Matrix_K<TK>* hsk_in,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
    hamilt::HContainer<TR>* hR_in,
    const UnitCell& ucell_in,
    Grid_Driver* gridD_in,
    const TwoCenterIntegrator* intor,
    const std::vector<double>& orb_cutoff)
    : hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in), intor_(intor), orb_cutoff_(orb_cutoff)
{
    this->cal_type = calculation_type::lcao_sc_lambda;
    this->ucell = &ucell_in;
    this->gridD = gridD_in;
#ifdef __DEBUG
    assert(this->ucell != nullptr);
    assert(this->gridD != nullptr);
#endif
    //set nspin
    this->nspin = PARAM.inp.nspin;
    this->spin_num = this->nspin == 2 ? 2 : 1;

    this->lambda_save.resize(this->ucell->nat * 3, 0.0);
    this->update_lambda_.resize(this->nspin, false);
}

// destructor
template <typename TK, typename TR>
hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, TR>>::~DeltaSpin()
{
    for (auto& hr : this->pre_hr)
    {
        if (hr != nullptr)
        {
            delete hr;
            hr = nullptr;
        }
    }
    this->pre_hr.clear();
    this->pre_hr.shrink_to_fit();
}

// simple functions to calculate the coefficients from lambda
inline void cal_coeff_lambda(const std::vector<double>& current_lambda, std::vector<double>& coefficients)
{
    coefficients[0] = current_lambda[0];
    coefficients[1] = -current_lambda[0];
}
inline void cal_coeff_lambda(const std::vector<double>& current_lambda, std::vector<std::complex<double>>& coefficients)
{// {\lambda^{I,3}, \lambda^{I,1}+i\lambda^{I,2}, \lambda^{I,1}-i\lambda^{I,2}, -\lambda^{I,3}}
    coefficients[0] = std::complex<double>(current_lambda[2], 0.0);
    coefficients[1] = std::complex<double>(current_lambda[0] , current_lambda[1]);
    coefficients[2] = std::complex<double>(current_lambda[0] , -1 * current_lambda[1]);
    coefficients[3] = std::complex<double>(-1 * current_lambda[2], 0.0);
}

template <typename TK, typename TR>
void hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    // if lambda has not changed, calculate the HR^I = lambda^I\sum_{lm}<phi_mu|alpha^I_{lm}><alpha^I_{lm}|phi_{nu,R}>
    // if lambda has changed, calculate the dHR^I = dlambda^I\sum_{lm}<phi_mu|alpha^I_{lm}><alpha^I_{lm}|phi_{nu,R}> 
    spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
    // there are three case for contributeHR 
    // 1. HR has not been calculated, reset lambda_save and calculate HR = lambda * pre_hr
    // 2. HR has been calculated, but lambda has changed, calculate dHR = dlambda * pre_hr
    // 3. HR has been calculated, and lambda has not changed, do nothing
    if(!this->hr_done)
    {
        // set the lambda_save to zero if lambda loop is started
        this->lambda_save.assign(this->ucell->nat * 3, 0.0);
    }
    else if(this->hr_done && !this->update_lambda_[this->current_spin])
    {
        return;
    }

    // calculate Hpre^I = \sum_{lm}<phi_mu|alpha^I_{lm}><alpha^I_{lm}|phi_{nu,R}>
    if(!this->initialized)
    {
        auto& constrain = sc.get_constrain();
        this->cal_constraint_atom_list(constrain);
        this->cal_pre_HR();
        this->initialized = true;
    }
    auto& lambda = sc.get_sc_lambda();
    
    for(int iat=0;iat<this->ucell->nat;iat++)
    {
        if(!this->constraint_atom_list[iat])
        {
            continue;
        }
        // calculate the delta lambda to update the real space Hamiltonian
        std::vector<double> current_lambda;
        if(this->nspin==4)
        {
            current_lambda = {lambda[iat].x - this->lambda_save[iat*3], 
            lambda[iat].y - this->lambda_save[iat*3+1], 
            lambda[iat].z - this->lambda_save[iat*3+2]};
        }
        else if(this->nspin==2)
        {
            current_lambda = {lambda[iat].z-this->lambda_save[iat*3+2], 0.0, 0.0};
        }
        std::vector<TR> coefficients(this->nspin);

        cal_coeff_lambda(current_lambda, coefficients);

        // magnetic moment = \sum_{\mu\nu,R} dmR * pre_hr
        for(int iap=0;iap<this->pre_hr[iat]->size_atom_pairs();iap++)
        {
            hamilt::AtomPair<TR>& tmp = this->pre_hr[iat]->get_atom_pair(iap);
            int iat1 = tmp.get_atom_i();
            int iat2 = tmp.get_atom_j();
            int row_size = tmp.get_row_size();
            int col_size = tmp.get_col_size();
            if(this->nspin==4)
            {
                this->pre_coeff_array(coefficients, row_size, col_size);
            }
            for(int ir = 0;ir < tmp.get_R_size(); ++ir )
            {
                const ModuleBase::Vector3<int> r_index = tmp.get_R_index(ir);
                const TR* pre_hr_data = tmp.get_pointer(ir);
                TR* dhr_data = this->hR->find_matrix(iat1, iat2, r_index[0], r_index[1], r_index[2])->get_pointer();
                // TR== double: axpy for dhr_data += current_lambda * pre_hr_data
                // TR!= double: call cal_lambda_hr_IJR
                if (this->nspin==2)
                {
                    //BlasConnector::axpy(row_size*col_size, coefficients[this->current_spin], pre_hr_data, dhr_data);
                    for(int i=0;i<tmp.get_size();i++)
                    {
                        dhr_data[i] += pre_hr_data[i] * coefficients[this->current_spin];
                    }
                }
                else
                {
                    for(int i=0;i<tmp.get_size();i++)
                    {
                        dhr_data[i] += pre_hr_data[i] * this->tmp_coeff_array[i];
                    }
                }
            }
        }
    }

    // save lambda to lambda_save or update the current_spin in NSPIN=2
    this->update_lambda_[this->current_spin] = false;
    if(this->current_spin == this->spin_num - 1)
    {
        for(int i=0;i<this->ucell->nat;i++)
        {
            if(this->constraint_atom_list[i])
            {
                for(int j=0;j<3;j++)
                {   
                    this->lambda_save[i*3+j] = lambda[i][j];
                }
            }
        }
    }
    return;
}

// cal_lambda_hr_IJR
template <typename TK, typename TR>
void hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, TR>>::pre_coeff_array(
    const std::vector<TR>& coeff, const int row_size, const int col_size)
{
    this->tmp_coeff_array.resize(row_size*col_size);
    for(int irow=0;irow<row_size;irow+=2)
    {
        for(int icol=0;icol<col_size;icol+=2)
        {
            this->tmp_coeff_array[irow*col_size+icol] = coeff[0];
            this->tmp_coeff_array[irow*col_size+icol+1] = coeff[1];
            this->tmp_coeff_array[(irow+1)*col_size+icol] = coeff[2];
            this->tmp_coeff_array[(irow+1)*col_size+icol+1] = coeff[3];
        }
    }
}

// cal_constraint_atom_list()
template <typename TK, typename TR>
void hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, TR>>::cal_constraint_atom_list(const std::vector<ModuleBase::Vector3<int>>& constraints)
{
    this->constraint_atom_list.clear();
    this->constraint_atom_list.resize(this->ucell->nat, false);
#ifdef __DEBUG
    assert(this->ucell->nat == constraints.size());
#endif
    for(int iat=0;iat<this->ucell->nat;iat++)
    {
        if(constraints[iat][0] + constraints[iat][1] + constraints[iat][2] == 0)
        {
            this->constraint_atom_list[iat] = false;
        }
        else
        {
            this->constraint_atom_list[iat] = true;
        }
    }
}

// cal_pre_HR()
template <typename TK, typename TR>
void hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, TR>>::cal_pre_HR()
{
    ModuleBase::TITLE("DeltaSpin", "cal_pre_HR");
    if(this->initialized)
    {
        return;
    }
    this->paraV = this->hR->get_paraV();
    ModuleBase::timer::tick("DeltaSpin", "cal_pre_HR");
    this->pre_hr.clear();
    this->pre_hr.resize(this->ucell->nat, nullptr);

    const int npol = this->ucell->get_npol();
    size_t memory_cost = 0;
    for(int iat=0;iat<this->ucell->nat;iat++)
    {
        if(!this->constraint_atom_list[iat])
        {
            continue;
        }
        
        auto tau0 = ucell->get_tau(iat);
        int T0, I0;
        this->ucell->iat2iait(iat, &I0, &T0);

        // first step: find the adjacent atoms and filter the real adjacent atoms
        AdjacentAtomInfo adjs;
        this->gridD->Find_atom(*ucell, tau0, T0, I0, &adjs);

        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat, iat1, R_index1).norm() * this->ucell->lat0
                < this->orb_cutoff_[T1] + PARAM.inp.onsite_radius)
            {
                is_adj[ad] = true;
            }
        }
        filter_adjs(is_adj, adjs);

        // second step: prepare the <IJR> pairs for each iat-atom
        this->pre_hr[iat] = new hamilt::HContainer<TR>(this->paraV);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = this->ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
            {
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int iat2 = this->ucell->itia2iat(T2, I2);
                ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                int r_vector[3] = {R_index2.x - R_index1.x, R_index2.y - R_index1.y, R_index2.z - R_index1.z};
                // keep the size of pre_hr for each atom less than this->hR
                if(this->hR->find_matrix(iat1, iat2, r_vector[0], r_vector[1], r_vector[2]) == nullptr)
                {
                    continue;
                }
                hamilt::AtomPair<TR> tmp(iat1,
                                         iat2,
                                         r_vector[0],
                                         r_vector[1],
                                         r_vector[2],
                                         paraV);
                this->pre_hr[iat]->insert_pair(tmp);
            }
        }
        this->pre_hr[iat]->allocate(nullptr, true);

        // third step: calculate the <phi|alpha> overlap integrals
        const int max_l_plus_1 = this->ucell->atoms[T0].nwl + 1;
        std::vector<std::unordered_map<int, std::vector<double>>> nlm_iat0(adjs.adj_num + 1);
        for(int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int iat1 = this->ucell->itia2iat(T1, I1);
            const Atom* atom1 = &ucell->atoms[T1];
            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];

            auto all_indexes = paraV->get_indexes_row(iat1);
            auto col_indexes = paraV->get_indexes_col(iat1);
            // insert col_indexes into all_indexes to get universal set with no repeat elements
            all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
            std::sort(all_indexes.begin(), all_indexes.end());
            all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());
            for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
            {
                const int iw1 = all_indexes[iw1l] / npol;
                // only first zeta orbitals in target L of atom iat0 are needed
                std::vector<double> nlm_target(max_l_plus_1 * max_l_plus_1);
                const int L1 = atom1->iw2l[ iw1 ];
                const int N1 = atom1->iw2n[ iw1 ];
                const int m1 = atom1->iw2m[ iw1 ];
                std::vector<std::vector<double>> nlm;
                // nlm is a vector of vectors, but size of outer vector is only 1 here
                // If we are calculating force, we need also to store the gradient
                // and size of outer vector is then 4
                // inner loop : all projectors (L0,M0)

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;
                ModuleBase::Vector3<double> dtau = tau0 - tau1;
                intor_->snap(T1, L1, N1, M1, T0, dtau * this->ucell->lat0, 0 /*cal_deri*/, nlm);

                // select the elements of nlm with target_L (0, 1, 2, 3 ...)
                int target_L = 0, index=0;
                for(int iw =0;iw < this->ucell->atoms[T0].nw; iw++)
                {
                    const int L0 = this->ucell->atoms[T0].iw2l[iw];
                    // only the first zeta of each l-orbital is needed
                    if(L0 == target_L)
                    {
                        for(int m = 0; m < 2*L0+1; m++)
                        {
                            nlm_target[index] = nlm[0][iw+m];
                            index++;
                        }
                        target_L++;
                    }
                }
                nlm_iat0[ad].insert({all_indexes[iw1l], nlm_target});
            }
        }

        // fourth step: calculate the <phi|alpha><alpha|phi>
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            const std::unordered_map<int,std::vector<double>>& nlm1 = nlm_iat0[ad1];
            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
            {
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int iat2 = ucell->itia2iat(T2, I2);
                const std::unordered_map<int,std::vector<double>>& nlm2 = nlm_iat0[ad2];
                ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                                  R_index2[1] - R_index1[1],
                                                  R_index2[2] - R_index1[2]);
                hamilt::BaseMatrix<TR>* tmp = this->pre_hr[iat]->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                // if not found , skip this pair of atoms
                if (tmp != nullptr)
                {
                    this->cal_HR_IJR(iat1, iat2, nlm1, nlm2, tmp->get_pointer());
                }
            }
        }
        memory_cost += this->pre_hr[iat]->get_memory_size();
    }
    ModuleBase::Memory::record("DeltaSpin:pre_HR", memory_cost);
    ModuleBase::timer::tick("DeltaSpin", "cal_pre_HR");
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(const int& iat1,
                    const int& iat2,
                    const std::unordered_map<int, std::vector<double>>& nlm1_all,
                    const std::unordered_map<int, std::vector<double>>& nlm2_all,
                    TR* data_pointer)
{
    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = this->paraV->get_indexes_row(iat1);
    auto col_indexes = this->paraV->get_indexes_col(iat2);
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol*npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = this->paraV->get_col_size(iat2) * is + is2;
        }
    }
    // calculate the local matrix
    const TR* tmp_d = nullptr;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            TR nlm_tmp = TR(0);
            for (int m1 = 0; m1 < nlm1.size(); m1++)
            {
                nlm_tmp += nlm1[m1] * nlm2[m1];
            }
            for (int is = 0; is < npol*npol; ++is)
            {
                data_pointer[step_trace[is]] += nlm_tmp;
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

// cal_moment
template <typename TK, typename TR>
std::vector<double> hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, TR>>::cal_moment(const HContainer<double>* dmR, const std::vector<ModuleBase::Vector3<int>>& constrain)
{
    const int mag_fold = this->nspin==4?3:1;
    std::vector<double> moment(this->ucell->nat * mag_fold, 0.0);
    if(dmR == nullptr)
    {
        return moment;
    }
    if (!this->initialized)
    {
        //spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        //auto& constrain = sc.get_constrain();
        this->cal_constraint_atom_list(constrain);
        this->cal_pre_HR();
        this->initialized = true;
    }
    for(int iat=0;iat<this->ucell->nat;iat++)
    {
        if(!this->constraint_atom_list[iat])
        {
            continue;
        }
        // magnetic moment = \sum_{\mu\nu,R} dmR * pre_hr
        for(int iap=0;iap<this->pre_hr[iat]->size_atom_pairs();iap++)
        {
            hamilt::AtomPair<TR>& tmp = this->pre_hr[iat]->get_atom_pair(iap);
            int iat1 = tmp.get_atom_i();
            int iat2 = tmp.get_atom_j();
            int row_size = tmp.get_row_size();
            int col_size = tmp.get_col_size();
            for(int ir = 0;ir < tmp.get_R_size(); ++ir )
            {
                const ModuleBase::Vector3<int> r_index = tmp.get_R_index(ir);
                double* dmr_data = dmR->find_matrix(iat1, iat2, r_index[0], r_index[1], r_index[2])->get_pointer();
                const TR* hr_data = tmp.get_pointer(ir);
                this->cal_moment_IJR(dmr_data, hr_data, row_size, col_size, &moment[iat*mag_fold]);
            }
        }
    }
#ifdef __MPI
    // sum up the magnetic moments
    Parallel_Reduce::reduce_all(moment.data(), moment.size());
#endif 
    return moment;
}

// cal_moment_IJR
template <typename TK, typename TR>
void hamilt::DeltaSpin<hamilt::OperatorLCAO<TK, TR>>::cal_moment_IJR(
    const double* dmR, 
    const TR* hr, 
    const int row_size,
    const int col_size, 
    double* moment
)
{
    // collinear spin case
    TR tmp_moment = TR(0);
    for(int i=0;i<row_size*col_size;i++)
    {
        tmp_moment += dmR[i] * hr[i];
    }
    moment[0] += tmp_moment;
}

template<>
void hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>::cal_moment_IJR(
    const double* dmR, 
    const std::complex<double>* hr, 
    const int row_size,
    const int col_size, 
    double* moment
)
{
    const int step_trace[4] = {0, 1, col_size, col_size+1};
    int index = 0;
    std::vector<std::complex<double>> tmp_moment(3, std::complex<double>(0.0, 0.0));
    for(int irow=0;irow<row_size;irow+=2)
    {
        for(int icol=0;icol<col_size;icol+=2)
        {
            tmp_moment[0] += dmR[index+step_trace[1]] * hr[index];
            tmp_moment[1] += dmR[index+step_trace[2]] * hr[index];
            tmp_moment[2] += dmR[index+step_trace[3]] * hr[index];
            index += 2;
        }
        index += col_size;
    }
    moment[0] += tmp_moment[0].real();
    moment[1] += tmp_moment[1].real();
    moment[2] += tmp_moment[2].real();
}

#include "dspin_force_stress.hpp"

template class hamilt::DeltaSpin<hamilt::OperatorLCAO<double, double>>;
template class hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;