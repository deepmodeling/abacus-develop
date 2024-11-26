#include "operator_lcao.h"

#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_hsolver/hsolver_lcao.h"

#include "module_parameter/parameter.h"

#ifdef __ELPA
#include "module_hsolver/diago_elpa.h"
#include "module_hsolver/diago_elpa_native.h"
#endif

namespace hamilt {

template <>
void OperatorLCAO<double, double>::get_hs_pointers() {
    ModuleBase::timer::tick("OperatorLCAO", "get_hs_pointers");
    this->hmatrix_k = this->hsk->get_hk();
    if ((this->new_e_iteration && ik == 0) || PARAM.inp.out_mat_hs[0])
    {
        if (this->smatrix_k == nullptr)
        {
            this->smatrix_k = new double[this->hsk->get_size()];
            this->allocated_smatrix = true;
        }
        const int inc = 1;
        BlasConnector::copy(this->hsk->get_size(), this->hsk->get_sk(), inc, this->smatrix_k, inc);
#ifdef __ELPA
        hsolver::DiagoElpa<double>::DecomposedState = 0;
        hsolver::DiagoElpaNative<double>::DecomposedState = 0;
#endif
        this->new_e_iteration = false;
    }
    ModuleBase::timer::tick("OperatorLCAO", "get_hs_pointers");
}

template<>
void OperatorLCAO<std::complex<double>, double>::get_hs_pointers()
{
    this->hmatrix_k = this->hsk->get_hk();
    this->smatrix_k = this->hsk->get_sk();
}

template<>
void OperatorLCAO<std::complex<double>, std::complex<double>>::get_hs_pointers()
{
    this->hmatrix_k = this->hsk->get_hk();
    this->smatrix_k = this->hsk->get_sk();
}

template<typename TK, typename TR>
void OperatorLCAO<TK, TR>::refresh_h()
{
    // Set the matrix 'H' to zero.
    this->hsk->set_zero_hk();
}

template <typename TK, typename TR>
void OperatorLCAO<TK, TR>::set_hr_done(bool hr_done_in) {
    this->hr_done = hr_done_in;
}

template<typename TK, typename TR>
void OperatorLCAO<TK, TR>::set_current_spin(const int current_spin_in)
{
    this->current_spin = current_spin_in;
    if(this->next_op != nullptr)
    {
        dynamic_cast<OperatorLCAO<TK, TR>*>(this->next_op)->set_current_spin(current_spin_in);
    }
    if(this->next_sub_op != nullptr)
    {
        dynamic_cast<OperatorLCAO<TK, TR>*>(this->next_sub_op)->set_current_spin(current_spin_in);
    }
}

template <typename TK, typename TR>
void OperatorLCAO<TK, TR>::init(const int ik_in) {
    ModuleBase::TITLE("OperatorLCAO", "init");
    ModuleBase::timer::tick("OperatorLCAO", "init");
    if (this->is_first_node) {
        // refresh HK
        this->refresh_h();
        if (!this->hr_done) {
            // refresh HR
            this->hR->set_zero();
        }
    }
    switch (this->cal_type) {
    case calculation_type::lcao_overlap: {
        // cal_type=lcao_overlap refer to overlap matrix operators, which are
        // only rely on stucture, and not changed during SCF

        if (!this->hr_done) {
            // update SR first
            // in cal_type=lcao_overlap, SR should be updated by each sub-chain
            // nodes
            OperatorLCAO<TK, TR>* last = this;
            while (last != nullptr) {
                last->contributeHR();
                last = dynamic_cast<OperatorLCAO<TK, TR>*>(last->next_sub_op);
            }
        }

        // update SK next
        // in cal_type=lcao_overlap, SK should be update here
        this->contributeHk(ik_in);

        break;
    }
    case calculation_type::lcao_fixed: {
        // cal_type=lcao_fixed refer to fixed matrix operators, which are only
        // rely on stucture, and not changed during SCF

        // update HR first
        if (!this->hr_done) {
            // in cal_type=lcao_fixed, HR should be updated by each sub-chain
            // nodes
            OperatorLCAO<TK, TR>* last = this;
            while (last != nullptr) {
                last->contributeHR();
                last = dynamic_cast<OperatorLCAO<TK, TR>*>(last->next_sub_op);
            }
        }

        // update HK next
        // in cal_type=lcao_fixed, HK will update in the last node with
        // OperatorLCAO::contributeHk()

        break;
    }
    case calculation_type::lcao_gint: {
        // cal_type=lcao_gint refer to grid integral operators, which are relied
        // on stucture and potential based on real space grids and should be
        // updated each SCF steps

        if (!this->hr_done) {
            OperatorLCAO<TK, TR>* last = this;
            while (last != nullptr) {
                // update HR first
                // in cal_type=lcao_gint, HR should be updated by every
                // sub-node.
                last->contributeHR();

                // update HK next
                // in cal_type=lcao_gint, HK will update in the last node with
                // OperatorLCAO::contributeHk()
                last = dynamic_cast<OperatorLCAO<TK, TR>*>(last->next_sub_op);
            }
        }

        break;
    }
#ifdef __DEEPKS
    case calculation_type::lcao_deepks: {
        // update HR first
        if (!this->hr_done) {
            // in cal_type=lcao_deepks, HR should be updated
            this->contributeHR();
        }

        // update H_V_delta_k next
        this->contributeHk(ik_in);

        break;
    }
#endif
        case calculation_type::lcao_dftu:
        {
            //only HK should be updated when cal_type=lcao_dftu
            //in cal_type=lcao_dftu, HK only need to update from one node
            if(!this->hr_done)
            {
                //in cal_type=lcao_deepks, HR should be updated
                this->contributeHR();
            }
            break;
        }
        case calculation_type::lcao_sc_lambda:
        {
            //update HR first
            this->contributeHR();
            //in cal_type=lcao_sc_mag, 
            //this->contributeHk(ik_in);
            break;
        }
        case calculation_type::lcao_exx:
        {
            //update HR first
            if (!this->hr_done)
            {
                this->contributeHR();
            }

            //update HK next
            //in cal_type=lcao_exx, HK only need to update from one node
            // this->contributeHk(ik_in);

        break;
    }
    case calculation_type::lcao_tddft_velocity: {
        if (!this->hr_done) {
            // in cal_type=lcao_fixed, HR should be updated by each sub-chain
            // nodes
            OperatorLCAO<TK, TR>* last = this;
            while (last != nullptr) {
                last->contributeHR();
                last = dynamic_cast<OperatorLCAO<TK, TR>*>(last->next_sub_op);
            }
        }
        this->contributeHk(ik_in);

        break;
    }
    default: {
        ModuleBase::WARNING_QUIT("OperatorLCAO::init", "unknown cal_type");
        break;
    }
    }
    if (this->next_op
        != nullptr) { // it is not the last node, loop next init() function
        // pass HR status to next node and than set HR status of this node to
        // done
        if (!this->hr_done) {
            dynamic_cast<OperatorLCAO<TK, TR>*>(this->next_op)->hr_done
                = this->hr_done;
        }
        // call init() function of next node
        this->next_op->init(ik_in);
    } else { // it is the last node, update HK with the current total HR
        OperatorLCAO<TK, TR>::contributeHk(ik_in);
    }

    // set HR status of this node to done
    this->hr_done = true;

    ModuleBase::timer::tick("OperatorLCAO", "init");
}

// contributeHk()
template <typename TK, typename TR>
void OperatorLCAO<TK, TR>::contributeHk(int ik) {
    ModuleBase::TITLE("OperatorLCAO", "contributeHk");
    ModuleBase::timer::tick("OperatorLCAO", "contributeHk");
    if(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
    {
        const int nrow = this->hsk->get_pv()->get_row_size();
        hamilt::folding_HR(*this->hR, this->hsk->get_hk(), this->kvec_d[ik], nrow, 1);
    }
    else
    {
        const int ncol = this->hsk->get_pv()->get_col_size();
        hamilt::folding_HR(*this->hR, this->hsk->get_hk(), this->kvec_d[ik], ncol, 0);
    }
    ModuleBase::timer::tick("OperatorLCAO", "contributeHk");
}

template class OperatorLCAO<double, double>;
template class OperatorLCAO<std::complex<double>, double>;
template class OperatorLCAO<std::complex<double>, std::complex<double>>;
} // namespace hamilt
