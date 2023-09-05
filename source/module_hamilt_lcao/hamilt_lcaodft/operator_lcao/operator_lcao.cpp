#include "operator_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#ifdef __ELPA
#include "module_hsolver/diago_elpa.h"
#endif

namespace hamilt
{

template<>
void OperatorLCAO<double>::get_hs_pointers()
{
    ModuleBase::timer::tick("OperatorLCAO", "get_hs_pointers");
    this->hmatrix_k = this->LM->Hloc.data();
    if ((this->new_e_iteration && ik == 0) || hsolver::HSolverLCAO::out_mat_hs)
    {
        if (this->smatrix_k == nullptr)
        {
            this->smatrix_k = new double[this->LM->Sloc.size()];
            this->allocated_smatrix = true;
        }
        const int inc = 1;
        BlasConnector::copy(this->LM->Sloc.size(), this->LM->Sloc.data(), inc, this->smatrix_k, inc);
#ifdef __ELPA
        hsolver::DiagoElpa::DecomposedState = 0;
#endif
        this->new_e_iteration = false;
    }
    ModuleBase::timer::tick("OperatorLCAO", "get_hs_pointers");
}

template<>
void OperatorLCAO<std::complex<double>>::get_hs_pointers()
{
    this->hmatrix_k = this->LM->Hloc2.data();
    this->smatrix_k = this->LM->Sloc2.data();
}

template<>
void OperatorLCAO<double>::refresh_h()
{
    // Set the matrix 'H' to zero.
    this->LM->zeros_HSgamma('H');
}

template<>
void OperatorLCAO<std::complex<double>>::refresh_h()
{
    // Set the matrix 'H' to zero.
    this->LM->zeros_HSk('H');
}

template<>
void OperatorLCAO<double>::folding_fixed(const int ik, const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("OperatorLCAO", "folding_fixed");
    ModuleBase::timer::tick("OperatorLCAO", "folding_fixed");
    //add T+VNL matrix.
	this->LM->update_Hloc();
    ModuleBase::timer::tick("OperatorLCAO", "folding_fixed");
}

template<>
void OperatorLCAO<std::complex<double>>::folding_fixed(const int ik, const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("OperatorLCAO", "folding_fixed");
    ModuleBase::timer::tick("OperatorLCAO", "folding_fixed");
    //-----------------------------------------
    // folding matrix here: S(k) (SlocR->Sloc2)
    // folding matrix here: T(k)+Vnl(k)
    // (Hloc_fixed->Hloc_fixed2)
    //-----------------------------------------
    this->LM->zeros_HSk('S');
    this->LM->zeros_HSk('T');
    this->LM->folding_fixedH(ik, kvec_d);

    //------------------------------------------
    // Add T(k)+Vnl(k)+Vlocal(k)
    // (Hloc2 += Hloc_fixed2), (std::complex matrix)
    //------------------------------------------
	this->LM->update_Hloc2(ik);
    ModuleBase::timer::tick("OperatorLCAO", "folding_fixed");
}

template<typename T>
void OperatorLCAO<T>::init(const int ik_in)
{
    ModuleBase::TITLE("OperatorLCAO", "init");
    ModuleBase::timer::tick("OperatorLCAO", "init");
    if(this->is_first_node)
    {
        this->refresh_h();
    }
    switch(this->cal_type)
    {
        case lcao_fixed:
        {
            //cal_type=lcao_fixed refer to fixed matrix operators, which are only rely on stucture, and not changed during SCF

            //update HR first
            //in cal_type=lcao_fixed, HR should be updated by each sub-chain nodes
            OperatorLCAO<T>* last = this;
            while(last != nullptr)
            {
                last->contributeHR();
                last = dynamic_cast<OperatorLCAO<T>*>(last->next_sub_op);
            }

            //update HK next
            //in cal_type=lcao_fixed, HK only need to update together in folding_fixed()

            break;
        }
        case lcao_gint:
        {
            //cal_type=lcao_gint refer to grid integral operators, which are relied on stucture and potential based on real space grids
            //and should be updated each SCF steps

            OperatorLCAO<T>* last = this;
            while(last != nullptr)
            {
                //update HR first
                //in cal_type=lcao_gint, HR should be updated by every sub-node.
                last->contributeHR();

                //update HK next
                //in cal_type=lcao_gint, HK should be updated by every sub-node.
                last->contributeHk(ik_in);

                last = dynamic_cast<OperatorLCAO<T>*>(last->next_sub_op);
            }

            break;
        }
        case lcao_deepks:
        {
            //update HR first
            //in cal_type=lcao_deepks, HR should be updated
            this->contributeHR();

            //update HK next
            //in cal_type=lcao_deepks, HK should be updated
            this->contributeHk(ik_in);

            break;

        }
        case lcao_dftu:
        {
            //only HK should be updated when cal_type=lcao_dftu
            //in cal_type=lcao_dftu, HK only need to update from one node
            //dftu is a special operator, it should be the last node of chain
            //Overlap matrix for ik is used by it, do folding first and then return
            this->folding_fixed(ik_in, this->kvec_d);
            this->contributeHk(ik_in);
            break;
        }
        case lcao_sc_lambda:
        {
            //update HK only
            //in cal_type=lcao_sc_mag, HK only need to update from one node
            // folding fixed should only be called once
            // with folding_fixed here, the lcao_dftu should not be invoked
            // only by comment out the folding_fixed here, the lcao_dftu could be invoked
            this->folding_fixed(ik_in, this->kvec_d);
            this->contributeHk(ik_in);
            return;
            break;
        }
        case lcao_exx:
        {
            //update HR first
            //in cal_type=lcao_exx, HR should be updated by most priority sub-chain nodes
            this->contributeHR();

            //update HK next
            //in cal_type=lcao_exx, HK only need to update from one node
            this->contributeHk(ik_in);

            break;
        }
        default:
        {
            ModuleBase::WARNING_QUIT("OperatorLCAO::init", "unknown cal_type");
            break;
        }
    }
    if(this->next_op != nullptr)
    {//it is not the last node, loop next init() function
        this->next_op->init(ik_in);
    }
    else
    {//it is the last node, do folding process together
        this->folding_fixed(ik_in, this->kvec_d);
    }

    ModuleBase::timer::tick("OperatorLCAO", "init");
}
template class OperatorLCAO<double>;
template class OperatorLCAO<std::complex<double>>;
}  // namespace hamilt
