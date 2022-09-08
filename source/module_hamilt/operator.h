#ifndef __OPERATOR
#define __OPERATOR

#include<complex>
#include "module_psi/psi.h"
#include "module_base/global_function.h"
#include "module_base/tool_quit.h"

namespace hamilt
{

// Basic class for operator module, 
// it is designed for "O|psi>" and "<psi|O|psi>"
// Operator "O" might have several different types, which should be calculated one by one.
// In basic class , function add() is designed for combine all operators together with a chain. 
template<typename T>
class Operator
{
    public:
    Operator(){};
    virtual ~Operator()
    { 
        if(this->hpsi != nullptr) delete this->hpsi;
        Operator* last = this->next_op;
        while(last != nullptr) 
        {
            Operator* node_delete = last;
            last = last->next_op;
            node_delete->next_op = nullptr;
            delete node_delete;
        } 
    }

    //this is the core function for Operator
    // do H|psi> from input |psi> , 
    // output of hpsi would be first member of the returned tuple 
    typedef std::tuple<const psi::Psi<T>*, const psi::Range, T*> hpsi_info;
    virtual hpsi_info hPsi(hpsi_info& input)const 
    {
        ModuleBase::WARNING_QUIT("Operator::hPsi", "hPsi error!");
        return hpsi_info(nullptr, 0, nullptr);
    }

    virtual void init(const int ik_in)
    {
        this->ik = ik_in;
        if(this->next_op != nullptr)
        {
            this->next_op->init(ik_in);
        }
    }

    virtual void add(Operator* next)
    {
        if(next==nullptr) return;
        if(next->next_op != nullptr) this->add(next->next_op);
        Operator* last = this;
        while(last->next_op != nullptr)
        {
            if(next->cal_type==last->cal_type)
            {
                last->add(next);
                return;
            }
            last = last->next_op;
        }
        last->next_op = next; 
    }

    protected:
    int ik = 0;

    mutable bool in_place = false;

    //calculation type, only different type can be in main chain table 
    int cal_type = 0;
    Operator* next_op = nullptr;

    //if this Operator is first node in chain table, hpsi would not be empty
    mutable psi::Psi<T>* hpsi = nullptr;

    /*This function would analyze hpsi_info and choose how to arrange hpsi storage
    In hpsi_info, if the third parameter hpsi_pointer is set, which indicates memory of hpsi is arranged by developer;
    if hpsi_pointer is not set(nullptr), which indicates memory of hpsi is arranged by Operator, this case is rare. 
    two cases would occurred:
    1. hpsi_pointer != nullptr && psi_pointer == hpsi_pointer , psi would be replaced by hpsi, hpsi need a temporary memory
    2. hpsi_pointer != nullptr && psi_pointer != hpsi_pointer , this is the commonly case 
    */
    T* get_hpsi(const hpsi_info& info)const
    {
        const int nbands_range = (std::get<1>(info).range_2 - std::get<1>(info).range_1 + 1);
        //in_place call of hPsi, hpsi inputs as new psi, 
        //create a new hpsi and delete old hpsi later
        T* hpsi_pointer = std::get<2>(info);
        const T* psi_pointer = std::get<0>(info)->get_pointer();
        if(this->hpsi != nullptr) 
        {
            delete this->hpsi;
            this->hpsi = nullptr;
        }
        if(!hpsi_pointer)
        {
            ModuleBase::WARNING_QUIT("Operator::hPsi", "hpsi_pointer can not be nullptr");
        }
        else if(hpsi_pointer == psi_pointer)
        {
            this->in_place = true;
            this->hpsi = new psi::Psi<T>(std::get<0>(info)[0], 1, nbands_range);
        }
        else
        {
            this->in_place = false;
            this->hpsi = new psi::Psi<T>(hpsi_pointer, std::get<0>(info)[0], 1, nbands_range);
        }
        
        hpsi_pointer = this->hpsi->get_pointer();
        size_t total_hpsi_size = nbands_range * this->hpsi->get_nbasis();
        ModuleBase::GlobalFunc::ZEROS(hpsi_pointer, total_hpsi_size);
        return hpsi_pointer;
    }
};

}//end namespace hamilt

#endif