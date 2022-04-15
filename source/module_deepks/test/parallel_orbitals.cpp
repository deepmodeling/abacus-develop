#include "parallel_orbitals.h"

namespace Test_Deepks
{

    Parallel_Orbitals::Parallel_Orbitals()
    {
        trace_loc_row = new int[1];
        trace_loc_col = new int[1];
    }

    Parallel_Orbitals::~Parallel_Orbitals()
    {
        delete[] trace_loc_row;
        delete[] trace_loc_col;
    }

    void Parallel_Orbitals::set_trace(void)
    {
        ModuleBase::TITLE("Parallel_Orbitals","set_trace");
        assert(GlobalV::NLOCAL>0);

        delete[] trace_loc_row;
        delete[] trace_loc_col;

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"trace_loc_row dimension",GlobalV::NLOCAL);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"trace_loc_col dimension",GlobalV::NLOCAL);

        trace_loc_row = new int[GlobalV::NLOCAL];
        trace_loc_col = new int[GlobalV::NLOCAL];
        // mohan update 2011-04-07
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
            trace_loc_row[i] = -1;
            trace_loc_col[i] = -1;
        }

        ModuleBase::Memory::record("Parallel_Orbitals","trace_loc_row",GlobalV::NLOCAL,"int");
        ModuleBase::Memory::record("Parallel_Orbitals","trace_loc_col",GlobalV::NLOCAL,"int");

        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            trace_loc_row[i] = i;
            trace_loc_col[i] = i;
        }
        this->nrow = GlobalV::NLOCAL;
        this->ncol = GlobalV::NLOCAL;
        this->nloc=nrow*ncol;

        return;
    }
}