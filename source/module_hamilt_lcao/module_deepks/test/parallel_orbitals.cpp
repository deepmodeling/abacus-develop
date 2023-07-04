#include "parallel_orbitals.h"

namespace Test_Deepks
{

    Parallel_Orbitals::Parallel_Orbitals()
    {
        global2local_row = nullptr;
        global2local_col = nullptr;
    }

    Parallel_Orbitals::~Parallel_Orbitals()
    {
        delete[] global2local_row;
        delete[] global2local_col;
    }

    void Parallel_Orbitals::set_global2local(void)
    {
        ModuleBase::TITLE("Parallel_Orbitals","set_trace");
        assert(GlobalV::NLOCAL>0);

        delete[] global2local_row;
        delete[] global2local_col;

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"global2local_row dimension",GlobalV::NLOCAL);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"global2local_col dimension",GlobalV::NLOCAL);

        global2local_row = new int[GlobalV::NLOCAL];
        global2local_col = new int[GlobalV::NLOCAL];
        // mohan update 2011-04-07
        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
            global2local_row[i] = -1;
            global2local_col[i] = -1;
        }

        ModuleBase::Memory::record("PO::global2local_row",sizeof(int) * GlobalV::NLOCAL);
        ModuleBase::Memory::record("PO::global2local_col",sizeof(int) * GlobalV::NLOCAL);

        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            global2local_row[i] = i;
            global2local_col[i] = i;
        }
        this->nrow = GlobalV::NLOCAL;
        this->ncol = GlobalV::NLOCAL;
        this->nloc=nrow*ncol;

        return;
    }
}