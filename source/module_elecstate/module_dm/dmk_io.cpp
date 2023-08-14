#include "dmk_io.h"
#include "module_base/timer.h"

// output the density matrix in k-space
// weiqing add 2023/8/9
void elecstate::write_dmk(
    const K_Vectors& kv, // k-vectors
    const int& ik, // index of k-vector
    const int& nlocal, // number of local orbitals
    const std::string &fn, // file name
    std::vector<ModuleBase::ComplexMatrix> &dm_k)
{
    ModuleBase::TITLE("elecstate","write_dmk");
    ModuleBase::timer::tick("elecstate","write_dmk");

    time_t start, end;
    std::ofstream ofs;

    start = time(NULL);
    ofs.open(fn.c_str());
    if (!ofs)
    {
        ModuleBase::WARNING("elecstate::write_dmk","Can't create DENSITY MATRIX File!");
    }
    ofs << kv.kvec_d[ik].x << " " << kv.kvec_d[ik].y << " " << kv.kvec_d[ik].z << std::endl;
    ofs << "\n  " << nlocal << " " << nlocal << std::endl;
    
    for(int i=0; i<nlocal; ++i)
    {
        for(int j=0; j<nlocal; ++j)
        {
            if(j%8==0) ofs << "\n";
            ofs << " " << dm_k[ik](i,j).real();
            //ofs << " " << DM[is][i][j];
        }
    }

    end = time(NULL);
    ModuleBase::GlobalFunc::OUT_TIME("write_dmk",start,end);
    ofs.close();
    
    ModuleBase::timer::tick("elecstate","write_dmk");

    return;
}

// read the density matrix in k-space
// weiqing add 2023/8/9
void elecstate::read_dmk(
    const K_Vectors& kv,
    const int& ik,
    const int& nlocal,
	const std::string &fn,
	std::vector<ModuleBase::ComplexMatrix> &dm_k)
{
    //
    bool quit_abacus = false;

    std::ifstream ifs;

    ifs.open(fn.c_str());
    if (!ifs)
    {
        quit_abacus = true;
    }
    else
    {
        // if the number is not match,
        // quit the program or not.
        bool quit=false;
            
        ModuleBase::CHECK_DOUBLE(ifs,kv.kvec_d[ik].x,quit);
        ModuleBase::CHECK_DOUBLE(ifs,kv.kvec_d[ik].y,quit);
        ModuleBase::CHECK_DOUBLE(ifs,kv.kvec_d[ik].z,quit);
        ModuleBase::CHECK_INT(ifs, nlocal);
        ModuleBase::CHECK_INT(ifs, nlocal);
    }// If file exist, read in data.
    // Finish reading the first part of density matrix.

    for(int i=0; i<nlocal; ++i)
    {
        for(int j=0; j<nlocal; ++j)
        {
            ifs >> dm_k[ik](i,j);
        }
    }

    ifs.close();

    return;
}
