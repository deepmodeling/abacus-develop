#include "module_io/output_radial.h"
#include "module_base/global_function.h"

void ModuleIO::OutputRadial::initialize(const std::string& filename)
{
    file_to_.open(filename);
}

void ModuleIO::OutputRadial::write_header()
{
    if(file_to_.good())
    {
        std::vector<std::string> sublayers = {"S", "P", "D", "F", "G", "H", "I", "J", "K"};
        for(int i = 0; i < 75; ++i)
        {
            file_to_ << "-";
        }
        file_to_ << "\n";
        // left aligned
        file_to_ << std::left << std::setw(28) << "Element" << symbol_ << "\n";
        file_to_ << std::left << std::setw(28) << "Energy Cutoff(Ry)" << std::to_string(int(ecut_)) << "\n";
        // rcut .1f, not scientific
        file_to_ << std::left << std::setw(28) << "Radius Cutoff(a.u.)" << std::fixed << std::setprecision(1) << rcut_ << "\n";
        file_to_ << std::left << std::setw(28) << "Lmax" << lmax_ << "\n";
        for(int i = 0; i < lmax_ + 1; ++i)
        {
            std::string title = "Number of " + sublayers[i] + "orbital-->";
            file_to_ << std::left << std::setw(28) << title << l_nchi_[i] << "\n";
        }
        for(int i = 0; i < 75; ++i)
        {
            file_to_ << "-";
        }
        file_to_ << "\n";
        file_to_ << "SUMMARY  END\n\n";
        file_to_ << std::left << std::setw(28) << "Mesh" << std::setprecision(0) << ngrid_ << "\n";
        file_to_ << std::left << std::setw(28) << "dr" << std::setprecision(2) << dr_ << "\n";
    }
    else
    {
        ModuleBase::WARNING_QUIT("ModuleIO::OutputRadial::write_header", "file is not good");
    }
}

void ModuleIO::OutputRadial::configure(const std::string& symbol,
                                       const double& ecut,
                                       const double& rcut,
                                       const int& lmax, 
                                       const int* l_nchi,
                                       const int& ngrid,
                                       const double& dr
                                       )
{
    symbol_ = symbol;
    ecut_ = ecut;
    rcut_ = rcut;
    lmax_ = lmax;
    l_nchi_.resize(lmax_ + 1);
    for(int i = 0; i < lmax_ + 1; ++i)
    {
        l_nchi_[i] = l_nchi[i];
    }
    ngrid_ = ngrid;
    dr_ = dr;
    write_header();
}

void ModuleIO::OutputRadial::push(int ngrid, 
                                  const double* rgrid, 
                                  const double* chi)
{
    if(ngrid == 0)
    {
        ModuleBase::WARNING("ModuleIO::OutputRadial::push", "ngrid is 0, use the configured one");
        ngrid = ngrid_;
    }
    if(rgrid == nullptr)
    {
        ModuleBase::WARNING("ModuleIO::OutputRadial::push", "rgrid is nullptr, use the configured one");
    }

    file_to_ << std::right << std::setw(20) << "Type" 
             << std::right << std::setw(20) << "L" 
             << std::right << std::setw(20) << "N" << "\n";
    file_to_ << std::right << std::setw(20) 
             << std::to_string(0) 
             << std::right <<std::setw(20) 
             << std::to_string(current_l_) 
             << std::right <<std::setw(20) 
             << std::to_string(current_chi_);

    for(int i = 0; i < ngrid; ++i)
    {
        if(i % 4 == 0)
        {
            file_to_ << "\n";
        }
        file_to_ << std::left << std::setw(22) << std::setprecision(14) << std::scientific << chi[i];
    }
    file_to_ << "\n";

    ++current_chi_;
    if(current_chi_ == l_nchi_[current_l_])
    {
        current_chi_ = 0;
        ++current_l_;
    }
    if(current_l_ > lmax_ + 1)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::OutputRadial::push", "current_l_ is out of range");
    }
}

void ModuleIO::OutputRadial::finalize()
{
    file_to_.close();
}