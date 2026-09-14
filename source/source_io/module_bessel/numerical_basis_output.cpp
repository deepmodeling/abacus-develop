#include "numerical_basis.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/global_variable.h"

#include <cstring>
#include <fstream>
#include <vector>

void Numerical_Basis::output_info(std::ofstream& ofs, const Bessel_Basis& bessel_basis, const K_Vectors& kv,
                                  const UnitCell& ucell)
{
    // only print out to the information by the first processor
    if (GlobalV::MY_RANK == 0)
    {
        ofs.precision(10);
        ofs << ucell.lat0 << std::endl;

        ofs << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << std::endl;
        ofs << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << std::endl;
        ofs << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << std::endl;

        ofs << ucell.ntype << " ntype" << std::endl;
        for (int it = 0; it < ucell.ntype; it++)
        {
            ofs << ucell.atoms[it].label << " label" << std::endl; // mohan add 2009-07-23
            ofs << ucell.atoms[it].na << " na" << std::endl;
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                ofs << ucell.atoms[it].tau[ia].x << " " << ucell.atoms[it].tau[ia].y << " " << ucell.atoms[it].tau[ia].z
                    << std::endl;
            }
        }
        // ecutwfc_jlq determine the jlq corresponding to plane wave calculation.
        ofs << PARAM.inp.ecutwfc << " ecutwfc" << std::endl; // mohan add 2009-09-08

        // this parameter determine the total number of jlq.
        ofs << bessel_basis.get_ecut() << " ecutwfc_jlq" << std::endl; // mohan modify 2009-09-08
        ofs << bessel_basis.get_rcut() << " rcut_Jlq" << std::endl;

        // mohan add 'smooth' and 'sigma' 2009-08-28
        ofs << bessel_basis.get_smooth() << " smooth" << std::endl;
        ofs << bessel_basis.get_sigma() << " sigma" << std::endl;

        ofs << bessel_basis.get_tolerence() << " tolerence" << std::endl;

        ofs << ucell.lmax << " lmax" << std::endl;
    }

    ofs << std::scientific;

    ofs << std::setprecision(8);
    // NOTICE: ofs_warning << "\n The precison may affect the optimize result.";

    if (GlobalV::MY_RANK == 0)
    {
        ofs << kv.get_nkstot() << " nks" << std::endl;
        ofs << PARAM.inp.nbands << " nbands" << std::endl;
        ofs << PARAM.globalv.nlocal << " nwfc" << std::endl;
        ofs << bessel_basis.get_ecut_number() << " ne " << std::endl;
    }
}

void Numerical_Basis::output_k(std::ofstream& ofs, const K_Vectors& kv)
{
    // (1)
    if (GlobalV::MY_RANK == 0)
    {
        ofs << "<WEIGHT_OF_KPOINTS>";
    }

    // only half of nkstot should be output in "NSPIN == 2" case, k_up and k_down has same k infomation
    int nkstot = kv.get_nkstot();

    // (2)
    for (int ik = 0; ik < nkstot; ik++)
    {
        double kx, ky, kz, wknow;
#ifdef __MPI
        // temprary restrict kpar=1 for NSPIN=2 case for generating_orbitals
        int pool = 0;
        if (PARAM.inp.nspin != 2) {
            pool = kv.para_k.whichpool[ik];
        }
        const int iknow = ik - kv.para_k.startk_pool[GlobalV::MY_POOL];
        if (GlobalV::RANK_IN_POOL == 0)
        {
            if (GlobalV::MY_POOL == 0)
            {
                if (pool == 0)
                {
                    kx = kv.kvec_c[ik].x;
                    ky = kv.kvec_c[ik].y;
                    kz = kv.kvec_c[ik].z;
                    wknow = kv.wk[ik];
                }
                else
                {

                    int startpro_pool = kv.para_k.get_startpro_pool(pool);
                    MPI_Status ierror;
                    MPI_Recv(&kx, 1, MPI_DOUBLE, startpro_pool, ik * 4, MPI_COMM_WORLD, &ierror);
                    MPI_Recv(&ky, 1, MPI_DOUBLE, startpro_pool, ik * 4 + 1, MPI_COMM_WORLD, &ierror);
                    MPI_Recv(&kz, 1, MPI_DOUBLE, startpro_pool, ik * 4 + 2, MPI_COMM_WORLD, &ierror);
                    MPI_Recv(&wknow, 1, MPI_DOUBLE, startpro_pool, ik * 4 + 3, MPI_COMM_WORLD, &ierror);
                }
            }
            else
            {
                if (GlobalV::MY_POOL == pool)
                {
                    MPI_Send(&kv.kvec_c[iknow].x, 1, MPI_DOUBLE, 0, ik * 4, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].y, 1, MPI_DOUBLE, 0, ik * 4 + 1, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].z, 1, MPI_DOUBLE, 0, ik * 4 + 2, MPI_COMM_WORLD);
                    MPI_Send(&kv.wk[iknow], 1, MPI_DOUBLE, 0, ik * 4 + 3, MPI_COMM_WORLD);
                }
            }
        }
        // this barrier is very important
        MPI_Barrier(MPI_COMM_WORLD);
#else
        if (GlobalV::MY_RANK == 0)
        {
            kx = kv.kvec_c[ik].x;
            ky = kv.kvec_c[ik].y;
            kz = kv.kvec_c[ik].z;
            wknow = kv.wk[ik];
        }
#endif

        if (GlobalV::MY_RANK == 0)
        {
            ofs << "\n" << kx << " " << ky << " " << kz;
            ofs << " " << wknow * 0.5;
        }
    }

    if (GlobalV::MY_RANK == 0)
    {
        ofs << "\n</WEIGHT_OF_KPOINTS>" << std::endl;
    }
}

void Numerical_Basis::output_overlap_Q(std::ofstream& ofs, const std::vector<ModuleBase::ComplexArray>& overlap_Q,
                                       const K_Vectors& kv)
{
    // (3)
    if (GlobalV::MY_RANK == 0)
    {
        ofs << "\n<OVERLAP_Q>";
    }

    // Copy to overlap_Q_k for Pkpoints.pool_collection temporaly.
    // It's better to refactor to Pkpoints.pool_collection(overlap_Q) in the future.
    // Peize Lin comments 2021.07.25
    assert(kv.get_nks() > 0);
    ModuleBase::ComplexArray overlap_Q_k(kv.get_nks(), overlap_Q[0].getBound1(), overlap_Q[0].getBound2(),
                                         overlap_Q[0].getBound3());
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        std::memcpy(overlap_Q_k.ptr + ik * overlap_Q[ik].getSize(), overlap_Q[ik].ptr,
                    overlap_Q[ik].getSize() * sizeof(std::complex<double>));
    }

    // only half of nkstot should be output in "NSPIN == 2" case, k_up and k_down has same k infomation
    int nkstot = kv.get_nkstot();
    int count = 0;
    for (int ik = 0; ik < nkstot; ik++)
    {
        ModuleBase::ComplexArray Qtmp(overlap_Q[ik].getBound1(), overlap_Q[ik].getBound2(), overlap_Q[ik].getBound3());
        Qtmp.zero_out();
        kv.para_k.pool_collection(Qtmp.ptr, overlap_Q_k, ik);
        if (GlobalV::MY_RANK == 0)
        {
            // begin data writing.
            const int dim = Qtmp.getSize();
            for (int i = 0; i < dim; i++)
            {
                if (count % 4 == 0) {
                    ofs << std::endl;
                }
                ofs << " " << Qtmp.ptr[i].real() << " " << Qtmp.ptr[i].imag();
                ++count;
            }
            // end data writing.
        }
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    // (5)
    if (GlobalV::MY_RANK == 0)
    {
        ofs << "\n</OVERLAP_Q>" << std::endl;
    }
}

void Numerical_Basis::output_overlap_Sq(const std::string& name, std::ofstream& ofs,
                                        const std::vector<ModuleBase::ComplexArray>& overlap_Sq, const K_Vectors& kv)
{
    if (GlobalV::MY_RANK == 0)
    {
        ofs << "\n<OVERLAP_Sq>";
        ofs.close();
    }

    // only half of nkstot should be output in "NSPIN == 2" case, k_up and k_down has same k infomation
    int ispin = 1;
    if (PARAM.inp.nspin == 2) {
        ispin = 2;
    }
    int nkstot = kv.get_nkstot() / ispin;
    int count = 0;
    for (int is = 0; is < ispin; is++)
    {
        for (int ik = 0; ik < nkstot; ik++)
        {
            if (GlobalV::MY_POOL == kv.para_k.whichpool[ik])
            {
                if (GlobalV::RANK_IN_POOL == 0)
                {
                    ofs.open(name.c_str(), std::ios::app);
                    const int ik_now = ik - kv.para_k.startk_pool[GlobalV::MY_POOL] + is * nkstot;

                    const int size = overlap_Sq[ik_now].getSize();
                    for (int i = 0; i < size; i++)
                    {
                        if (count % 2 == 0) {
                            ofs << std::endl;
                        }
                        ofs << " " << overlap_Sq[ik_now].ptr[i].real() << " " << overlap_Sq[ik_now].ptr[i].imag();
                        ++count;
                    }
                    ofs.flush();
                    ofs.close();
                }
#ifdef __MPI
                MPI_Barrier(MPI_COMM_WORLD);
#endif
            }
            else
            {
#ifdef __MPI
                MPI_Barrier(MPI_COMM_WORLD);
#endif
            }
        }
    }
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(name.c_str(), std::ios::app);
        ofs << "\n</OVERLAP_Sq>" << std::endl;
    }
}

// Peize Lin add 2020.04.23
void Numerical_Basis::output_overlap_V(std::ofstream& ofs, const ModuleBase::matrix& overlap_V)
{
    if (GlobalV::MY_RANK == 0)
    {
        ofs << "\n<OVERLAP_V>" << std::endl;
        overlap_V.print(ofs);
        ofs << "</OVERLAP_V>" << std::endl;
    }
}
