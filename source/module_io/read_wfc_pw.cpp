#include "read_wfc_pw.h"

#include "binstream.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_common.h"
#include "module_base/parallel_global.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"

bool ModuleIO::read_wfc_pw(const std::string& filename,
                           const ModulePW::PW_Basis_K* pw_wfc,
                           const int& ik,
                           const int& nkstot,
                           ModuleBase::ComplexMatrix& wfc)
{
    ModuleBase::TITLE("ModuleIO", "read_wfc_pw");
    ModuleBase::timer::tick("ModuleIO", "read_wfc_pw");

    const int nx = pw_wfc->nx;
    const int ny = pw_wfc->ny;
    const int nz = pw_wfc->nz;
    const int npwk_max = pw_wfc->npwk_max;

    int npwtot, ikstot, max_dim;

    // get npwtot and ikstot
#ifdef __MPI
    MPI_Allreduce(&pw_wfc->npwk[ik], &npwtot, 1, MPI_INT, MPI_SUM, POOL_WORLD);
    MPI_Allreduce(&npwk_max, &max_dim, 1, MPI_INT, MPI_MAX, POOL_WORLD);
    int nkp = nkstot / GlobalV::KPAR;
    int rem = nkstot % GlobalV::KPAR;
    if (GlobalV::MY_POOL < rem)
    {
        ikstot = GlobalV::MY_POOL * nkp + GlobalV::MY_POOL + ik;
    }
    else
    {
        ikstot = GlobalV::MY_POOL * nkp + rem + ik;
    }
#else
    max_dim = npwk_max;
    npwtot = pw_wfc->npwk[ik];
    ikstot = ik;
#endif

    npwtot *= GlobalV::NPOL;

    Binstream rfs;
    std::ifstream ifs;
    bool error = false;
    int size = 0;
    std::string msg = "";
    std::string filetype = filename.substr(filename.length() - 3, 3);

    // whether can open the file
    if (GlobalV::RANK_IN_POOL == 0)
    {
        if (filetype == "txt")
        {
            ifs.open(filename);
            if (!ifs)
            {
                error = true;
                msg = "Can't open file " + filename;
            }
        }
        else if (filetype == "dat")
        {
            rfs.open(filename, "r");
            if (!rfs)
            {
                error = true;
                msg = "Can't open file " + filename;
            }
        }
        else
        {
            error = true;
            msg = "Unknown file type " + filetype;
        }
    }

#ifdef __MPI
    // bcast msg to root process
    int ip = 0;
    if (error)
    {
        ip = GlobalV::MY_RANK;
        size = msg.size();
    }
    MPI_Allreduce(MPI_IN_PLACE, &ip, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Bcast(&size, 1, MPI_INT, ip, MPI_COMM_WORLD);
    char* swap = new char[size + 1];
    if (error)
    {
        strcpy(swap, msg.c_str());
    }
    MPI_Bcast(swap, size + 1, MPI_CHAR, ip, MPI_COMM_WORLD);
    msg = static_cast<std::string>(swap);
    delete[] swap;

    MPI_Bcast(&error, 1, MPI_C_BOOL, ip, MPI_COMM_WORLD);
#endif

    if (error)
    {
        ModuleBase::WARNING("ModuleIO::read_wfc_pw", msg);
        return false;
    }

    // read in some information
    int ikstot_in, nkstot_in, npwtot_in, nbands_in;
    double kvec[3], weight, ecutwfc_in, lat0_in, tpiba_in;

    if (GlobalV::RANK_IN_POOL == 0)
    {
        if (filetype == "txt")
        {
        }
        else if (filetype == "dat")
        {
            rfs >> size >> ikstot_in >> nkstot_in >> kvec[0] >> kvec[1] >> kvec[2] >> weight >> npwtot_in >> nbands_in
                >> ecutwfc_in >> lat0_in >> tpiba_in >> size;
        }
    }

#ifdef __MPI
    MPI_Bcast(&ikstot_in, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(&nkstot_in, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(kvec, 3, MPI_DOUBLE, 0, POOL_WORLD);
    MPI_Bcast(&weight, 1, MPI_DOUBLE, 0, POOL_WORLD);
    MPI_Bcast(&npwtot_in, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(&nbands_in, 1, MPI_INT, 0, POOL_WORLD);
    MPI_Bcast(&ecutwfc_in, 1, MPI_DOUBLE, 0, POOL_WORLD);
    MPI_Bcast(&lat0_in, 1, MPI_DOUBLE, 0, POOL_WORLD);
    MPI_Bcast(&tpiba_in, 1, MPI_DOUBLE, 0, POOL_WORLD);
#endif

    if (ikstot_in != ikstot + 1 || nkstot_in != nkstot || npwtot_in != npwtot || nbands_in != GlobalV::NBANDS)
    {
        std::cout << "ikstot_in = " << ikstot_in << std::endl;
        std::cout << "ikstot = " << ikstot + 1 << std::endl;
        std::cout << "nkstot_in = " << nkstot_in << std::endl;
        std::cout << "nkstot = " << nkstot << std::endl;
        std::cout << "npwtot_in = " << npwtot_in << std::endl;
        std::cout << "npwtot = " << npwtot << std::endl;
        std::cout << "nbands_in = " << nbands_in << std::endl;
        std::cout << "nbands = " << GlobalV::NBANDS << std::endl;
        ModuleBase::WARNING_QUIT(
            "ModuleIO::read_wfc_pw",
            "ikstot_in != ikstot || nkstot_in != nkstot || npwtot_in != npwtot || nbands_in != GlobalV::NBANDS");
    }

    if (kvec[0] != pw_wfc->kvec_c[ik].x || kvec[1] != pw_wfc->kvec_c[ik].y || kvec[2] != pw_wfc->kvec_c[ik].z)
    {
        std::cout << "kvec_in[" << ikstot << "] = " << kvec[0] << " " << kvec[1] << " " << kvec[2] << std::endl;
        std::cout << "kvec[" << ikstot << "] = " << pw_wfc->kvec_c[ik].x << " " << pw_wfc->kvec_c[ik].y << " "
                  << pw_wfc->kvec_c[ik].z << std::endl;
        ModuleBase::WARNING_QUIT("ModuleIO::read_wfc_pw", "k vector in file is not the same as the one in memory");
    }

    if (lat0_in != pw_wfc->lat0 || tpiba_in != pw_wfc->tpiba)
    {
        std::cout << "lat0_in = " << lat0_in << std::endl;
        std::cout << "lat0 = " << pw_wfc->lat0 << std::endl;
        std::cout << "tpiba_in = " << tpiba_in << std::endl;
        std::cout << "tpiba = " << pw_wfc->tpiba << std::endl;
        ModuleBase::WARNING_QUIT("ModuleIO::read_wfc_pw", "lat0_in != pw_wfc->lat0 || tpiba_in != pw_wfc->tpiba");
    }

    // read in G
    ModuleBase::Vector3<double> G_in[3];
    if (GlobalV::RANK_IN_POOL == 0)
    {
        if (filetype == "txt")
        {
        }
        else if (filetype == "dat")
        {
            rfs >> size >> G_in[0].x >> G_in[0].y >> G_in[0].z >> G_in[1].x >> G_in[1].y >> G_in[1].z >> G_in[2].x
                >> G_in[2].y >> G_in[2].z >> size;
        }
    }

#ifdef __MPI
    MPI_Bcast(G_in, 3 * 3, MPI_DOUBLE, 0, POOL_WORLD);
#endif

    if (G_in[0].x != pw_wfc->G.e11 || G_in[0].y != pw_wfc->G.e12 || G_in[0].z != pw_wfc->G.e13
        || G_in[1].x != pw_wfc->G.e21 || G_in[1].y != pw_wfc->G.e22 || G_in[1].z != pw_wfc->G.e23
        || G_in[2].x != pw_wfc->G.e31 || G_in[2].y != pw_wfc->G.e32 || G_in[2].z != pw_wfc->G.e33)
    {
        std::cout << "G_in[0] = " << G_in[0].x << " " << G_in[0].y << " " << G_in[0].z << std::endl;
        std::cout << "G_in[1] = " << G_in[1].x << " " << G_in[1].y << " " << G_in[1].z << std::endl;
        std::cout << "G_in[2] = " << G_in[2].x << " " << G_in[2].y << " " << G_in[2].z << std::endl;
        std::cout << "G = " << pw_wfc->G.e11 << " " << pw_wfc->G.e12 << " " << pw_wfc->G.e13 << std::endl;
        std::cout << "G = " << pw_wfc->G.e21 << " " << pw_wfc->G.e22 << " " << pw_wfc->G.e23 << std::endl;
        std::cout << "G = " << pw_wfc->G.e31 << " " << pw_wfc->G.e32 << " " << pw_wfc->G.e33 << std::endl;
        ModuleBase::WARNING_QUIT("ModuleIO::read_wfc_pw", "G_in != G");
    }

    // read in miller index
    ModuleBase::Vector3<int>* miller = new ModuleBase::Vector3<int>[npwtot_in];
    int* glo_order = nullptr;
    if (GlobalV::RANK_IN_POOL == 0)
    {
        if (filetype == "txt")
        {
        }
        else if (filetype == "dat")
        {
            rfs >> size;
            for (int i = 0; i < npwtot_in; ++i)
            {
                rfs >> miller[i].x >> miller[i].y >> miller[i].z;
            }
            rfs >> size;
        }

        // map global index to read ordering for plane waves
        glo_order = new int[nx * ny * nz];
        for (int i = 0; i < nx * ny * nz; i++)
        {
            glo_order[i] = -1;
        }
        for (int i = 0; i < npwtot_in / GlobalV::NPOL; ++i)
        {
            int index = (miller[i].x * ny + miller[i].y) * nz + miller[i].z;
            glo_order[index] = i;
        }
    }

    // map local to global index for plane waves
    int* l2g_pw = new int[pw_wfc->npwk[ik]];
    for (int i = 0; i < pw_wfc->npwk[ik]; ++i)
    {
        int isz = pw_wfc->igl2isz_k[ik * npwk_max + i];
        int iz = isz % nz;
        int is = isz / nz;
        int ixy = pw_wfc->is2fftixy[is];
        int index = ixy * nz + iz;
        l2g_pw[i] = index;
    }

    // read in wfc
    std::complex<double>* wfc_in = new std::complex<double>[npwtot_in];
    for (int ib = 0; ib < nbands_in; ib++)
    {
        if (GlobalV::RANK_IN_POOL == 0)
        {
            if (filetype == "txt")
            {
            }
            else if (filetype == "dat")
            {
                rfs >> size;
                for (int i = 0; i < npwtot_in; ++i)
                {
                    rfs >> wfc_in[i];
                }
                rfs >> size;
            }
        }

        // distribute wave functions to processers
#ifdef __MPI
        for (int ip = 0; ip < GlobalV::NPROC_IN_POOL; ++ip)
        {
            if (ip != 0)
            {
                if (GlobalV::RANK_IN_POOL == ip)
                {
                    MPI_Send(l2g_pw, pw_wfc->npwk[ik], MPI_INT, 0, ip, POOL_WORLD);
                    MPI_Recv(&wfc(ib, 0),
                             pw_wfc->npwk[ik],
                             MPI_DOUBLE_COMPLEX,
                             0,
                             ip + GlobalV::NPROC_IN_POOL,
                             POOL_WORLD,
                             MPI_STATUS_IGNORE);
                    if (GlobalV::NPOL == 2)
                    {
                        MPI_Recv(&wfc(ib, npwk_max),
                                 pw_wfc->npwk[ik],
                                 MPI_DOUBLE_COMPLEX,
                                 0,
                                 ip + 2 * GlobalV::NPROC_IN_POOL,
                                 POOL_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                }
                if (GlobalV::RANK_IN_POOL == 0)
                {
                    int* ig_ip = new int[max_dim];
                    std::complex<double>* wfc_ip = new std::complex<double>[max_dim];

                    MPI_Status wfc_status;
                    MPI_Recv(ig_ip, max_dim, MPI_INT, ip, ip, POOL_WORLD, &wfc_status);
                    MPI_Get_count(&wfc_status, MPI_INT, &size);

                    for (int i = 0; i < size; i++)
                    {
                        wfc_ip[i] = wfc_in[glo_order[ig_ip[i]]];
                    }
                    MPI_Send(wfc_ip, size, MPI_DOUBLE_COMPLEX, ip, ip + GlobalV::NPROC_IN_POOL, POOL_WORLD);
                    if (GlobalV::NPOL == 2)
                    {
                        for (int i = 0; i < size; i++)
                        {
                            wfc_ip[i] = wfc_in[glo_order[ig_ip[i]] + npwtot_in / 2];
                        }
                        MPI_Send(wfc_ip, size, MPI_DOUBLE_COMPLEX, ip, ip + 2 * GlobalV::NPROC_IN_POOL, POOL_WORLD);
                    }
                    delete[] ig_ip;
                    delete[] wfc_ip;
                }
            }
            else
            {
                if (GlobalV::RANK_IN_POOL == 0)
                {
                    for (int i = 0; i < pw_wfc->npwk[ik]; ++i)
                    {
                        wfc(ib, i) = wfc_in[glo_order[l2g_pw[i]]];
                    }
                    if (GlobalV::NPOL == 2)
                    {
                        for (int i = 0; i < pw_wfc->npwk[ik]; ++i)
                        {
                            wfc(ib, i + npwk_max) = wfc_in[glo_order[l2g_pw[i]] + npwtot_in / 2];
                        }
                    }
                }
            }
            MPI_Barrier(POOL_WORLD);
        }
#else
        for (int i = 0; i < pw_wfc->npwk[ik]; ++i)
        {
            wfc(ib, i) = wfc_in[glo_order[l2g_pw[i]]];
        }
        if (GlobalV::NPOL == 2)
        {
            for (int i = 0; i < pw_wfc->npwk[ik]; ++i)
            {
                wfc(ib, i + npwk_max) = wfc_in[glo_order[l2g_pw[i]] + npwtot_in / 2];
            }
        }
#endif
    }

    delete[] l2g_pw;
    delete[] miller;
    delete[] wfc_in;

    if (GlobalV::RANK_IN_POOL == 0)
    {
        delete[] glo_order;
        ifs.close();
    }

    ModuleBase::timer::tick("ModuleIO", "read_wfc_pw");
    return true;
}