#include "to_wannier90.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/math_sphbes.h"
#include "module_base/math_ylmreal.h"

toWannier90::toWannier90(int num_kpts, ModuleBase::Matrix3 recip_lattice)
{
    this->num_kpts = num_kpts;
    this->recip_lattice = recip_lattice;
    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
        this->cal_num_kpts = this->num_kpts;
    else if (GlobalV::NSPIN == 2)
        this->cal_num_kpts = this->num_kpts / 2;
}

toWannier90::toWannier90(int num_kpts, ModuleBase::Matrix3 recip_lattice, std::complex<double> ***wfc_k_grid_in)
{
    this->wfc_k_grid = wfc_k_grid_in;
    this->num_kpts = num_kpts;
    this->recip_lattice = recip_lattice;
    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
        this->cal_num_kpts = this->num_kpts;
    else if (GlobalV::NSPIN == 2)
        this->cal_num_kpts = this->num_kpts / 2;
}

toWannier90::~toWannier90()
{
    if (num_exclude_bands > 0)
        delete[] exclude_bands;
    if (GlobalV::BASIS_TYPE == "lcao")
        delete unk_inLcao;
    delete[] R_centre;
    delete[] L;
    delete[] m;
    delete[] rvalue;
    delete[] alfa;
    delete[] tag_cal_band;
}

void toWannier90::init_wannier_pw(const ModuleBase::matrix& ekb,
    const ModulePW::PW_Basis* rhopw,
    const ModulePW::PW_Basis_K* wfcpw,
    const ModulePW::PW_Basis_Big* bigpw,
    const K_Vectors& kv,
    const psi::Psi<std::complex<double>>* psi)
{
    this->read_nnkp(kv);

    if (GlobalV::NSPIN == 2)
    {
        wannier_spin = INPUT.wannier_spin;
        if (wannier_spin == "up")
            start_k_index = 0;
        else if (wannier_spin == "down")
            start_k_index = num_kpts / 2;
        else
        {
            ModuleBase::WARNING_QUIT("toWannier90::init_wannier", "Error wannier_spin set,is not \"up\" or \"down\" ");
        }
    }

    writeUNK(wfcpw, *psi, bigpw);
    outEIG(ekb);
    cal_Mmn(*psi, rhopw, wfcpw);
    cal_Amn(*psi, wfcpw);

    /*
    if(GlobalV::MY_RANK==0)
    {
        if(GlobalV::BASIS_TYPE == "pw")
        {
            cal_Amn(psi.evc);
            cal_Mmn(psi.evc);
            writeUNK(psi.evc);
            outEIG();
        }
        else if(GlobalV::BASIS_TYPE == "lcao")
        {
            getUnkFromLcao(wfcpw, kv, wfcpw->npwk_max);
            cal_Amn(this->unk_inLcao);
            cal_Mmn(this->unk_inLcao);
            writeUNK(this->unk_inLcao);
            outEIG();
        }
    }
    */
}

#ifdef __LCAO
void toWannier90::init_wannier_lcao(const Grid_Technique& gt,
                                    const ModuleBase::matrix& ekb,
                                    const ModulePW::PW_Basis* rhopw,
                                    const ModulePW::PW_Basis_K* wfcpw,
                                    const ModulePW::PW_Basis_Big* bigpw,
                                    const Structure_Factor& sf,
                                    const K_Vectors& kv,
                                    const psi::Psi<std::complex<double>>* psi)
{
    this->gridt = &gt;
    this->read_nnkp(kv);

    if (GlobalV::NSPIN == 2)
    {
        wannier_spin = INPUT.wannier_spin;
        if (wannier_spin == "up")
            start_k_index = 0;
        else if (wannier_spin == "down")
            start_k_index = num_kpts / 2;
        else
        {
            ModuleBase::WARNING_QUIT("toWannier90::init_wannier", "Error wannier_spin set,is not \"up\" or \"down\" ");
        }
    }

    getUnkFromLcao(wfcpw, sf, kv, wfcpw->npwk_max);
    cal_Amn(this->unk_inLcao[0], wfcpw);
    cal_Mmn(this->unk_inLcao[0], rhopw, wfcpw);
    writeUNK(wfcpw, this->unk_inLcao[0], bigpw);
    outEIG(ekb);
}
#endif

void toWannier90::read_nnkp(const K_Vectors& kv)
{
    // read *.nnkp file

    wannier_file_name = INPUT.nnkpfile;
    wannier_file_name = wannier_file_name.substr(0, wannier_file_name.length() - 5);

    GlobalV::ofs_running << "reading the " << wannier_file_name << ".nnkp file." << std::endl;

    std::ifstream nnkp_read(INPUT.nnkpfile.c_str(), ios::in);

    if (!nnkp_read)
        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error during readin parameters.");

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "real_lattice"))
    {
        ModuleBase::Matrix3 real_lattice_nnkp;
        nnkp_read >> real_lattice_nnkp.e11 >> real_lattice_nnkp.e12 >> real_lattice_nnkp.e13 >> real_lattice_nnkp.e21
            >> real_lattice_nnkp.e22 >> real_lattice_nnkp.e23 >> real_lattice_nnkp.e31 >> real_lattice_nnkp.e32
            >> real_lattice_nnkp.e33;

        real_lattice_nnkp = real_lattice_nnkp / GlobalC::ucell.lat0_angstrom;

        if (std::abs(real_lattice_nnkp.e11 - GlobalC::ucell.latvec.e11) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
        if (std::abs(real_lattice_nnkp.e12 - GlobalC::ucell.latvec.e12) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
        if (std::abs(real_lattice_nnkp.e13 - GlobalC::ucell.latvec.e13) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
        if (std::abs(real_lattice_nnkp.e21 - GlobalC::ucell.latvec.e21) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
        if (std::abs(real_lattice_nnkp.e22 - GlobalC::ucell.latvec.e22) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
        if (std::abs(real_lattice_nnkp.e23 - GlobalC::ucell.latvec.e23) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
        if (std::abs(real_lattice_nnkp.e31 - GlobalC::ucell.latvec.e31) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
        if (std::abs(real_lattice_nnkp.e32 - GlobalC::ucell.latvec.e32) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
        if (std::abs(real_lattice_nnkp.e33 - GlobalC::ucell.latvec.e33) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error real_lattice in *.nnkp file");
    }

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "recip_lattice"))
    {
        ModuleBase::Matrix3 recip_lattice_nnkp;
        nnkp_read >> recip_lattice_nnkp.e11 >> recip_lattice_nnkp.e12 >> recip_lattice_nnkp.e13
            >> recip_lattice_nnkp.e21 >> recip_lattice_nnkp.e22 >> recip_lattice_nnkp.e23 >> recip_lattice_nnkp.e31
            >> recip_lattice_nnkp.e32 >> recip_lattice_nnkp.e33;

        const double tpiba_angstrom = ModuleBase::TWO_PI / GlobalC::ucell.lat0_angstrom;
        recip_lattice_nnkp = recip_lattice_nnkp / tpiba_angstrom;

        if (std::abs(recip_lattice_nnkp.e11 - GlobalC::ucell.G.e11) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
        if (std::abs(recip_lattice_nnkp.e12 - GlobalC::ucell.G.e12) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
        if (std::abs(recip_lattice_nnkp.e13 - GlobalC::ucell.G.e13) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
        if (std::abs(recip_lattice_nnkp.e21 - GlobalC::ucell.G.e21) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
        if (std::abs(recip_lattice_nnkp.e22 - GlobalC::ucell.G.e22) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
        if (std::abs(recip_lattice_nnkp.e23 - GlobalC::ucell.G.e23) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
        if (std::abs(recip_lattice_nnkp.e31 - GlobalC::ucell.G.e31) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
        if (std::abs(recip_lattice_nnkp.e32 - GlobalC::ucell.G.e32) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
        if (std::abs(recip_lattice_nnkp.e33 - GlobalC::ucell.G.e33) > 1.0e-4)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error recip_lattice in *.nnkp file");
    }

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "kpoints"))
    {
        int numkpt_nnkp;
        ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, numkpt_nnkp);
        if ((GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4) && numkpt_nnkp != kv.nkstot)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error kpoints in *.nnkp file");
        else if (GlobalV::NSPIN == 2 && numkpt_nnkp != (kv.nkstot / 2))
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error kpoints in *.nnkp file");

        ModuleBase::Vector3<double> *kpoints_direct_nnkp = new ModuleBase::Vector3<double>[numkpt_nnkp];
        for (int ik = 0; ik < numkpt_nnkp; ik++)
        {
            nnkp_read >> kpoints_direct_nnkp[ik].x >> kpoints_direct_nnkp[ik].y >> kpoints_direct_nnkp[ik].z;
            if (std::abs(kpoints_direct_nnkp[ik].x - kv.kvec_d[ik].x) > 1.0e-4)
                ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error kpoints in *.nnkp file");
            if (std::abs(kpoints_direct_nnkp[ik].y - kv.kvec_d[ik].y) > 1.0e-4)
                ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error kpoints in *.nnkp file");
            if (std::abs(kpoints_direct_nnkp[ik].z - kv.kvec_d[ik].z) > 1.0e-4)
                ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "Error kpoints in *.nnkp file");
        }

        delete[] kpoints_direct_nnkp;

        ModuleBase::Vector3<double> my_gamma_point(0.0, 0.0, 0.0);
        // if( (kv.nkstot == 1) && (kv.kvec_d[0] == my_gamma_point) ) gamma_only_wannier = true;
    }

    if (GlobalV::NSPIN != 4)
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "projections"))
        {
            ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, num_wannier);
            // test
            // GlobalV::ofs_running << "num_wannier = " << num_wannier << std::endl;
            // test
            if (num_wannier < 0)
            {
                ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "wannier number is lower than 0");
            }

            R_centre = new ModuleBase::Vector3<double>[num_wannier];
            L = new int[num_wannier];
            m = new int[num_wannier];
            rvalue = new int[num_wannier];
            ModuleBase::Vector3<double> *z_axis = new ModuleBase::Vector3<double>[num_wannier];
            ModuleBase::Vector3<double> *x_axis = new ModuleBase::Vector3<double>[num_wannier];
            alfa = new double[num_wannier];

            for (int count = 0; count < num_wannier; count++)
            {
                nnkp_read >> R_centre[count].x >> R_centre[count].y >> R_centre[count].z;
                nnkp_read >> L[count] >> m[count];
                ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, rvalue[count]);
                nnkp_read >> z_axis[count].x >> z_axis[count].y >> z_axis[count].z;
                nnkp_read >> x_axis[count].x >> x_axis[count].y >> x_axis[count].z;
                ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, alfa[count]);
            }
            delete[] z_axis;
            delete[] x_axis;
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "noncolin spin is not done yet");
    }

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "nnkpts"))
    {
        ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, nntot);
        nnlist.resize(kv.nkstot);
        nncell.resize(kv.nkstot);
        for (int ik = 0; ik < kv.nkstot; ik++)
        {
            nnlist[ik].resize(nntot);
            nncell[ik].resize(nntot);
        }

        int numkpt_nnkp;
        if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
            numkpt_nnkp = kv.nkstot;
        else if (GlobalV::NSPIN == 2)
            numkpt_nnkp = kv.nkstot / 2;
        else
            throw std::runtime_error("numkpt_nnkp uninitialized in " + ModuleBase::GlobalFunc::TO_STRING(__FILE__)
                                     + " line " + ModuleBase::GlobalFunc::TO_STRING(__LINE__));

        for (int ik = 0; ik < numkpt_nnkp; ik++)
        {
            for (int ib = 0; ib < nntot; ib++)
            {
                int ik_nnkp;
                nnkp_read >> ik_nnkp;
                if (ik_nnkp != (ik + 1))
                    ModuleBase::WARNING_QUIT("toWannier90::read_nnkp", "error nnkpts in *.nnkp file");
                nnkp_read >> nnlist[ik][ib];
                nnkp_read >> nncell[ik][ib].x >> nncell[ik][ib].y >> nncell[ik][ib].z;
                nnlist[ik][ib]--; // this is c++ , begin from 0
            }
        }
    }

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(nnkp_read, "exclude_bands"))
    {
        ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, num_exclude_bands);
        if (num_exclude_bands > 0)
            exclude_bands = new int[num_exclude_bands];
        else if (num_exclude_bands < 0)
            ModuleBase::WARNING_QUIT("toWannier90::read_nnkp",
                                     "the exclude bands is wrong , please check *.nnkp file.");

        if (num_exclude_bands > 0)
        {
            for (int i = 0; i < num_exclude_bands; i++)
            {
                ModuleBase::GlobalFunc::READ_VALUE(nnkp_read, exclude_bands[i]);
                exclude_bands[i]--; // this is c++ , begin from 0
            }
        }
    }

    // test by jingan
    // GlobalV::ofs_running << "num_exclude_bands = " << num_exclude_bands << std::endl;
    // for(int i = 0; i < num_exclude_bands; i++)
    //{
    //	GlobalV::ofs_running << "exclude_bands : " << exclude_bands[i] << std::endl;
    //}
    // test by jingan

    nnkp_read.close();

    for (int i = 0; i < num_wannier; i++)
    {
        R_centre[i] = R_centre[i] * GlobalC::ucell.latvec;
        m[i] = m[i] - 1;
    }

    // test by jingan
    // GlobalV::ofs_running << "num_wannier is " << num_wannier << std::endl;
    // for(int i = 0; i < num_wannier; i++)
    //{
    //	GlobalV::ofs_running << "num_wannier" << std::endl;
    //	GlobalV::ofs_running << L[i] << " " << m[i] << " " << rvalue[i] << " " << alfa[i] << std::endl;
    //}
    // test by jingan

    tag_cal_band = new bool[GlobalV::NBANDS];
    if (GlobalV::NBANDS <= num_exclude_bands)
        ModuleBase::WARNING_QUIT("toWannier90::read_nnkp",
                                 "you set the band numer is not enough, please add bands number.");
    if (num_exclude_bands == 0)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            tag_cal_band[ib] = true;
    }
    else
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            tag_cal_band[ib] = true;
            for (int ibb = 0; ibb < num_exclude_bands; ibb++)
            {
                if (exclude_bands[ibb] == ib)
                {
                    tag_cal_band[ib] = false;
                    break;
                }
            }
        }
    }

    if (num_exclude_bands < 0)
        num_bands = GlobalV::NBANDS;
    else
        num_bands = GlobalV::NBANDS - num_exclude_bands;
}

void toWannier90::outEIG(const ModuleBase::matrix& ekb)
{
    if (GlobalV::MY_RANK == 0)
    {
        std::string fileaddress = GlobalV::global_out_dir + wannier_file_name + ".eig";
        std::ofstream eig_file(fileaddress.c_str());
        for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
        {
            int index_band = 0;
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                if (!tag_cal_band[ib])
                    continue;
                index_band++;
                eig_file << std::setw(5) << index_band << std::setw(5) << ik + 1 - start_k_index << std::setw(18)
                         << showpoint << fixed << std::setprecision(12)
                         << ekb(ik, ib) * ModuleBase::Ry_to_eV << std::endl;
            }
        }

        eig_file.close();
    }
}

void toWannier90::writeUNK(const ModulePW::PW_Basis_K* wfcpw,
                           const psi::Psi<std::complex<double>>& psi_pw,
                           const ModulePW::PW_Basis_Big* bigpw)
{

/*
    std::complex<double> *porter = new std::complex<double>[wfcpw->nrxx];

    for(int ik = start_k_index; ik < (cal_num_kpts+start_k_index); ik++)
    {
        std::stringstream name;
        if(GlobalV::NSPIN==1 || GlobalV::NSPIN==4)
        {
            name << GlobalV::global_out_dir << "UNK" << std::setw(5) << setfill('0') << ik+1 << ".1" ;
        }
        else if(GlobalV::NSPIN==2)
        {
            if(wannier_spin=="up") name << GlobalV::global_out_dir << "UNK" << std::setw(5) << setfill('0') <<
   ik+1-start_k_index << ".1" ; else if(wannier_spin=="down") name << GlobalV::global_out_dir << "UNK" << std::setw(5)
   << setfill('0') << ik+1-start_k_index << ".2" ;
        }

        std::ofstream unkfile(name.str());

        unkfile << std::setw(12) << rhopw->nx << std::setw(12) << rhopw->ny << std::setw(12) <<
   rhopw->nz << std::setw(12) << ik+1 << std::setw(12) << num_bands << std::endl;

        for(int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            if(!tag_cal_band[ib]) continue;
            //std::complex<double> *porter = GlobalC::UFFT.porter;
            //  u_k in real space
            ModuleBase::GlobalFunc::ZEROS(porter, GlobalC::rhopw->nrxx);
            for (int ig = 0; ig < psi_pw.get_ngk(ik); ig++)
            {
                porter[sf.ig2fftw[psi.igk(ik, ig)]] = psi_pw[ik](ib, ig);
            }
            sf.FFT_wfc.FFT3D(porter, 1);

            for(int k=0; k<rhopw->nz; k++)
            {
                for(int j=0; j<rhopw->ny; j++)
                {
                    for(int i=0; i<rhopw->nx; i++)
                    {
                        if(!gamma_only_wannier)
                        {
                            unkfile << std::setw(20) << std::setprecision(9) << std::setiosflags(ios::scientific) <<
   porter[i*rhopw->ny*rhopw->nz + j*rhopw->nz + k].real()
                                    << std::setw(20) << std::setprecision(9) << std::setiosflags(ios::scientific) <<
   porter[i*rhopw->ny*rhopw->nz + j*rhopw->nz + k].imag()
                                    //jingan test
                                    //<< "       " << std::setw(12) << std::setprecision(9) <<
   std::setiosflags(ios::scientific) << std::abs(porter[i*rhopw->ny*rhopw->nz + j*rhopw->nz + k])
                                    << std::endl;
                        }
                        else
                        {
                            double zero = 0.0;
                            unkfile << std::setw(20) << std::setprecision(9) << std::setiosflags(ios::scientific) <<
   std::abs( porter[i*rhopw->ny*rhopw->nz + j*rhopw->nz + k] )
                                    << std::setw(20) << std::setprecision(9) << std::setiosflags(ios::scientific) <<
   zero
                                    //jingan test
                                    //<< "       " << std::setw(12) << std::setprecision(9) <<
   std::setiosflags(ios::scientific) << std::abs(porter[i*rhopw->ny*rhopw->nz + j*rhopw->nz + k])
                                    << std::endl;
                        }
                    }
                }
            }


        }


        unkfile.close();

    }

    delete[] porter;
*/
#ifdef __MPI
    // num_z: how many planes on processor 'ip'
    int *num_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
    for (int iz = 0; iz < bigpw->nbz; iz++)
    {
        int ip = iz % GlobalV::NPROC_IN_POOL;
        num_z[ip] += bigpw->bz;
    }

    // start_z: start position of z in
    // processor ip.
    int *start_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
    for (int ip = 1; ip < GlobalV::NPROC_IN_POOL; ip++)
    {
        start_z[ip] = start_z[ip - 1] + num_z[ip - 1];
    }

    // which_ip: found iz belongs to which ip.
    int* which_ip = new int[wfcpw->nz];
    ModuleBase::GlobalFunc::ZEROS(which_ip, wfcpw->nz);
    for (int iz = 0; iz < wfcpw->nz; iz++)
    {
        for (int ip = 0; ip < GlobalV::NPROC_IN_POOL; ip++)
        {
            if (iz >= start_z[GlobalV::NPROC_IN_POOL - 1])
            {
                which_ip[iz] = GlobalV::NPROC_IN_POOL - 1;
                break;
            }
            else if (iz >= start_z[ip] && iz < start_z[ip + 1])
            {
                which_ip[iz] = ip;
                break;
            }
        }
    }

    // only do in the first pool.
    std::complex<double>* porter = new std::complex<double>[wfcpw->nrxx];
    int nxy = wfcpw->nx * wfcpw->ny;
    std::complex<double> *zpiece = new std::complex<double>[nxy];

    if (GlobalV::MY_POOL == 0)
    {
        for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
        {
            std::ofstream unkfile;

            if (GlobalV::MY_RANK == 0)
            {
                std::stringstream name;
                if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
                {
                    name << GlobalV::global_out_dir << "UNK" << std::setw(5) << setfill('0') << ik + 1 << ".1";
                }
                else if (GlobalV::NSPIN == 2)
                {
                    if (wannier_spin == "up")
                        name << GlobalV::global_out_dir << "UNK" << std::setw(5) << setfill('0')
                             << ik + 1 - start_k_index << ".1";
                    else if (wannier_spin == "down")
                        name << GlobalV::global_out_dir << "UNK" << std::setw(5) << setfill('0')
                             << ik + 1 - start_k_index << ".2";
                }

                unkfile.open(name.str(), ios::out);

                unkfile << std::setw(12) << wfcpw->nx << std::setw(12) << wfcpw->ny << std::setw(12) << wfcpw->nz
                        << std::setw(12) << ik + 1 << std::setw(12) << num_bands << std::endl;
            }

            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                if (!tag_cal_band[ib])
                    continue;

                wfcpw->recip2real(&psi_pw(ik, ib, 0), porter, ik);

                // save the rho one z by one z.
                for (int iz = 0; iz < wfcpw->nz; iz++)
                {
                    // tag must be different for different iz.
                    ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
                    int tag = iz;
                    MPI_Status ierror;

                    // case 1: the first part of rho in processor 0.
                    if (which_ip[iz] == 0 && GlobalV::RANK_IN_POOL == 0)
                    {
                        for (int ir = 0; ir < nxy; ir++)
                        {
                            zpiece[ir] = porter[ir * wfcpw->nplane + iz - wfcpw->startz_current];
                        }
                    }
                    // case 2: > first part rho: send the rho to
                    // processor 0.
                    else if (which_ip[iz] == GlobalV::RANK_IN_POOL)
                    {
                        for (int ir = 0; ir < nxy; ir++)
                        {
                            zpiece[ir] = porter[ir * wfcpw->nplane + iz - wfcpw->startz_current];
                        }
                        MPI_Send(zpiece, nxy, MPI_DOUBLE_COMPLEX, 0, tag, POOL_WORLD);
                    }

                    // case 2: > first part rho: processor 0 receive the rho
                    // from other processors
                    else if (GlobalV::RANK_IN_POOL == 0)
                    {
                        MPI_Recv(zpiece, nxy, MPI_DOUBLE_COMPLEX, which_ip[iz], tag, POOL_WORLD, &ierror);
                    }

                    // write data
                    if (GlobalV::MY_RANK == 0)
                    {
                        for (int iy = 0; iy < wfcpw->ny; iy++)
                        {
                            for (int ix = 0; ix < wfcpw->nx; ix++)
                            {
                                unkfile << std::setw(20) << std::setprecision(9) << std::setiosflags(ios::scientific)
                                        << zpiece[ix * wfcpw->ny + iy].real() << std::setw(20) << std::setprecision(9)
                                        << std::setiosflags(ios::scientific) << zpiece[ix * wfcpw->ny + iy].imag()
                                        << std::endl;
                            }
                        }
                    }
                } // end iz
                MPI_Barrier(POOL_WORLD);
            }

            if (GlobalV::MY_RANK == 0)
            {
                unkfile.close();
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    delete[] num_z;
    delete[] start_z;
    delete[] which_ip;
    delete[] porter;
    delete[] zpiece;

#endif
}

void toWannier90::cal_Amn(const psi::Psi<std::complex<double>>& psi_pw, const ModulePW::PW_Basis_K* wfcpw)
{
    const int pwNumberMax = wfcpw->npwk_max;

    std::ofstream Amn_file;

    if (GlobalV::MY_RANK == 0)
    {
        time_t time_now = time(NULL);
        std::string fileaddress = GlobalV::global_out_dir + wannier_file_name + ".amn";
        Amn_file.open(fileaddress.c_str(), ios::out);
        Amn_file << " Created on " << ctime(&time_now);
        Amn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << num_wannier
                 << std::endl;
    }

    ModuleBase::ComplexMatrix *trial_orbitals = new ModuleBase::ComplexMatrix[cal_num_kpts];
    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        trial_orbitals[ik].create(num_wannier, pwNumberMax);
        produce_trial_in_pw(psi_pw, ik, wfcpw, trial_orbitals[ik]);
    }

    // test by jingan
    // GlobalV::ofs_running << __FILE__ << __LINE__ << "start_k_index = " << start_k_index << "  cal_num_kpts = " <<
    // cal_num_kpts << std::endl;
    // test by jingan

    for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
    {
        for (int iw = 0; iw < num_wannier; iw++)
        {
            int index_band = 0;
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                if (!tag_cal_band[ib])
                    continue;
                index_band++;
                std::complex<double> amn(0.0, 0.0);
                std::complex<double> amn_tem(0.0, 0.0);
                for (int ig = 0; ig < pwNumberMax; ig++)
                {
                    int cal_ik = ik - start_k_index;
                    amn_tem = amn_tem + conj(psi_pw(ik, ib, ig)) * trial_orbitals[cal_ik](iw, ig);
                }
#ifdef __MPI
                MPI_Allreduce(&amn_tem, &amn, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#else
                amn = amn_tem;
#endif
                if (GlobalV::MY_RANK == 0)
                {
                    Amn_file << std::setw(5) << index_band << std::setw(5) << iw + 1 << std::setw(5)
                             << ik + 1 - start_k_index << std::setw(18) << showpoint << fixed << std::setprecision(12)
                             << amn.real() << std::setw(18) << showpoint << fixed << std::setprecision(12)
                             << amn.imag()
                             // jingan test
                             //<< "   " << std::setw(18) << std::setprecision(13) << std::abs(amn)
                             << std::endl;
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0)
        Amn_file.close();

    delete[] trial_orbitals;
}

void toWannier90::cal_Mmn(const psi::Psi<std::complex<double>>& psi_pw,
                          const ModulePW::PW_Basis* rhopw,
                          const ModulePW::PW_Basis_K* wfcpw)
{
    // test by jingan
    // GlobalV::ofs_running << __FILE__ << __LINE__ << " cal_num_kpts = " << cal_num_kpts << std::endl;
    // test by jingan

    std::ofstream mmn_file;

    if (GlobalV::MY_RANK == 0)
    {
        std::string fileaddress = GlobalV::global_out_dir + wannier_file_name + ".mmn";
        mmn_file.open(fileaddress.c_str(), ios::out);

        time_t time_now = time(NULL);
        mmn_file << " Created on " << ctime(&time_now);
        mmn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << nntot << std::endl;
    }

    /*
    ModuleBase::ComplexMatrix Mmn(GlobalV::NBANDS,GlobalV::NBANDS);
    if(gamma_only_wannier)
    {
        for(int ib = 0; ib < nntot; ib++)
        {
            ModuleBase::Vector3<double> phase_G = nncell[0][ib];
            for(int m = 0; m < GlobalV::NBANDS; m++)
            {
                if(!tag_cal_band[m]) continue;
                for(int n = 0; n <= m; n++)
                {
                    if(!tag_cal_band[n]) continue;
                    std::complex<double> mmn_tem = gamma_only_cal(m,n,psi_pw,phase_G);
                    Mmn(m,n) = mmn_tem;
                    if(m!=n) Mmn(n,m) = Mmn(m,n);
                }
            }
        }
    }
    */

    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        for (int ib = 0; ib < nntot; ib++)
        {
            int ikb = nnlist[ik][ib];

            ModuleBase::Vector3<double> phase_G = nncell[ik][ib];

            if (GlobalV::MY_RANK == 0)
            {
                mmn_file << std::setw(5) << ik + 1 << std::setw(5) << ikb + 1 << std::setw(5) << int(phase_G.x)
                         << std::setw(5) << int(phase_G.y) << std::setw(5) << int(phase_G.z) << std::endl;
            }

            for (int m = 0; m < GlobalV::NBANDS; m++)
            {
                if (!tag_cal_band[m])
                    continue;
                for (int n = 0; n < GlobalV::NBANDS; n++)
                {
                    if (!tag_cal_band[n])
                        continue;
                    std::complex<double> mmn(0.0, 0.0);

                    if (!gamma_only_wannier)
                    {
                        int cal_ik = ik + start_k_index;
                        int cal_ikb = ikb + start_k_index;
                        // test by jingan
                        // GlobalV::ofs_running << __FILE__ << __LINE__ << "cal_ik = " << cal_ik << "cal_ikb = " <<
                        // cal_ikb << std::endl;
                        // test by jingan
                        // std::complex<double> *unk_L_r = new std::complex<double>[wfcpw->nrxx];
                        // ToRealSpace(cal_ik,n,psi_pw,unk_L_r,phase_G);
                        // mmn = unkdotb(unk_L_r,cal_ikb,m,psi_pw);
                        mmn = unkdotkb(rhopw, wfcpw, cal_ik, cal_ikb, n, m, phase_G, psi_pw);
                        // delete[] unk_L_r;
                    }
                    else
                    {
                        // GlobalV::ofs_running << "gamma only test" << std::endl;
                        // mmn = Mmn(n,m);
                    }

                    if (GlobalV::MY_RANK == 0)
                    {
                        mmn_file << std::setw(18) << std::setprecision(12) << showpoint << fixed << mmn.real()
                                 << std::setw(18) << std::setprecision(12) << showpoint << fixed
                                 << mmn.imag()
                                 // jingan test
                                 //<< "    " << std::setw(12) << std::setprecision(9) << std::abs(mmn)
                                 << std::endl;
                    }
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0)
        mmn_file.close();
}

void toWannier90::produce_trial_in_pw(const psi::Psi<std::complex<double>>& psi_pw,
                                      const int& ik,
                                      const ModulePW::PW_Basis_K* wfcpw,
                                      ModuleBase::ComplexMatrix& trial_orbitals_k)
{

    for (int i = 0; i < num_wannier; i++)
    {
        if (L[i] < -5 || L[i] > 3)
            std::cout << "toWannier90::produce_trial_in_pw() your L angular momentum is wrong , please check !!! "
                      << std::endl;

        if (L[i] >= 0)
        {
            if (m[i] < 0 || m[i] > 2 * L[i])
                std::cout << "toWannier90::produce_trial_in_pw() your m momentum is wrong , please check !!! "
                          << std::endl;
        }
        else
        {
            if (m[i] < 0 || m[i] > -L[i])
                std::cout << "toWannier90::produce_trial_in_pw() your m momentum is wrong , please check !!! "
                          << std::endl;
        }
    }

    const int npw = wfcpw->npwk[ik];
    const int npwx = wfcpw->npwk_max;
    const int total_lm = 16;
    ModuleBase::matrix ylm(total_lm, npw);

    double bs2, bs3, bs6, bs12;
    bs2 = 1.0 / sqrt(2.0);
    bs3 = 1.0 / sqrt(3.0);
    bs6 = 1.0 / sqrt(6.0);
    bs12 = 1.0 / sqrt(12.0);

    ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = wfcpw->getgpluskcar(ik,ig);
    }

    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk, ylm);

    // test by jingan
    // GlobalV::ofs_running << "the mathzone::ylm_real is successful!" << std::endl;
    // GlobalV::ofs_running << "produce_trial_in_pw: num_wannier is " << num_wannier << std::endl;
    // test by jingan

    const int mesh_r = 333;
    const double dx = 0.025;
    const double x_min = -6.0;
    ModuleBase::matrix r(num_wannier, mesh_r);
    ModuleBase::matrix dr(num_wannier, mesh_r);
    ModuleBase::matrix psi(num_wannier, mesh_r);
    ModuleBase::matrix psir(num_wannier, mesh_r);
    ModuleBase::matrix psik(num_wannier, npw);

    for (int i = 0; i < num_wannier; i++)
    {
        double x = 0;
        for (int ir = 0; ir < mesh_r; ir++)
        {
            x = x_min + ir * dx;
            r(i, ir) = exp(x) / alfa[i];
            dr(i, ir) = dx * r(i, ir);
        }
    }

    for (int i = 0; i < num_wannier; i++)
    {
        double alfa32 = pow(alfa[i], 3.0 / 2.0);
        double alfa_new = alfa[i];
        int wannier_index = i;

        if (rvalue[i] == 1)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi(wannier_index, ir) = 2.0 * alfa32 * exp(-alfa_new * r(wannier_index, ir));
            }
        }

        if (rvalue[i] == 2)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi(wannier_index, ir) = 1.0 / sqrt(8.0) * alfa32 * (2.0 - alfa_new * r(wannier_index, ir))
                                         * exp(-alfa_new * r(wannier_index, ir) * 0.5);
            }
        }

        if (rvalue[i] == 3)
        {
            for (int ir = 0; ir < mesh_r; ir++)
            {
                psi(wannier_index, ir)
                    = sqrt(4.0 / 27.0) * alfa32
                      * (1.0 - 2.0 / 3.0 * alfa_new * r(wannier_index, ir)
                         + 2.0 / 27.0 * pow(alfa_new, 2.0) * r(wannier_index, ir) * r(wannier_index, ir))
                      * exp(-alfa_new * r(wannier_index, ir) * 1.0 / 3.0);
            }
        }
    }

    for (int i = 0; i < num_wannier; i++)
    {
        for (int ir = 0; ir < mesh_r; ir++)
        {
            psir(i, ir) = psi(i, ir) * r(i, ir);
        }
    }

    for (int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
    {
        if (L[wannier_index] >= 0)
        {
            get_trial_orbitals_lm_k(wannier_index,
                                    L[wannier_index],
                                    m[wannier_index],
                                    ylm,
                                    dr,
                                    r,
                                    psir,
                                    mesh_r,
                                    gk,
                                    npw,
                                    npwx,
                                    trial_orbitals_k);
        }
        else
        {
            if (L[wannier_index] == -1 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs2 * tem_array[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array;
            }
            else if (L[wannier_index] == -1 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs2 * tem_array[ig] - bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array;
            }
            else if (L[wannier_index] == -2 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -2 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] - bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -2 && m[wannier_index] == 2)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array[ig] + 2 * bs6 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array;
            }
            else if (L[wannier_index] == -3 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = 0.5
                          * (tem_array_1[ig] + tem_array_2[ig] + tem_array_3[ig] + trial_orbitals_k(wannier_index, ig));
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -3 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = 0.5
                          * (tem_array_1[ig] + tem_array_2[ig] - tem_array_3[ig] - trial_orbitals_k(wannier_index, ig));
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -3 && m[wannier_index] == 2)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = 0.5
                          * (tem_array_1[ig] - tem_array_2[ig] + tem_array_3[ig] - trial_orbitals_k(wannier_index, ig));
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -3 && m[wannier_index] == 3)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = 0.5
                          * (tem_array_1[ig] - tem_array_2[ig] - tem_array_3[ig] + trial_orbitals_k(wannier_index, ig));
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - bs6 * tem_array_2[ig] - bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 2)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs3 * tem_array_1[ig] - 2 * bs6 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 3)
            {
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs2 * tem_array_1[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
            }
            else if (L[wannier_index] == -4 && m[wannier_index] == 4)
            {
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = -1.0 * bs2 * tem_array_1[ig] + bs2 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 0)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig) = bs6 * tem_array_1[ig] - bs2 * tem_array_2[ig]
                                                          - bs12 * tem_array_3[ig]
                                                          + 0.5 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 1)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 1, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig) = bs6 * tem_array_1[ig] + bs2 * tem_array_2[ig]
                                                          - bs12 * tem_array_3[ig]
                                                          + 0.5 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 2)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig) = bs6 * tem_array_1[ig] - bs2 * tem_array_2[ig]
                                                          - bs12 * tem_array_3[ig]
                                                          - 0.5 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 3)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 2, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_3 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_3[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 3, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig) = bs6 * tem_array_1[ig] + bs2 * tem_array_2[ig]
                                                          - bs12 * tem_array_3[ig]
                                                          - 0.5 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
                delete[] tem_array_3;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 4)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs6 * tem_array_1[ig] - bs2 * tem_array_2[ig] + bs3 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
            else if (L[wannier_index] == -5 && m[wannier_index] == 5)
            {
                get_trial_orbitals_lm_k(wannier_index, 0, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_1 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_1[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 1, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                std::complex<double> *tem_array_2 = new std::complex<double>[npwx];
                for (int ig = 0; ig < npwx; ig++)
                {
                    tem_array_2[ig] = trial_orbitals_k(wannier_index, ig);
                }
                get_trial_orbitals_lm_k(wannier_index, 2, 0, ylm, dr, r, psir, mesh_r, gk, npw, npwx, trial_orbitals_k);
                for (int ig = 0; ig < npwx; ig++)
                {
                    trial_orbitals_k(wannier_index, ig)
                        = bs6 * tem_array_1[ig] + bs2 * tem_array_2[ig] + bs3 * trial_orbitals_k(wannier_index, ig);
                }
                delete[] tem_array_1;
                delete[] tem_array_2;
            }
        }
    }
    delete[] gk;
}

void toWannier90::get_trial_orbitals_lm_k(const int wannier_index,
                                          const int orbital_L,
                                          const int orbital_m,
                                          ModuleBase::matrix &ylm,
                                          ModuleBase::matrix &dr,
                                          ModuleBase::matrix &r,
                                          ModuleBase::matrix &psir,
                                          const int mesh_r,
                                          ModuleBase::Vector3<double> *gk,
                                          const int npw,
                                          const int npwx,
                                          ModuleBase::ComplexMatrix &trial_orbitals_k)
{

    double *psik = new double[npw];
    double *psir_tem = new double[mesh_r];
    double *r_tem = new double[mesh_r];
    double *dr_tem = new double[mesh_r];
    double *psik_tem = new double[GlobalV::NQX];
    ModuleBase::GlobalFunc::ZEROS(psir_tem, mesh_r);
    ModuleBase::GlobalFunc::ZEROS(r_tem, mesh_r);
    ModuleBase::GlobalFunc::ZEROS(dr_tem, mesh_r);

    for (int ir = 0; ir < mesh_r; ir++)
    {
        psir_tem[ir] = psir(wannier_index, ir);
        r_tem[ir] = r(wannier_index, ir);
        dr_tem[ir] = dr(wannier_index, ir);
    }

    toWannier90::integral(mesh_r, psir_tem, r_tem, dr_tem, orbital_L, psik_tem);

    for (int ig = 0; ig < npw; ig++)
    {
        psik[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(psik_tem,
                                                                 GlobalV::NQX,
                                                                 GlobalV::DQ,
                                                                 gk[ig].norm() * GlobalC::ucell.tpiba);
    }

    std::complex<double> *sk = new std::complex<double>[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        const double arg = (gk[ig] * R_centre[wannier_index]) * ModuleBase::TWO_PI;
        sk[ig] = std::complex<double>(cos(arg), -sin(arg));
    }

    double *wannier_ylm = new double[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        int index = orbital_L * orbital_L + orbital_m;
        if (index == 2 || index == 3 || index == 5 || index == 6 || index == 14 || index == 15)
        {
            wannier_ylm[ig] = -1 * ylm(index, ig);
        }
        else
        {
            wannier_ylm[ig] = ylm(index, ig);
        }
    }

    std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, orbital_L);
    for (int ig = 0; ig < npwx; ig++)
    {
        if (ig < npw)
        {
            trial_orbitals_k(wannier_index, ig) = lphase * sk[ig] * wannier_ylm[ig] * psik[ig];
        }
        else
            trial_orbitals_k(wannier_index, ig) = std::complex<double>(0.0, 0.0);
    }

    std::complex<double> anorm(0.0, 0.0);
    for (int ig = 0; ig < npwx; ig++)
    {
        anorm = anorm + conj(trial_orbitals_k(wannier_index, ig)) * trial_orbitals_k(wannier_index, ig);
    }

    std::complex<double> anorm_tem(0.0, 0.0);
#ifdef __MPI
    MPI_Allreduce(&anorm, &anorm_tem, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#else
    anorm_tem = anorm;
#endif

    for (int ig = 0; ig < npwx; ig++)
    {
        trial_orbitals_k(wannier_index, ig) = trial_orbitals_k(wannier_index, ig) / sqrt(anorm_tem);
    }

    delete[] psik;
    delete[] psir_tem;
    delete[] r_tem;
    delete[] dr_tem;
    delete[] psik_tem;
    delete[] sk;
    delete[] wannier_ylm;

    return;
}

void toWannier90::integral(const int meshr,
                           const double *psir,
                           const double *r,
                           const double *rab,
                           const int &l,
                           double *table)
{
    const double pref = ModuleBase::FOUR_PI / sqrt(GlobalC::ucell.omega);

    double *inner_part = new double[meshr];
    for (int ir = 0; ir < meshr; ir++)
    {
        inner_part[ir] = psir[ir] * psir[ir];
    }

    double unit = 0.0;
    ModuleBase::Integral::Simpson_Integral(meshr, inner_part, rab, unit);
    delete[] inner_part;

    double *aux = new double[meshr];
    double *vchi = new double[meshr];
    for (int iq = 0; iq < GlobalV::NQX; iq++)
    {
        const double q = GlobalV::DQ * iq;
        ModuleBase::Sphbes::Spherical_Bessel(meshr, r, q, l, aux);
        for (int ir = 0; ir < meshr; ir++)
        {
            vchi[ir] = psir[ir] * aux[ir] * r[ir];
        }

        double vqint = 0.0;
        ModuleBase::Integral::Simpson_Integral(meshr, vchi, rab, vqint);

        table[iq] = vqint * pref;
    }
    delete[] aux;
    delete[] vchi;
    return;
}

/*
void toWannier90::ToRealSpace(const int &ik,
                              const int &ib,
                              const ModuleBase::ComplexMatrix *evc,
                              std::complex<double> *psir,
                              const ModuleBase::Vector3<double> G)
{
    // (1) set value
    std::complex<double> *phase = GlobalC::UFFT.porter;
    ModuleBase::GlobalFunc::ZEROS(psir, wfcpw->nrxx);
    ModuleBase::GlobalFunc::ZEROS(phase, wfcpw->nrxx);

    for (int ig = 0; ig < psi->get_ngk(ik); ig++)
    {
        psir[wfcpw->ng2fftw[psi.igk(ik, ig)]] = evc[ik](ib, ig);
    }

    // get the phase value in realspace
    for (int ig = 0; ig < wfcpw->ngmw; ig++)
    {
        if (wfcpw->ndirect[ig] == G)
        {
            phase[wfcpw->ng2fftw[ig]] = std::complex<double>(1.0, 0.0);
            break;
        }
    }
    // (2) fft and get value
    wfcpw->nFT_wfc.FFT3D(psir, 1);
    wfcpw->nFT_wfc.FFT3D(phase, 1);

    for (int ir = 0; ir < wfcpw->nrxx; ir++)
    {
        psir[ir] = psir[ir] * phase[ir];
    }
    return;
}

std::complex<double> toWannier90::unkdotb(const std::complex<double> *psir,
                                          const int ikb,
                                          const int bandindex,
                                          const ModuleBase::ComplexMatrix *psi_pw)
{
    std::complex<double> result(0.0, 0.0);
    int knumber = psi->get_ngk(ikb);
    std::complex<double> *porter = GlobalC::UFFT.porter;
    ModuleBase::GlobalFunc::ZEROS(porter, wfcpw->nrxx);
    for (int ir = 0; ir < wfcpw->nrxx; ir++)
    {
        porter[ir] = psir[ir];
    }
    wfcpw->nFT_wfc.FFT3D(porter, -1);

    for (int ig = 0; ig < knumber; ig++)
    {
        result = result + conj(porter[wfcpw->ng2fftw[psi.igk(ikb, ig)]]) * psi_pw[ikb](bandindex, ig);
    }
    return result;
}
*/
std::complex<double> toWannier90::unkdotkb(const ModulePW::PW_Basis* rhopw,
                                           const ModulePW::PW_Basis_K* wfcpw,
                                           const int& ik,
                                           const int& ikb,
                                           const int& iband_L,
                                           const int& iband_R,
                                           const ModuleBase::Vector3<double> G,
                                           const psi::Psi<std::complex<double>>& psi_pw)
{
    // (1) set value
    std::complex<double> result(0.0, 0.0);
    std::complex<double> *psir = new std::complex<double>[wfcpw->nmaxgr];
    std::complex<double>* phase = new std::complex<double>[rhopw->nmaxgr];

    // get the phase value in realspace
    for (int ig = 0; ig < rhopw->npw; ig++)
    {
        if (rhopw->gdirect[ig] == G) // It should be used carefully. We cannot judge if two double are equal.
        {
            phase[ig] = std::complex<double>(1.0, 0.0);
            break;
        }
    }

    // (2) fft and get value
    rhopw->recip2real(phase, phase);
    wfcpw->recip2real(&psi_pw(ik, iband_L, 0), psir, ik);

    for (int ir = 0; ir < wfcpw->nrxx; ir++)
    {
        psir[ir] *= phase[ir];
    }

    wfcpw->real2recip(psir, psir, ik);

    std::complex<double> result_tem(0.0, 0.0);

    for (int ig = 0; ig < psi_pw.get_ngk(ikb); ig++)
    {
        result_tem = result_tem + conj(psir[ig]) * psi_pw(ikb, iband_R, ig);
    }
#ifdef __MPI
    MPI_Allreduce(&result_tem, &result, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#else
    result = result_tem;
#endif
    delete[] psir;
    delete[] phase;
    return result;
}

/*
std::complex<double> toWannier90::gamma_only_cal(const int &ib_L,
                                                   const int &ib_R,
                                                   const ModuleBase::ComplexMatrix *psi_pw,
                                                   const ModuleBase::Vector3<double> G)
{
    std::complex<double> *phase = new std::complex<double>[wfcpw->nrxx];
    std::complex<double> *psir = new std::complex<double>[wfcpw->nrxx];
    std::complex<double> *psir_2 = new std::complex<double>[wfcpw->nrxx];
    ModuleBase::GlobalFunc::ZEROS(phase, wfcpw->nrxx);
    ModuleBase::GlobalFunc::ZEROS(psir, wfcpw->nrxx);
    ModuleBase::GlobalFunc::ZEROS(psir_2, wfcpw->nrxx);

    for (int ig = 0; ig < psi->get_ngk(0); ig++)
    {
        // psir[ wfcpw->ng2fftw[ psi.igk(0,ig) ] ] = psi_pw[0](ib_L, ig);
        psir[wfcpw->ng2fftw[psi.igk(0, ig)]] = std::complex<double>(std::abs(psi_pw[0](ib_L, ig)), 0.0);
    }

    // get the phase value in realspace
    for (int ig = 0; ig < wfcpw->ngmw; ig++)
    {
        if (wfcpw->ndirect[ig] == G)
        {
            phase[wfcpw->ng2fftw[ig]] = std::complex<double>(1.0, 0.0);
            break;
        }
    }
    // (2) fft and get value
    wfcpw->nFT_wfc.FFT3D(psir, 1);
    wfcpw->nFT_wfc.FFT3D(phase, 1);

    for (int ir = 0; ir < wfcpw->nrxx; ir++)
    {
        psir_2[ir] = conj(psir[ir]) * phase[ir];
    }

    for (int ir = 0; ir < wfcpw->nrxx; ir++)
    {
        psir[ir] = psir[ir] * phase[ir];
    }

    wfcpw->nFT_wfc.FFT3D(psir, -1);
    wfcpw->nFT_wfc.FFT3D(psir_2, -1);

    std::complex<double> result(0.0, 0.0);

    for (int ig = 0; ig < psi->get_ngk(0); ig++)
    {
        // result = result + conj(psir_2[ wfcpw->ng2fftw[psi.igk(0,ig)] ]) * psi_pw[0](ib_R,ig) + psir[
wfcpw->ng2fftw[ psi.igk(0,ig)] ] * conj(psi_pw[0](ib_R,ig));
// std::complex<double> tem = std::complex<double>( std::abs(psi_pw[0](ib_R,ig)), 0.0 );
result = result + conj(psir[wfcpw->ng2fftw[psi.igk(0, ig)]]); // * tem;
    }

    delete[] phase;
    delete[] psir;
    delete[] psir_2;

    return result;
}
*/

#ifdef __LCAO
void toWannier90::lcao2pw_basis(const int ik,
                                const ModulePW::PW_Basis_K* wfcpw,
                                const Structure_Factor& sf,
                                ModuleBase::ComplexMatrix& orbital_in_G)
{
    this->table_local.create(GlobalC::ucell.ntype, GlobalC::ucell.nmax_total, GlobalV::NQX);
    Wavefunc_in_pw::make_table_q(GlobalC::ORB.orbital_file, this->table_local);
    Wavefunc_in_pw::produce_local_basis_in_pw(ik, wfcpw, sf, orbital_in_G, this->table_local);
}

void toWannier90::getUnkFromLcao(const ModulePW::PW_Basis_K* wfcpw,
                                 const Structure_Factor& sf,
                                 const K_Vectors& kv,
                                 const int npwx)
{
    std::complex<double> ***lcao_wfc_global = new std::complex<double> **[num_kpts];
    for (int ik = 0; ik < num_kpts; ik++)
    {
        lcao_wfc_global[ik] = new std::complex<double> *[GlobalV::NBANDS];
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            lcao_wfc_global[ik][ib] = new std::complex<double>[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(lcao_wfc_global[ik][ib], GlobalV::NLOCAL);
        }
    }

    if (this->unk_inLcao != nullptr)
    {
        delete this->unk_inLcao;
    }
    this->unk_inLcao = new psi::Psi<std::complex<double>>(num_kpts, GlobalV::NBANDS, npwx, nullptr);
    ModuleBase::ComplexMatrix *orbital_in_G = new ModuleBase::ComplexMatrix[num_kpts];

    for (int ik = 0; ik < num_kpts; ik++)
    {

        get_lcao_wfc_global_ik(lcao_wfc_global[ik], this->wfc_k_grid[ik]);

        int npw = kv.ngk[ik];
        orbital_in_G[ik].create(GlobalV::NLOCAL, npw);
        this->lcao2pw_basis(ik, wfcpw, sf, orbital_in_G[ik]);
    }

    for (int ik = 0; ik < num_kpts; ik++)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            for (int ig = 0; ig < kv.ngk[ik]; ig++)
            {
                for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
                {
                    unk_inLcao[0](ik, ib, ig) += orbital_in_G[ik](iw, ig) * lcao_wfc_global[ik][ib][iw];
                }
            }
        }
    }

    for (int ik = 0; ik < num_kpts; ik++)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            std::complex<double> anorm(0.0, 0.0);
            for (int ig = 0; ig < kv.ngk[ik]; ig++)
            {
                anorm = anorm + conj(unk_inLcao[0](ik, ib, ig)) * unk_inLcao[0](ik, ib, ig);
            }

            std::complex<double> anorm_tem(0.0, 0.0);
#ifdef __MPI
            MPI_Allreduce(&anorm, &anorm_tem, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif

            for (int ig = 0; ig < kv.ngk[ik]; ig++)
            {
                unk_inLcao[0](ik, ib, ig) = unk_inLcao[0](ik, ib, ig) / sqrt(anorm_tem);
            }
        }
    }

    for (int ik = 0; ik < kv.nkstot; ik++)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            delete[] lcao_wfc_global[ik][ib];
        }
        delete[] lcao_wfc_global[ik];
    }
    delete[] lcao_wfc_global;

    delete[] orbital_in_G;

#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    return;
}

void toWannier90::get_lcao_wfc_global_ik(std::complex<double> **ctot, std::complex<double> **cc)
{
    std::complex<double> *ctot_send = new std::complex<double>[GlobalV::NBANDS * GlobalV::NLOCAL];

#ifdef __MPI
    MPI_Status status;
#endif

    for (int i = 0; i < GlobalV::DSIZE; i++)
    {
        if (GlobalV::DRANK == 0)
        {
            if (i == 0)
            {
                // get the wave functions from 'ctot',
                // save them in the matrix 'c'.
                for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
                {
                    const int mu_local = this->gridt->trace_lo[iw];
                    if (mu_local >= 0)
                    {
                        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
                        {
                            // ctot[ib][iw] = cc[ib][mu_local];
                            ctot_send[ib * GlobalV::NLOCAL + iw] = cc[ib][mu_local];
                        }
                    }
                }
            }
            else
            {
                int tag;
                // receive lgd2
                int lgd2 = 0;
                tag = i * 3;
#ifdef __MPI
                MPI_Recv(&lgd2, 1, MPI_INT, i, tag, DIAG_WORLD, &status);
#endif
                if (lgd2 == 0)
                {
                }
                else
                {
                    // receive trace_lo2
                    tag = i * 3 + 1;
                    int *trace_lo2 = new int[GlobalV::NLOCAL];
#ifdef __MPI
                    MPI_Recv(trace_lo2, GlobalV::NLOCAL, MPI_INT, i, tag, DIAG_WORLD, &status);
#endif

                    // receive crecv
                    std::complex<double> *crecv = new std::complex<double>[GlobalV::NBANDS * lgd2];
                    ModuleBase::GlobalFunc::ZEROS(crecv, GlobalV::NBANDS * lgd2);
                    tag = i * 3 + 2;
#ifdef __MPI
                    MPI_Recv(crecv, GlobalV::NBANDS * lgd2, MPI_DOUBLE_COMPLEX, i, tag, DIAG_WORLD, &status);
#endif
                    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
                    {
                        for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
                        {
                            const int mu_local = trace_lo2[iw];
                            if (mu_local >= 0)
                            {
                                // ctot[ib][iw] = crecv[mu_local*GlobalV::NBANDS+ib];
                                ctot_send[ib * GlobalV::NLOCAL + iw] = crecv[mu_local * GlobalV::NBANDS + ib];
                            }
                        }
                    }

                    delete[] crecv;
                    delete[] trace_lo2;
                }
            }
        } // end GlobalV::DRANK=0
        else if (i == GlobalV::DRANK)
        {
            int tag;

            // send this->gridt->lgd
            tag = GlobalV::DRANK * 3;
#ifdef __MPI
            MPI_Send(&this->gridt->lgd, 1, MPI_INT, 0, tag, DIAG_WORLD);
#endif

            if (this->gridt->lgd != 0)
            {
                // send trace_lo
                tag = GlobalV::DRANK * 3 + 1;
#ifdef __MPI
                MPI_Send(this->gridt->trace_lo, GlobalV::NLOCAL, MPI_INT, 0, tag, DIAG_WORLD);
#endif

                // send cc
                std::complex<double> *csend = new std::complex<double>[GlobalV::NBANDS * this->gridt->lgd];
                ModuleBase::GlobalFunc::ZEROS(csend, GlobalV::NBANDS * this->gridt->lgd);

                for (int ib = 0; ib < GlobalV::NBANDS; ib++)
                {
                    for (int mu = 0; mu < this->gridt->lgd; mu++)
                    {
                        csend[mu * GlobalV::NBANDS + ib] = cc[ib][mu];
                    }
                }

                tag = GlobalV::DRANK * 3 + 2;
#ifdef __MPI
                MPI_Send(csend, GlobalV::NBANDS * this->gridt->lgd, MPI_DOUBLE_COMPLEX, 0, tag, DIAG_WORLD);
#endif

                delete[] csend;
            }
        } // end i==GlobalV::DRANK
#ifdef __MPI
        MPI_Barrier(DIAG_WORLD);
#endif
    }
#ifdef __MPI
    MPI_Bcast(ctot_send, GlobalV::NBANDS * GlobalV::NLOCAL, MPI_DOUBLE_COMPLEX, 0, DIAG_WORLD);
#endif

    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
        {
            ctot[ib][iw] = ctot_send[ib * GlobalV::NLOCAL + iw];
        }
    }

    delete[] ctot_send;

    return;
}

#endif
