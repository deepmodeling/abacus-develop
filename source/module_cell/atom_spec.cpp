#include "atom_spec.h"
#include "module_parameter/parameter.h"
#include <cstdlib>

Atom::Atom()
{
    na = 0;
    label = "\0";
    label_orb = "\0";
    nw = 0;
    nwl = 0;
    Rcut = 0.0; // pengfei Li 16-2-29
    type = 0;
    stapos_wf = 0;
    mass = 0.0;
}

Atom::~Atom()
{
}

void Atom::set_index(void)
{
    assert(nw != 0);
    this->iw2l.resize(nw, 0);
    this->iw2n.resize(nw, 0);
    this->iw2m.resize(nw, 0);
    this->iw2_ylm.resize(nw, 0);
    this->iw2_new.resize(nw, false); // bool array to check if the local orbital is new

    int index = 0;
    for (int L = 0; L <= nwl; L++)
    {
        assert(l_nchi[L] >= 0);
        for (int N = 0; N < l_nchi[L]; N++)
        {
            for (int m = 0; m < 2 * L + 1; m++)
            {
                iw2l[index] = L;
                iw2n[index] = N;
                iw2m[index] = m;
                iw2_ylm[index] = L * L + m;
                if (m == 0)
                {
                    iw2_new[index] = true;
                }
                else
                {
                    iw2_new[index] = false;
                }
                ++index;
            }
        }
    }
    return;
}

void Atom::print_Atom(std::ofstream& ofs)
{
    // OUT(ofs,"print_Atom()");
    ModuleBase::GlobalFunc::OUT(ofs, "label", label);
    ModuleBase::GlobalFunc::OUT(ofs, "type", type);
    ModuleBase::GlobalFunc::OUT(ofs, "na", na);
    ModuleBase::GlobalFunc::OUT(ofs, "nwl", nwl);
    ModuleBase::GlobalFunc::OUT(ofs, "Rcut", Rcut); // pengfei Li 16-2-29
    ModuleBase::GlobalFunc::OUT(ofs, "nw", nw);
    ModuleBase::GlobalFunc::OUT(ofs, "stapos_wf", stapos_wf);
    ModuleBase::GlobalFunc::OUT(ofs, "mass", mass);
    ofs << std::endl;

    output::printv31_d(ofs, "atom_position(cartesian)", tau.data(), na);
    /*
    for (int i = 0;i < na;i++)
    {
        ofs << std::setw(15) << this->tau[i].x
            << std::setw(15) << this->tau[i].y
            << std::setw(15) << this->tau[i].z << std::endl;
    }
    */
    ofs << std::endl;

    return;
}

#include "module_base/parallel_common.h"
#ifdef __MPI
void Atom::bcast_atom(void)
{
    Parallel_Common::bcast_int(type);
    Parallel_Common::bcast_int(na);
    Parallel_Common::bcast_int(nwl);
    Parallel_Common::bcast_double(Rcut); // pengfei Li 16-2-29
    Parallel_Common::bcast_int(nw);
    Parallel_Common::bcast_int(stapos_wf);
    Parallel_Common::bcast_string(label);
    Parallel_Common::bcast_bool(coulomb_potential);
    if (GlobalV::MY_RANK != 0)
    {
        this->l_nchi.resize(nwl + 1, 0);
    }
    Parallel_Common::bcast_int(l_nchi.data(), nwl + 1);
    Parallel_Common::bcast_bool(this->flag_empty_element);
    Parallel_Common::bcast_double(mass);

    if (GlobalV::MY_RANK != 0)
    {
        assert(na != 0);
        this->tau.resize(na, ModuleBase::Vector3<double>(0, 0, 0));
        this->dis.resize(na, ModuleBase::Vector3<double>(0, 0, 0));
        this->taud.resize(na, ModuleBase::Vector3<double>(0, 0, 0));
        this->vel.resize(na, ModuleBase::Vector3<double>(0, 0, 0));
        this->mag.resize(na, 0);
        this->angle1.resize(na, 0);
        this->angle2.resize(na, 0);
        this->m_loc_.resize(na, ModuleBase::Vector3<double>(0, 0, 0));
        this->mbl.resize(na, ModuleBase::Vector3<int>(0, 0, 0));
        this->lambda.resize(na, ModuleBase::Vector3<double>(0, 0, 0));
        this->constrain.resize(na, ModuleBase::Vector3<int>(0, 0, 0));
    }

    for (int i = 0; i < na; i++)
    {
        Parallel_Common::bcast_double(tau[i].x);
        Parallel_Common::bcast_double(tau[i].y);
        Parallel_Common::bcast_double(tau[i].z);
        Parallel_Common::bcast_double(taud[i].x);
        Parallel_Common::bcast_double(taud[i].y);
        Parallel_Common::bcast_double(taud[i].z);
        Parallel_Common::bcast_double(dis[i].x);
        Parallel_Common::bcast_double(dis[i].y);
        Parallel_Common::bcast_double(dis[i].z);
        Parallel_Common::bcast_double(vel[i].x);
        Parallel_Common::bcast_double(vel[i].y);
        Parallel_Common::bcast_double(vel[i].z);
        Parallel_Common::bcast_double(mag[i]);
        Parallel_Common::bcast_double(angle1[i]);
        Parallel_Common::bcast_double(angle2[i]);
        Parallel_Common::bcast_double(m_loc_[i].x);
        Parallel_Common::bcast_double(m_loc_[i].y);
        Parallel_Common::bcast_double(m_loc_[i].z);
        Parallel_Common::bcast_int(mbl[i].x);
        Parallel_Common::bcast_int(mbl[i].y);
        Parallel_Common::bcast_int(mbl[i].z);
        Parallel_Common::bcast_double(lambda[i].x);
        Parallel_Common::bcast_double(lambda[i].y);
        Parallel_Common::bcast_double(lambda[i].z);
        Parallel_Common::bcast_int(constrain[i].x);
        Parallel_Common::bcast_int(constrain[i].y);
        Parallel_Common::bcast_int(constrain[i].z);
    }

    return;
}

void Atom::bcast_atom2()
{
    this->ncpp.bcast_atom_pseudo();
}

#endif
