#include "lbfgs.h"
#include "source_pw/hamilt_pwdft/global.h"
#include "source_base/matrix3.h"
#include "module_parameter/parameter.h"
#include "ions_move_basic.h"
#include "source_cell/update_cell.h"
#include "source_cell/print_cell.h" // mohan add 2025-06-19  

void LBFGS::allocate(const int _size) // initialize H0、H、pos0、force0、force
{
    alpha=70;//default value in ase is 70
    maxstep=PARAM.inp.relax_bfgs_rmax;
    size=_size;
    memory=100;
    iteration=0;
    H = std::vector<std::vector<double>>(3*size, std::vector<double>(3*size, 0.0));
    H0=1/alpha;
    pos = std::vector<std::vector<double>> (size, std::vector<double>(3, 0.0)); 
    pos0 = std::vector<double>(3*size, 0.0);
    pos_taud = std::vector<std::vector<double>> (size, std::vector<double>(3, 0.0)); 
    pos_taud0 = std::vector<double>(3*size, 0.0);
    dpos = std::vector<std::vector<double>>(size, std::vector<double>(3, 0.0));
    force0 = std::vector<double>(3*size, 0.0);
    force = std::vector<std::vector<double>>(size, std::vector<double>(3, 0.0));
    steplength = std::vector<double>(size, 0.0);
    //l_search.init_line_search();
}

void LBFGS::relax_step(const ModuleBase::matrix _force,UnitCell& ucell,const double &etot)

{
    get_pos(ucell,pos);  
    get_pos_taud(ucell,pos_taud);
    //solver=p_esolver;
    ucell.ionic_position_updated = true;
    for(int i = 0; i < _force.nr; i++)
    {
        for(int j=0;j<_force.nc;j++)
        {
            force[i][j]=_force(i,j)*ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A;
        }
    }
    int k=0;
    for(int i=0;i<ucell.ntype;i++)
    {
        for(int j=0;j<ucell.atoms[i].na;j++)
        {
            if(ucell.atoms[i].mbl[j].x==0)
            {
                force[k+j][0]=0;
            }
            if(ucell.atoms[i].mbl[j].y==0)
            {
                force[k+j][1]=0;
            }
            if(ucell.atoms[i].mbl[j].z==0)
            {
                force[k+j][2]=0;
            }
        }
        k+=ucell.atoms[i].na;
    }
    this->prepare_step(force,pos,H,pos0,force0,dpos,ucell,etot);
    this->determine_step(steplength,dpos,maxstep);
    this->update_pos(ucell);
    this->calculate_largest_grad(_force,ucell);
    this->is_restrain(dpos);  
    // mohan add 2025-06-22
    unitcell::print_tau(ucell.atoms,ucell.Coordinate,ucell.ntype,ucell.lat0,GlobalV::ofs_running);
}

void LBFGS::get_pos(UnitCell& ucell,std::vector<std::vector<double>>& pos)
{
    int k=0;
    for(int i=0;i<ucell.ntype;i++)
    {
        for(int j=0;j<ucell.atoms[i].na;j++)
        {
            pos[k+j][0]=ucell.atoms[i].tau[j].x*ModuleBase::BOHR_TO_A*ucell.lat0;
            pos[k+j][1]=ucell.atoms[i].tau[j].y*ModuleBase::BOHR_TO_A*ucell.lat0;
            pos[k+j][2]=ucell.atoms[i].tau[j].z*ModuleBase::BOHR_TO_A*ucell.lat0; 
        }
        k+=ucell.atoms[i].na;
    }
}

void LBFGS::get_pos_taud(UnitCell& ucell,std::vector<std::vector<double>>& pos_taud)
{
    int k=0;
    for(int i=0;i<ucell.ntype;i++)
    {
        for(int j=0;j<ucell.atoms[i].na;j++)
        {
            pos_taud[k+j][0]=ucell.atoms[i].taud[j].x;
            pos_taud[k+j][1]=ucell.atoms[i].taud[j].y;
            pos_taud[k+j][2]=ucell.atoms[i].taud[j].z;
        }
        k+=ucell.atoms[i].na;
    }
}

void LBFGS::prepare_step(std::vector<std::vector<double>>& force,
                         std::vector<std::vector<double>>& pos,
                         std::vector<std::vector<double>>& H,
                         std::vector<double>& pos0,
                         std::vector<double>& force0,
                         std::vector<std::vector<double>>& dpos,
                         UnitCell& ucell,
                         const double &etot)
{
    std::vector<double> changedforce = ReshapeMToV(force);
    std::vector<double> changedpos = ReshapeMToV(pos);
    this->update(pos_taud,pos_taud0,changedforce,force0,ucell,iteration,memory,s,y,rho);
    std::vector<double> q=DotInVAndFloat(changedforce,-1);
    int loopmax=std::min(memory,iteration);
    std::vector<double> a(loopmax);
    for(int i=loopmax-1;i>=0;i--)
    {
        a[i]=rho[i]*DotInVAndV(s[i],q);
        std::vector<double> temp=DotInVAndFloat(y[i],a[i]);
        q=VSubV(q,temp);
    }
    std::vector<double> z=DotInVAndFloat(q,H0);
    for(int i=0;i<loopmax;i++)
    {
        double b=rho[i]*DotInVAndV(y[i],z);
        std::vector<double> temp=DotInVAndFloat(s[i],a[i]-b);
        z=VAddV(z,temp);
    }
    std::vector<double> temp0=DotInVAndFloat(z,-1);
    dpos=ReshapeVToM(temp0);
    std::vector<double> temp1=DotInVAndFloat(changedforce,-1);
    std::vector<std::vector<double>> g=ReshapeVToM(temp1);
    energy=etot;
    //alpha_k=l_search.line_search(ucell,pos,g,energy,maxstep,size,dpos,pos,solver);
    //std::vector<double> temp2=DotInVAndFloat(temp0,alpha_k);
    std::vector<double> temp2=DotInVAndFloat(temp0,1);
    dpos=ReshapeVToM(temp2);
    for(int i = 0; i < size; i++)
    {
        double k = 0;
        for(int j = 0; j < 3; j++)
        {
            k += dpos[i][j] * dpos[i][j];
        }
        steplength[i] = sqrt(k);
    }
    iteration+=1;
    pos0 = ReshapeMToV(pos);
    pos_taud0=ReshapeMToV(pos_taud);
    force0 = changedforce;
}
void LBFGS::update(std::vector<std::vector<double>>& pos_taud, 
                   std::vector<double>& pos_taud0, 
                   std::vector<double>& force,
                   std::vector<double>& force0, 
                   UnitCell& ucell,
                   int iteration,
                   int memory,
                   std::vector<std::vector<double>>& s,
                   std::vector<std::vector<double>>& y,
                   std::vector<double>& rho)
{
    if(iteration>0)
    {
        std::vector<double> term=ReshapeMToV(pos_taud);
        std::vector<double> dpos =VSubV(term, pos_taud0);
        for(int i=0;i<3*size;i++)
        {
            double shortest_move = dpos[i];
            for (int cell = -1; cell <= 1; ++cell)
            {
                const double now_move = dpos[i] + cell;
                if (std::abs(now_move) < std::abs(shortest_move))
                {
                    shortest_move = now_move;
                }
            }
            dpos[i]=shortest_move;
        }
        std::vector<std::vector<double>> c=ReshapeVToM(dpos);
        for(int iat=0; iat<size; iat++)
        {
            //Cartesian coordinate
            //convert from Angstrom to unit of latvec (Bohr)

            //convert unit
            ModuleBase::Vector3<double> move_ion_cart;
            move_ion_cart.x = c[iat][0] *ModuleBase::BOHR_TO_A * ucell.lat0;
            move_ion_cart.y = c[iat][1] * ModuleBase::BOHR_TO_A * ucell.lat0;
            move_ion_cart.z = c[iat][2] * ModuleBase::BOHR_TO_A * ucell.lat0;

            //convert pos
            ModuleBase::Vector3<double> move_ion_dr = move_ion_cart* ucell.latvec;
            int it = ucell.iat2it[iat];
            int ia = ucell.iat2ia[iat];
            Atom* atom = &ucell.atoms[it];
            if(atom->mbl[ia].x == 1)
            {
                dpos[iat * 3] = move_ion_dr.x;
            }
            if(atom->mbl[ia].y == 1)
            {
                dpos[iat * 3 + 1] = move_ion_dr.y ;
            }
            if(atom->mbl[ia].z == 1)
            {
                dpos[iat * 3 + 2] = move_ion_dr.z ;
            }
        }
        std::vector<double> dforce =VSubV(force0, force);
        double rho0=1.0/DotInVAndV(dpos,dforce);
        s.push_back(dpos);
        y.push_back(dforce);
        rho.push_back(rho0);
    }

    if(iteration>memory)
    {
        s.erase(s.begin());
        y.erase(y.begin());
        rho.erase(rho.begin());
    }
}
void LBFGS::determine_step(std::vector<double>& steplength,std::vector<std::vector<double>>& dpos,double& maxstep)
{
    std::vector<double>::iterator maxsteplength = max_element(steplength.begin(), steplength.end());
    double a = *maxsteplength;
    if(a >= maxstep)
    {
        double scale = maxstep / a;
        for(int i = 0; i < size; i++)
        {
            for(int j=0;j<3;j++)
            {
                dpos[i][j]*=scale;
            }
        }
    }
}
void LBFGS::update_pos(UnitCell& ucell)
{
    double a[3*size];
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            a[i*3+j]=pos[i][j]+dpos[i][j];
            a[i*3+j]/=ModuleBase::BOHR_TO_A;
        }
    }
    unitcell::update_pos_tau(ucell.lat,a,ucell.ntype,ucell.nat,ucell.atoms);
}

void LBFGS::is_restrain(std::vector<std::vector<double>>& dpos)
{
    Ions_Move_Basic::converged = Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / 0.529177<PARAM.inp.force_thr_ev;
}

void LBFGS::calculate_largest_grad(const ModuleBase::matrix& _force,UnitCell& ucell)
{
    std::vector<double> grad= std::vector<double>(3*size, 0.0);
    int iat = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        Atom *atom = &ucell.atoms[it];
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            for (int ik = 0; ik < 3; ++ik)
            {
                if (atom->mbl[ia][ik])
                {
                    grad[3 * iat + ik] = -_force(iat, ik) * ucell.lat0;
                }
            }
            ++iat;
        }
    }
    Ions_Move_Basic::largest_grad = 0.0;
    for (int i = 0; i < 3*size; i++)
    {
        if (Ions_Move_Basic::largest_grad < std::abs(grad[i]))
        {
            Ions_Move_Basic::largest_grad = std::abs(grad[i]);
        }
    }
    Ions_Move_Basic::largest_grad /= ucell.lat0;
    if (PARAM.inp.out_level == "ie")
    {
        std::cout << " LARGEST GRAD (eV/A)  : " << Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / 0.5291772109
                  << std::endl;
    }

}

