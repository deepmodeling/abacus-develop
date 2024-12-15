#include "bfgs.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/matrix3.h"
#include "module_parameter/parameter.h"
#include "ions_move_basic.h"

//! initialize H0、H、pos0、force0、force
void BFGS::allocate(const int _size) 
{
    alpha=70;//default value in ase is 70
    maxstep=PARAM.inp.relax_bfgs_rmax;
    size=_size;
    sign =true;
    
    H = std::vector<std::vector<double>>(3*size, std::vector<double>(3*size, 0.0));
    
    for (int i = 0; i < 3*size; ++i) 
    {
        H[i][i] = alpha;  
    }
    
    pos = std::vector<std::vector<double>> (size, std::vector<double>(3, 0.0)); 
    pos0 = std::vector<double>(3*size, 0.0);
    pos_taud = std::vector<std::vector<double>> (size, std::vector<double>(3, 0.0)); 
    pos_taud0 = std::vector<double>(3*size, 0.0);
    dpos = std::vector<std::vector<double>>(size, std::vector<double>(3, 0.0));
    force0 = std::vector<double>(3*size, 0.0);
    force = std::vector<std::vector<double>>(size, std::vector<double>(3, 0.0));
    steplength = std::vector<double>(size, 0.0);  
}


void BFGS::relax_step(const ModuleBase::matrix& _force,UnitCell& ucell) 

{
    GetPos(ucell,pos);  
    GetPostaud(ucell,pos_taud);
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
    
    this->PrepareStep(force,pos,H,pos0,force0,steplength,dpos,ucell);
    this->DetermineStep(steplength,dpos,maxstep);
    
    /*std::cout<<"force"<<std::endl;
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<force[i][j]<<' ';
        }
        std::cout<<std::endl;
    }
    std::cout<<"dpos"<<std::endl;
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<dpos[i][j]<<' ';
        }
        std::cout<<std::endl;
    }
    std::cout<<"pos"<<std::endl;
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<pos[i][j]<<' ';
        }
        std::cout<<std::endl;
    }*/
    
    this->UpdatePos(ucell);
    this->CalculateLargestGrad(_force,ucell);
    this->IsRestrain(dpos);  
}

void BFGS::GetPos(UnitCell& ucell,std::vector<std::vector<double>>& pos)
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

void BFGS::GetPostaud(UnitCell& ucell,
                      std::vector<std::vector<double>>& pos_taud)
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

void BFGS::PrepareStep(std::vector<std::vector<double>>& force,
                       std::vector<std::vector<double>>& pos,
                       std::vector<std::vector<double>>& H,
                       std::vector<double>& pos0,
                       std::vector<double>& force0,
                       std::vector<double>& steplength,
                       std::vector<std::vector<double>>& dpos,
                       UnitCell& ucell)
{
    std::vector<double> changedforce = this->ReshapeMToV(force);
    std::vector<double> changedpos = this->ReshapeMToV(pos);
    this->Update(changedpos, changedforce,H,ucell);
    
    //! call dysev
    std::vector<double> omega(3*size);
    std::vector<double> work(3*size*3*size);
    int lwork=3*size*3*size;
    int info=0;
    std::vector<double> H_flat;
    
    for(const auto& row : H)
    {
        H_flat.insert(H_flat.end(), row.begin(), row.end());
    }
    
    int value=3*size;
    int* ptr=&value;
    dsyev_("V","U",ptr,H_flat.data(),ptr,omega.data(),work.data(),&lwork,&info);
    std::vector<std::vector<double>> V(3*size, std::vector<double>(3*size, 0.0));
    for(int i = 0; i < 3*size; i++)
    {
        for(int j = 0; j < 3*size; j++)
        {
            V[j][i] = H_flat[3*size*i + j];
        }
    }
    std::vector<double> a=this->DotInMAndV2(V, changedforce);
    for(int i = 0; i < a.size(); i++)
    {
        a[i]/=std::abs(omega[i]);    
    }
    std::vector<double> tmpdpos = this->DotInMAndV1(V, a);
    dpos = this->ReshapeVToM(tmpdpos);
    /*std::cout<<"dpos0"<<std::endl;
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<dpos[i][j]<<' ';
        }
        std::cout<<std::endl;
    }*/
    for(int i = 0; i < size; i++)
    {
        double k = 0;
        for(int j = 0; j < 3; j++)
        {
            k += dpos[i][j] * dpos[i][j];
        }
        steplength[i] = sqrt(k);
    }
    pos0 = this->ReshapeMToV(pos);
    pos_taud0=this->ReshapeMToV(pos_taud);
    force0 = this->ReshapeMToV(force);
}

void BFGS::Update(std::vector<double>& pos, 
                  std::vector<double>& force,
                  std::vector<std::vector<double>>& H,
                  UnitCell& ucell)
{
    if(sign)
    {
        sign=false;
        return;
    }
    //std::vector<double> dpos=this->VSubV(pos,pos0);
    auto term=this->ReshapeMToV(pos_taud);
    std::vector<double> dpos = this->VSubV(term, pos_taud0);
    for(int i=0;i<3*size;i++)
    {
        double shortest_move = dpos[i];
        //dpos[i]/=ModuleBase::BOHR_TO_A;
        //dpos[i]/=ucell.lat0;
        for (int cell = -1; cell <= 1; ++cell)
        {
            const double now_move = dpos[i] + cell;
            if (std::abs(now_move) < std::abs(shortest_move))
            {
                shortest_move = now_move;
            }
        }
        //shortest_move=shortest_move*ModuleBase::BOHR_TO_A*ucell.lat0;
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
    /*std::cout<<"Printpos"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        std::cout<<pos[i]<<' ';
    }
    std::cout<<std::endl;
    std::cout<<"Printpos0"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        std::cout<<pos0[i]<<' ';
    }
    std::cout<<std::endl;*/
    /*std::cout<<"PrintDpos"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        std::cout<<dpos[i]<<' ';
    }
    std::cout<<std::endl;*/
    if(*max_element(dpos.begin(), dpos.end()) < 1e-7)
    {
        return;
    }
    
    std::vector<double> dforce = this->VSubV(force, force0);
    double a = this->DotInVAndV(dpos, dforce);
    std::vector<double> dg = this->DotInMAndV1(H, dpos);
    double b = this->DotInVAndV(dpos, dg);
    
    /*std::cout<<"a"<<std::endl;
    std::cout<<a<<std::endl;
    std::cout<<"b"<<std::endl;
    std::cout<<b<<std::endl;*/
    auto term1=this->OuterVAndV(dforce, dforce);
    auto term2=this->OuterVAndV(dg, dg);
    auto term3=this->MPlus(term1, a);
    auto term4=this->MPlus(term2, b);
    H = this->MSubM(H, term3);
    H = this->MSubM(H, term4);
}

void BFGS::DetermineStep(std::vector<double>& steplength,
                         std::vector<std::vector<double>>& dpos,
                         double& maxstep)
{
    auto maxsteplength = max_element(steplength.begin(), steplength.end());
    double a = *maxsteplength;
    /*std::cout<<"maxstep"<<std::endl;
    std::cout<<maxstep<<std::endl;
    std::cout<<"maxsteplength"<<std::endl;
    std::cout<<a<<std::endl;*/
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

void BFGS::UpdatePos(UnitCell& ucell)
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
    std::cout<<std::endl;
    int k=0;
    ucell.update_pos_tau(a);
    /*double move_ion[3*size];
    ModuleBase::zeros(move_ion, size*3);

    for(int iat=0; iat<size; iat++)
    {
        //Cartesian coordinate
        //convert from Angstrom to unit of latvec (Bohr)

        //convert unit
        ModuleBase::Vector3<double> move_ion_cart;
        move_ion_cart.x = dpos[iat][0] / ModuleBase::BOHR_TO_A / ucell.lat0;
        move_ion_cart.y = dpos[iat][1] / ModuleBase::BOHR_TO_A / ucell.lat0;
        move_ion_cart.z = dpos[iat][2] / ModuleBase::BOHR_TO_A / ucell.lat0;

        //convert to Direct coordinate
        //note here the old GT is used

        //convert pos
        ModuleBase::Vector3<double> move_ion_dr = move_ion_cart* ucell.GT;

        int it = ucell.iat2it[iat];
        int ia = ucell.iat2ia[iat];
        Atom* atom = &ucell.atoms[it];

        if(atom->mbl[ia].x == 1)
        {
            move_ion[iat * 3] = move_ion_dr.x;
        }
        if(atom->mbl[ia].y == 1)
        {
            move_ion[iat * 3 + 1] = move_ion_dr.y ;
        }
        if(atom->mbl[ia].z == 1)
        {
            move_ion[iat * 3 + 2] = move_ion_dr.z ;
        }
    }
	ucell.update_pos_taud(move_ion);
    pos = this->MAddM(pos, dpos);*/
}

void BFGS::IsRestrain(std::vector<std::vector<double>>& dpos)
{
    Ions_Move_Basic::converged = Ions_Move_Basic::largest_grad 
        * ModuleBase::Ry_to_eV / 0.529177<PARAM.inp.force_thr_ev;
}

void BFGS::CalculateLargestGrad(const ModuleBase::matrix& _force,UnitCell& ucell)
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
        std::cout << " LARGEST GRAD (eV/A)  : " << Ions_Move_Basic::largest_grad 
            * ModuleBase::Ry_to_eV / 0.5291772109
                  << std::endl;
    }

}
// matrix methods

std::vector<double> BFGS::ReshapeMToV(std::vector<std::vector<double>>& matrix) 
{
    int size = matrix.size();
    std::vector<double> result;
    result.reserve(3*size);
    for (const auto& row : matrix) {
        result.insert(result.end(), row.begin(), row.end());
    }
    return result;
}

std::vector<std::vector<double>> BFGS::MAddM(std::vector<std::vector<double>>& a, 
                                             std::vector<std::vector<double>>& b) 
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(a[0].size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

std::vector<double> BFGS::VSubV(std::vector<double>& a, std::vector<double>& b) 
{
    std::vector<double> result = std::vector<double>(a.size(), 0.0);
    for(int i = 0; i < a.size(); i++)
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<std::vector<double>> BFGS::ReshapeVToM(std::vector<double>& matrix) 
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(matrix.size() / 3, std::vector<double>(3));
    for(int i = 0; i < result.size(); i++)
    {
        for(int j = 0; j < 3; j++)
        {
            result[i][j] = matrix[i*3 + j];
        }
    }
    return result;
}

std::vector<double> BFGS::DotInMAndV1(std::vector<std::vector<double>>& matrix, std::vector<double>& vec) 
{
    std::vector<double> result(matrix.size(), 0.0);
    for(int i = 0; i < result.size(); i++)
    {
        for(int j = 0; j < vec.size(); j++)
        {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}
std::vector<double> BFGS::DotInMAndV2(std::vector<std::vector<double>>& matrix, std::vector<double>& vec) 
{
    std::vector<double> result(matrix.size(), 0.0);
    for(int i = 0; i < result.size(); i++)
    {
        for(int j = 0; j < vec.size(); j++)
        {
            result[i] += matrix[j][i] * vec[j];
        }
    }
    return result;
}

double BFGS::DotInVAndV(std::vector<double>& vec1, std::vector<double>& vec2) 
{
    double result = 0.0;
    for(int i = 0; i < vec1.size(); i++)
    {
        result += vec1[i] * vec2[i];
    }
    return result;
}

std::vector<std::vector<double>> BFGS::OuterVAndV(std::vector<double>& a, std::vector<double>& b) 
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(b.size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < b.size(); j++)
        {
            result[i][j] = a[i] * b[j];
        }
    }
    return result;
}

std::vector<std::vector<double>> BFGS::MPlus(std::vector<std::vector<double>>& a, double& b)
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(a[0].size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            result[i][j] = a[i][j] / b;
        }
    }
    return result;
}

std::vector<std::vector<double>> BFGS::MSubM(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b)
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(a[0].size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            result[i][j] = a[i][j] - b[i][j];
        }
    }
    return result;
}
