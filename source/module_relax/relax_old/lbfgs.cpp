#include "lbfgs.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/matrix3.h"
#include "module_parameter/parameter.h"
#include "ions_move_basic.h"

void LBFGS::allocate(const int _size) // initialize H0、H、pos0、force0、force
{
    alpha=70;//default value in ase is 70
    maxstep=PARAM.inp.relax_bfgs_rmax;
    size=_size;
    sign =true;
    memory=100;
    damping=70.0;
    iteration=0;
    xtol=1e-14;
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
    isave = std::vector<int>(2, 0);
    dsave = std::vector<double>(13, 0.0);
    old_stp=0;
    task=0;
}

void LBFGS::relax_step(const ModuleBase::matrix _force,UnitCell& ucell,const double &etot,ModuleESolver::ESolver* p_esolver) 
{
    GetPos(ucell,pos);  
    GetPostaud(ucell,pos_taud);
    solver=p_esolver;
    ucell.ionic_position_updated = true;
    for(int i = 0; i < _force.nr; i++)
    {
        for(int j=0;j<_force.nc;j++)
        {
            force[i][j]=_force(i,j)*ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A;
        }
    }
    /*std::cout<<"force"<<std::endl;
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<force[i][j]<<' ';
        }
        std::cout<<std::endl;
    }*/
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
    //std::cout<<"enterstep1"<<std::endl;
    this->PrepareStep(force,pos,H,pos0,force0,dpos,ucell,etot);
    this->DetermineStep(steplength,dpos,maxstep);
    //std::cout<<"enterstep2"<<std::endl;
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

void LBFGS::GetPos(UnitCell& ucell,std::vector<std::vector<double>>& pos)
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

void LBFGS::GetPostaud(UnitCell& ucell,std::vector<std::vector<double>>& pos_taud)
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

void LBFGS::PrepareStep
(std::vector<std::vector<double>>& force,
std::vector<std::vector<double>>& pos,
std::vector<std::vector<double>>& H,
std::vector<double>& pos0,
std::vector<double>& force0,
std::vector<std::vector<double>>& dpos,
UnitCell& ucell,
const double &etot)
{
    /*std::cout<<"force"<<std::endl;
    for(int i=0;i<force.size();i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<force[i][j]<<' ';
        }
        std::cout<<std::endl;
    }*/
    std::vector<double> changedforce = this->ReshapeMToV(force);
    std::vector<double> changedpos = this->ReshapeMToV(pos);
    this->Update(pos_taud,pos_taud0,changedforce,force0,ucell,iteration,memory,s,y,rho);
    /*std::cout<<'s'<<std::endl;
    for(int i=0;i<s.size();i++)
    {
        for(int j=0;j<s[i].size();j++)
        {
            std::cout<<s[i][j]<<' ';
        }
    std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<'y'<<std::endl;
    for(int i=0;i<y.size();i++)
    {
        for(int j=0;j<y[i].size();j++)
        {
            std::cout<<y[i][j]<<' ';
        }
    std::cout<<std::endl;
    }
    std::cout<<std::endl;*/
    std::vector<double> q=DotInVAndFloat(changedforce,-1);
    int loopmax=std::min(memory,iteration);
    std::vector<double> a(loopmax);
    for(int i=loopmax-1;i>=0;i--)
    {
        a[i]=rho[i]*this->DotInVAndV(s[i],q);
        auto temp=this->DotInVAndFloat(y[i],a[i]);
        q=this->VSubV(q,temp);
    }
    /*std::cout<<'q'<<std::endl;
    for(int i=0;i<q.size();i++)
    {
        std::cout<<q[i]<<' ';
    }
    std::cout<<std::endl;*/
    std::vector<double> z=this->DotInVAndFloat(q,H0);
    for(int i=0;i<loopmax;i++)
    {
        double b=rho[i]*this->DotInVAndV(y[i],z);
        auto temp=DotInVAndFloat(s[i],a[i]-b);
        z=VAddV(z,temp);
    }
    auto temp0=DotInVAndFloat(z,-1);
    dpos=ReshapeVToM(temp0);
    auto temp1=DotInVAndFloat(changedforce,-1);
    std::vector<std::vector<double>> g=ReshapeVToM(temp1);
    energy=etot;

    //alpha_k=LineSearch(ucell,pos,g,energy);


    //auto temp2=DotInVAndFloat(temp0,alpha_k);
    auto temp2=DotInVAndFloat(temp0,1);
    dpos=ReshapeVToM(temp2);
    /*std::cout<<"dpos"<<std::endl;
    for(int i=0;i<force.size();i++)
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
    iteration+=1;
    pos0 = this->ReshapeMToV(pos);
    pos_taud0=this->ReshapeMToV(pos_taud);
    force0 = changedforce;
}
void LBFGS::Update(std::vector<std::vector<double>>& pos_taud, 
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
        std::vector<double> dforce = this->VSubV(force0, force);
        double rho0=1.0/this->DotInVAndV(dpos,dforce);
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
void LBFGS::DetermineStep(std::vector<double>& steplength,
                         std::vector<std::vector<double>>& dpos,
                         double& maxstep)
{
    auto maxsteplength = max_element(steplength.begin(), steplength.end());
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
double LBFGS::LineSearch(UnitCell& ucell,std::vector<std::vector<double>>& r,std::vector<std::vector<double>>& g,double e)
{
    std::vector<double> tempp=ReshapeMToV(dpos);
    double p_size=0;
    for(int i=0;i<tempp.size();i++)
    {
        p_size+=tempp[i]*tempp[i];
    }
    p_size=sqrt(p_size);
    if(p_size<sqrt(size*1e-10))
    {
        for(int i=0;i<tempp.size();i++)
        {
            tempp[i]/=(p_size/sqrt(size*1e-10));
        }
    }
    std::vector<double> tempr=ReshapeMToV(r);
    std::vector<double> tempg=ReshapeMToV(g);
    double c1=0.23;
    double c2=0.46;
    stpmax=50;
    stpmin=1e-8;
    xtrapl=1.1;
    xtrapu=4;
    no_update=false;
    double phi0=e;
    double derphi0=DotInVAndV(tempg,tempp);
    double gms=maxstep*sqrt(3*size);
    double alpha1=1;
    bool no_update=false;
    double fval=e;
    double stp=0;
    std::vector<double> gval=tempg;
    while(true)
    {
        stp=Step(alpha1,phi0,derphi0,c1,c2,xtol,isave,dsave);
        if(task==1)
        {
            alpha1=stp;
            fval=GetEnergy(ucell,stp);
            gval=GetForce(ucell,stp);
            phi0=fval;
            auto temp=ReshapeMToV(dpos);
            derphi0=DotInVAndV(gval,temp);
            old_stp=alpha1;
            if(no_update)
            {
                break;
            }
        }
        else
        {
            break;
        }
    }
    return stp;

}
double LBFGS::Step(double stp,double f,double g,double c1,double c2,double xtol,std::vector<int>& isave,std::vector<double>& dsave)
{
    if(task==0)
    {
        if(stp<stpmin)
        {
            std::cout<<"ERROR: STP LT minstep"<<std::endl;
            task=2;
        }
        if(stp>stpmax)
        {
            std::cout<<"ERROR: STP GT maxstep"<<std::endl;
            task=2;
        }
        if(g>=0)
        {
            std::cout<<"ERROR: INITIAL G >= 0"<<std::endl;
            task=2;
        }
        if(c1<0)
        {
            std::cout<<"ERROR: c1 LT 0"<<std::endl;
            task=2;
        }
        if(c2<0)
        {
            std::cout<<"ERROR: c2 LT 0"<<std::endl;
            task=2;
        }
        if(xtol<0)
        {
            std::cout<<"ERROR: XTOL LT 0"<<std::endl;
            task=2;
        }
        if(stpmin<0)
        {
            std::cout<<"ERROR: minstep LT 0"<<std::endl;
            task=2;
        }
        if(stpmax<stpmin)
        {
            std::cout<<"ERROR: maxstep LT minstep"<<std::endl;
            task=2;
        }
        if(task==2)
        {
            return stp;
        }
        bracket=false;
        int stage=1;
        double finit=f;
        double ginit=g;
        double gtest=c1*ginit;
        double width=stpmax-stpmin;
        double width1=width/0.5;
        double stx=0;
        double fx=finit;
        double gx=ginit;
        double sty=0;
        double fy=finit;
        double gy=ginit;
        double stmin=0;
        double stmax=stp+xtrapu*stp;
        task=1;
        Save(stage,ginit,gtest,gx,gy,finit,fx,fy,stx,sty,stmin,stmax,width,width1);
        stp=DetermineStep(stp);
    }
    else
    {
        if(isave[0]==1)
        {
            bracket=true;
        }
        else
        {
            bracket=false;
        }
        int stage=isave[1];
        double ginit=dsave[0];
        double gtest=dsave[1];
        double gx=dsave[2];
        double gy=dsave[3];
        double finit=dsave[4];
        double fx=dsave[5];
        double fy=dsave[6];
        double stx=dsave[7];
        double sty=dsave[8];
        double stmin=dsave[9];
        double stmax=dsave[10];
        double width=dsave[11];
        double width1=dsave[12];

        double ftest=finit+stp*gtest;
        if(stage==1 && f<ftest && g>=0)
        {
            stage=2;
        }
        if(bracket && (stp<=stmin || stp>=stmax))
        {
            std::cout<<"WARNING: ROUNDING ERRORS PREVENT PROGRESS"<<std::endl;
            task=2;
        }
        if(bracket&&(stmax-stmin)<=xtol*stmax)
        {
            std::cout<<"WARNING: XTOL TEST SATISFIED"<<std::endl;
            task=2;
        }
        if(stp==stpmax&&f<=ftest&&g<=gtest)
        {
            std::cout<<"WARNING: STP=maxstep"<<std::endl;
            task=2;
        }
        if(stp==stpmin&&(f>ftest||g>=gtest))
        {
            std::cout<<"WARNING::STP=minstep"<<std::endl;
            task=2;
        }
        if(f<=ftest&&std::abs(g)<=c2*(-ginit))
        {
            task=3;
        }     
        if(task==2||task==3)
        {
            Save(stage,ginit,gtest,gx,gy,finit,fx,fy,stx,sty,stmin,stmax,width,width1);
            return stp;
        }
        UpdateLineSearch(stx,fx,gx,sty,fy,gy,stp,f,g,stmin,stmax);
        if(bracket)
        {
            if(std::abs(sty-stx)>=0.66*width1)
            {
                stp=stx+0.5*(sty-stx);
            }
            width1=width;
            width=std::abs(sty-stx);
            stmin=std::min(stx,sty);
            stmax=std::max(stx,sty);
        }
        else
        {
            stmin=stp+xtrapl*(stp-stx);
            stmax=stp+xtrapu*(stp-stx);
        }
        stp=std::max(stp,stpmin);
        stp=std::min(stp,stpmax);
        if(stx==stp&&stp==stpmax&&stmin>stpmax)
        {
            no_update=true;
        }
        if((bracket&&stp<stmin)||stp>=stmax||(bracket&&(stmax-stmin)<xtol*stmax))
        {
            stp=stx;
        }
        task=1;
        Save(stage,ginit,gtest,gx,gy,finit,fx,fy,stx,sty,stmin,stmax,width,width1);
    }
    return stp;
}

void LBFGS::UpdateLineSearch(double& stx,double& fx,double& gx,double& sty,double& fy,double& gy,double& stp,double& fp,double& gp,double& stpmin,double& stpmax)
{
    double sign=gp*(gx/std::abs(gx));
    double stpf=0;
    double theta=0;
    double s=0;
    double gamma=0;
    double p=0;
    double q=0;
    double r=0;
    double stpc=0;
    double stpq=0;
    if(fp>fx)
    {
        theta=3*(fx-fp)/(stp-stx)+gx+gp;
        s=std::max({std::abs(theta),std::abs(gx),std::abs(gp)});
        gamma=s*sqrt((theta*theta-gx*gp)/(s*s));
        if(stp<stx)
        {
            gamma=-gamma;
        }
        p=gamma-gx+theta;
        q=gamma-gx+gamma+gp;
        r=p/q;
        stpc=stx+r*(stp-stx);
        stpq=stx+((gx/(((fx-fp)/(stp-stx))+gx))/2.0)*(stp-stx);
        if(abs(stpc-stx)<abs(stpq-stx))
        {
            stpf=stpc;
        }
        else
        {
            stpf = stpc + (stpq - stpc) / 2.0;
        }
        bracket = true;
    }
    else if (sign < 0) 
    {
        theta = 3.0 * (fx - fp) / (stp - stx) + gx + gp;
        s = std::max({std::abs(theta), std::abs(gx), std::abs(gp)});
        gamma = s * sqrt(std::pow(theta / s, 2.0) - (gx / s) * (gp / s));
        if (stp > stx) 
        {
            gamma = -gamma;
        }
        p = (gamma - gp) + theta;
        q = ((gamma - gp) + gamma) + gx;
        r = p / q;
        stpc = stp + r * (stx - stp);
        stpq = stp + (gp / (gp - gx)) * (stx - stp);
        if (abs(stpc - stp) > abs(stpq - stp))
        {
            stpf = stpc;
        }       
        else
        {
            stpf = stpq;
        }
        bracket = true;
    }
    else if (abs(gp) < abs(gx)) 
    {
        theta = 3.0 * (fx - fp) / (stp - stx) + gx + gp;
        s = std::max({std::abs(theta), std::abs(gx), std::abs(gp)});
        gamma = s * sqrt(std::max(0.0, std::pow(theta / s, 2.0) - (gx / s) * (gp / s)));
        if (stp > stx)
        {
            gamma = -gamma;
        }
        p = (gamma - gp) + theta;
        q = (gamma + (gx - gp)) + gamma;
        r = p / q;
        if (r < 0.0 && gamma != 0.0)
        {
            stpc = stp + r * (stx - stp);
        }      
        else if (stp > stx)
        {
            stpc = stpmax;
        }           
        else
        {
            stpc = stpmin;
        } 
        stpq = stp + (gp / (gp - gx)) * (stx - stp);
        if (bracket) 
        {
            if (abs(stpc - stp) < abs(stpq - stp))
            {
                stpf = stpc;
            }           
            else
            {
                stpf = stpq;
            }
            if (stp > stx)
            {
                stpf = std::min(stp + 0.66 * (sty - stp), stpf);
            }
            else
            {
                stpf = std::max(stp + 0.66 * (sty - stp), stpf);
            }
        }
        else
        {
            if (abs(stpc - stp) > abs(stpq - stp))
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
            stpf = std::min(stpmax, stpf);
            stpf = std::max(stpmin, stpf);
        }
    }
    else
    {
        if(bracket) 
        {
            theta = 3.0 * (fp - fy) / (sty - stp) + gy + gp;
            s = std::max({std::abs(theta), std::abs(gy), std::abs(gp)});
            gamma = s * sqrt(pow(theta / s, 2.0) - (gy / s) * (gp / s));
            if (stp > sty)
            {
                gamma = -gamma;
            }
            p = (gamma - gp) + theta;
            q = ((gamma - gp) + gamma) + gy;
            r = p / q;
            stpc = stp + r * (sty - stp);
            stpf = stpc;
        }
        else if (stp > stx)
        {
            stpf = stpmax;
        }
        else
        {
            stpf = stpmin;
        }
    }
    if (fp > fx)
    {
        sty = stp;
        fy = fp;
        gy = gp;
    }
    else
    {
        if (sign < 0)
        {
            sty = stx;
            fy = fx;
            gy = gx;
        }
        stx = stp;
        fx = fp;
        gx = gp;
    }
    stp = DetermineStep(stpf);
}

double LBFGS::DetermineStep(double stp)
{
    double dr=stp-old_stp;
    for(int i = 0; i < size; i++)
    {
        double k = 0;
        for(int j = 0; j < 3; j++)
        {
            k += dpos[i][j] * dpos[i][j]*dr*dr;
        }
        steplength[i] = sqrt(k);
    }
    auto maxsteplength = max_element(steplength.begin(), steplength.end());
    double a = *maxsteplength;
    if(a >= maxstep)
    {
        double scale = maxstep / a;
        dr*=scale;
    }
    double k=old_stp+dr;
    return k;
}
void LBFGS::Save(int a,double b,double c,double d,double e,double f,double g,double h,double i,double j,double k,double l,double m,double n)
{
    if(bracket)
    {
        isave[0]=1;
    }
    else
    {
        isave[0]=0;
    }
    isave[1]=a;
    dsave[0]=b;
    dsave[1]=c;
    dsave[2]=d;
    dsave[3]=e;
    dsave[4]=f;
    dsave[5]=g;
    dsave[6]=h;
    dsave[7]=i;
    dsave[8]=j;
    dsave[9]=k;
    dsave[10]=l;
    dsave[11]=m;
    dsave[12]=n;
}



double LBFGS::GetEnergy(UnitCell& ucell,double stp)
{
    double a[3*size];
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            a[i*3+j]=pos[i][j]+stp*dpos[i][j];
            a[i*3+j]/=ModuleBase::BOHR_TO_A;
        }
    }
    int k=0;
    ucell.update_pos_tau(a);
    return solver->cal_energy();
}
std::vector<double> LBFGS::GetForce(UnitCell& ucell,double stp)
{
    double a[3*size];
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            a[i*3+j]=pos[i][j]+stp*dpos[i][j];
            a[i*3+j]/=ModuleBase::BOHR_TO_A;
        }
    }
    int k=0;
    ucell.update_pos_tau(a);
    ModuleBase::matrix b;
    solver->cal_force(ucell,b);
    std::vector<double> c=std::vector<double>(3*size, 0.0);
    for(int i = 0; i < b.nr; i++)
    {
        for(int j=0;j<b.nc;j++)
        {
            c[i*b.nr+j]=b(i,j)*ModuleBase::Ry_to_eV/ModuleBase::BOHR_TO_A;
        }
    }
    return c;

}
void LBFGS::UpdatePos(UnitCell& ucell)
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
    /*std::cout<<"a"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        std::cout<<a[i]<<std::endl;
    }*/
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

void LBFGS::IsRestrain(std::vector<std::vector<double>>& dpos)
{
    Ions_Move_Basic::converged = Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / 0.529177<PARAM.inp.force_thr_ev;
}

void LBFGS::CalculateLargestGrad(const ModuleBase::matrix& _force,UnitCell& ucell)
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


// matrix methods
std::vector<double> LBFGS::ReshapeMToV(std::vector<std::vector<double>>& matrix) 
{
    int size = matrix.size();
    std::vector<double> result;
    result.reserve(3*size);
    for (const auto& row : matrix) {
        result.insert(result.end(), row.begin(), row.end());
    }
    return result;
}

std::vector<std::vector<double>> LBFGS::MAddM(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b) 
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

std::vector<double> LBFGS::VSubV(std::vector<double>& a, std::vector<double>& b) 
{
    std::vector<double> result = std::vector<double>(a.size(), 0.0);
    for(int i = 0; i < a.size(); i++)
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<double> LBFGS::VAddV(std::vector<double>& a, std::vector<double>& b) 
{
    std::vector<double> result = std::vector<double>(a.size(), 0.0);
    for(int i = 0; i < a.size(); i++)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<std::vector<double>> LBFGS::ReshapeVToM(std::vector<double>& matrix) 
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

std::vector<double> LBFGS::DotInMAndV1(std::vector<std::vector<double>>& matrix, std::vector<double>& vec) 
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
std::vector<double> LBFGS::DotInMAndV2(std::vector<std::vector<double>>& matrix, std::vector<double>& vec) 
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

double LBFGS::DotInVAndV(std::vector<double>& vec1, std::vector<double>& vec2) 
{
    double result = 0.0;
    for(int i = 0; i < vec1.size(); i++)
    {
        result += vec1[i] * vec2[i];
    }
    return result;
}

std::vector<std::vector<double>> LBFGS::OuterVAndV(std::vector<double>& a, std::vector<double>& b) 
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

std::vector<std::vector<double>> LBFGS::MPlus(std::vector<std::vector<double>>& a, double b)
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

std::vector<std::vector<double>> LBFGS::MSubM(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b)
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

std::vector<double> LBFGS::DotInVAndFloat(std::vector<double>& vec, double b) 
{
    std::vector<double> result(vec.size(), 0.0);
    for(int i = 0; i < vec.size(); i++)
    {
        result[i] = vec[i] * b;
    }
    return result;
}