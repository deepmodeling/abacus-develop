#include "line_search.h"
#include <cmath>
#include <algorithm>

#include "module_hamilt_pw/hamilt_pwdft/global.h"

bool Line_Search::line_search(
    const bool restart,
    const double x, //current point
    const double y, //function value at current point
    const double f, //derivative at current point
    double & xnew,  //the next point that we want to try
    const double conv_thr)
{
    if(restart) ls_step = 0;

    if(ls_step == 0) //first point: make a trial step into trial direction
    {
        return this->first_order(x,y,f,xnew);
    }
    if(ls_step == 1) //second point: third order extrapolation/interpolation
    {
        return this->third_order(x,y,f,xnew,conv_thr);
    }

    //if still not converging after 2 steps: start Brent
    if(ls_step == 2)
    {
        this->init_brent(x,y,f);
    }
    if(ls_step >= 3)
    {
        this->update_brent(x,y,f);
    }

    if(ls_step >= 2)
    {
        return this->brent(x,y,f,xnew,conv_thr);
    }
    ModuleBase::WARNING_QUIT("line_search","ls_step <0");
    __builtin_unreachable();
}

bool Line_Search::first_order(
    const double x,
    const double y,
    const double f,
    double & xnew)
{
    xa = x; // set the first endpoint
    ya = y;
    fa = f;
    xnew = x + 1; // move one step towards trial direction
    ls_step ++;
    return false;
}

bool Line_Search::third_order(
    const double x,
    const double y,
    const double f,
    double & xnew,
    const double conv_thr)
{
    double dmove, dmoveh, dmove1, dmove2;

    xb = x; // set the second endpoint
    yb = y;
    fb = f;

    double ffin = yb - ya;
    double fab = (fa + fb)/2.0;

    double k3 = 3.0*(fb    +fa-2.0*ffin);
    double k2 =      fb+2.0*fa-3.0*ffin;
    double k1 =      fa;

    double tmp = k1*k3/k2/k2;

    //cubic extrapolation
    if( std::abs(k3/k1)<1.0e-2 || std::abs(k3) < 2.4e-5 || tmp>1.0 ) //harmonic case
    {
        dmove =  -fa/(fab-fa)/2.0;
        if(dmove<0)
        {
            dmove  = 4.0;
        }
    }
    else //anharmonic case
    {
        dmove1 = k2/k3*(1.0-std::sqrt(1.0-tmp));
        dmove2 = k2/k3*(1.0+std::sqrt(1.0-tmp));

        double dy1, dy2;
        dy1 = -(k1-(k2-k3*dmove1/3.0)*dmove1)*dmove1;
        dy2 = -(k1-(k2-k3*dmove2/3.0)*dmove2)*dmove2;

        if(dy1>dy2)
        {
            dmove = dmove1;
        }
        else
        {
            dmove = dmove2;
        }

        dmoveh = -fa/(fab-fa)/2.0;
        if(dmoveh<0) dmoveh = 4.0;

        if(dmove > 2.0*dmoveh || dmoveh > 2.0*dmove ||
            (fa*fb > 0 && dmove < 1.0) ||
            (fa*fb < 0 && dmove > 1.0) )
        {
            dmove = dmoveh;
        }
    } //end anharmonic case

    if(dmove > 4.0) dmove = 4.0;
    xnew  = dmove +xa;

    double dy = (fb+(fab-fb)/(xa-xb)*(dmove-xb))*(dmove-xb);
    if(std::abs(dy)<conv_thr) return true;

    ls_step ++;
    return false;
}

void Line_Search::init_brent(
    const double x,
    const double y,
    const double f)
{
    bracked = true;

    if(x>xb) // x > b, start interval [b,x]
    {
        xa = xb;
        ya = yb;
        fa = fb;
        xb = x;
        yb = y;
        fb = f;
        if(fa*fb>0) bracked = false;
        fstart = fa;
    }
    else // x < b
    {
        if(fa*f<=0) // minimum between [a,x]
        {
            xb = x;
            yb = y;
            fb = f;
        }
        else if (fb*f<=0) // minimum between [x,b]
        {
            xa = xb;
            ya = yb;
            fa = fb;
            xb = x;
            yb = y;
            fb = f;
        }
        else //problematic case, no minimum between [a,b]
        {
            xa = xb;
            ya = yb;
            fa = fb;
            xb = x;
            yb = y;
            fb = f;
            bracked = false;
        }
    }

    xc = xb;
    fc = fb;
}

void Line_Search::update_brent(
    const double x,
    const double y,
    const double f)
{
    xb = x;
    yb = y;
    fb = f;
    if(!bracked && fstart*f<0)
    {
        bracked = true;
        xc = xb;
        fc = fb;
    }
}

bool Line_Search::brent(
    const double x,
    const double y,
    const double f,
    double & xnew,
    const double conv_thr)
{
    ls_step ++;

    double xd,xe,xm;
    // if no zero is between xa and xb
    if(!bracked)
    {
        if(std::abs(fc)<=std::abs(fb) || (xa-xc)<(xb-xa))
        {
            xc = xa;
            fc = fa;
            xd = xb-xa;
            xe = xb-xa;
        }
        double tol1 = 2.0*e8*std::abs(xb)+0.5*e8;
        xm = 0.5*(xc-xb);

        if(!(xc<=xa && xa<=xb))
        {
            ModuleBase::WARNING_QUIT("Brent","something wrong with Brent line search!");
        }
        if(std::abs(xm)<=tol1 || fb==0.0)
        {
            return true;
        }
        if(std::abs(xe)>=tol1 && std::abs(fa)>std::abs(fb))
        {
            double s=fb/fa;
            double p,qq;
            if(xa==xc)
            {
                p = 2.0*xm*s;
                qq = 1.0 - s;
            }
            else
            {
                qq = fa/fc;
                double r = fb/fc;
                p=s*(2.0*xm*qq*(qq-r)-(xb-xa)*(r-1.0));
                qq=(qq-1.0)*(r-1.0)*(s-1.0);
            }
            if(p>0.0) qq=-qq;
            p = std::abs(p);

            if( p < std::min(2.0*(xb-xa)*qq-std::abs(tol1*qq) , std::abs(xe*qq)/2.0) )
            {
                xe=xd;
                xd=p/qq;
            }
            else
            {
                xd=2.0*(xb-xa);
                xe=xd;
            }
        }
        else
        {
            xd = 2.0*(xb-xa);
            xe = xd;
        }

        double fab = (fa+fb)/2.0;
        double dy = (fb+(fab-fb)/(xa-xb)*xd)*xd;

        xc = xa;
        fc = fa;
        xa = xb;
        fa = fb;

        if(std::abs(xd)>tol1)
        {
            xb = xb+xd;
        }
        else
        {
            xb = xb + tol1;
        }

        xnew = xb;
        if(std::abs(dy)<conv_thr) return true;
        if(ls_step == 4) //I'm not sure if this is a good choice, but the idea is there should not be so many line search steps
                         //I feel if the line search does not converge, we'd better change the direction and restart line search
        {
            GlobalV::ofs_running << "Too many Brent steps, let's do next CG step" << std::endl;
            return true;
            //ModuleBase::WARNING_QUIT("Brent","too many steps in line search, something wrong");
        }

        return false;
    }//end bracked
    else
    {
        if(!( (xa<=xb && xb<=xc) || (xc<=xb && xb<=xa) ))
        {
            ModuleBase::WARNING_QUIT("Brent","something wrong with Brent line search!");
        }

        if((fb>0 && fc>0) || (fb<0 && fc<0))
        {
            xc = xa;
            fc = fa;
            xd = xb-xa;
            xe = xd;
        }
        if(std::abs(fc)<std::abs(fb)) //rearrange b,c so that abs(fc)>=abs(fb)
        {
            xa = xb;
            xb = xc;
            xc = xa;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        double tol1 = 2.0*e8*std::abs(xb)+0.5*e8;
        xm = 0.5*(xc-xb);
        if(std::abs(xm)<=tol1 || fb==0.0)
        {
            return true;
        }

        if(std::abs(xe)>=tol1 && std::abs(fa)>std::abs(fb))
        {
            double s=fb/fa;
            double p,
            qq;
            if(xa==xc)
            {
                p = 2.0*xm*s;
                qq = 1.0 - s;
            }
            else
            {
                qq = fa/fc;
                double r = fb/fc;
                p=s*(2.0*xm*qq*(qq-r)-(xb-xa)*(r-1.0));
                qq=(qq-1.0)*(r-1.0)*(s-1.0);
            }
            if(p>0.0) qq=-qq;
            p = std::abs(p);

            if(2.0*p < std::min( 3.0*xm*qq-std::abs(tol1*qq), std::abs(xe*qq) ) )
            {
                xe = xd;
                xd = p/qq;
            }
            else
            {
                xd = xm;
                xe = xd;
            }
        }
        else
        {
            xd = xm;
            xe = xd;
        }

        double fab = (fa+fb)/2.0;
        double dy = (fb+(fab-fb)/(xa-xb)*xd)*xd;

        xa = xb;
        fa = fb;
        ya = yb;

        if(std::abs(xc)>tol1)
        {
            xb = xb + xd;
        }
        else
        {
            if(xm > 0)
            {
                xb = xb + tol1;
            }
            else
            {
                xb = xb - tol1;
            }
        }

        xnew = xb;
        if(std::abs(dy)<conv_thr) return true;
        if(ls_step == 4)
        {
            GlobalV::ofs_running << "Too many Brent steps, let's do next CG step" << std::endl;
            return true;
            //ModuleBase::WARNING_QUIT("Brent","too many steps in line search, something wrong");
        }
        return false;
    }//end ibrack
}