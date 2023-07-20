#include "paw_element.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"

// ptilde_l(q) = int_0^{rc} dr r ptilde_l(r) j_l(qr)
// and 2nd derivative d2ptilde_l(q) for spline
void Paw_Element::transform_ptilde()
{
    // Defining the q grid; I will be following Quantum Espresso here
    double dq = 0.01;
    int    nq = int( (std::sqrt(ecutwfc) / dq + 4) * cell_factor );

    const double pi = 3.141592653589793238462643383279502884197;
    const double twopi = 2.0 * pi;

    std::vector<double> integrand;
    integrand.resize(nr);

    qgrid.resize(nq);
    for (int iq = 0; iq < nq; iq++)
    {
        qgrid[iq] = double(iq) * dq;
    }

    ptilde_q.resize(nstates);
    d2ptilde_q.resize(nstates);
    for(int istate = 0; istate < nstates; istate ++)
    {
        ptilde_q[istate].resize(nq);
        d2ptilde_q[istate].resize(nq);
        int l = lstate[istate]; // the l quantum number

        for (int iq = 0; iq < nq; iq++)
        {
            ptilde_q[istate][iq] = this->spherical_bessel_transform(l, ptilde_r[istate], qgrid[iq]);
        }

        // end point condition for spline
        double yp1 = 0.0, ypn;
        if(l == 1)
        {
            for(int ir = 0; ir < nr; ir ++)
            {
                integrand[ir] = twopi * ptilde_r[istate][ir] * std::pow(rr[ir],3);
            }
            yp1 = this -> simpson_integration(integrand) / 3.0;
        }

        double bes, besp;
        for(int ir = 0; ir < nr; ir ++)
        {
            double x = rr[ir] * double(nq - 1) * dq;

            this -> spherical_bessel_function(l,x,bes,besp,1);
            integrand[ir] = twopi * besp * ptilde_r[istate][ir] * std::pow(rr[ir],3);
        }

        ypn = this -> simpson_integration(integrand);

        this -> spline(qgrid, ptilde_q[istate], d2ptilde_q[istate],yp1,ypn);
    }
}

// int_0^{rc} dr r fr(r) j_l(qr)
double Paw_Element::spherical_bessel_transform(const int l, std::vector<double> & fr, const double q) const
{

    if(std::abs(q) < 1e-8 && l != 0) return 0.0;

    assert(fr.size() == nr);
    std::vector<double> integrand;
    integrand.resize(nr);

    if(std::abs(q) < 1e-8)
    {
        for(int ir = 0; ir < nr; ir++)
        {
            integrand[ir] = fr[ir] * rr[ir] * rr[ir]; // r j_l(qr) ptilde(r)
        }        
    }
    else
    {
        for(int ir = 0; ir < nr; ir++)
        {
            double x = rr[ir] * q;
            double sph_bes, tmp;
            this-> spherical_bessel_function(l,x,sph_bes,tmp,0);
            integrand[ir] = sph_bes * fr[ir] * rr[ir] * rr[ir]; // r j_l(qr) ptilde(r)
        }
    }

    return this->simpson_integration(integrand);
}

void Paw_Element::spherical_bessel_function(const int l, const double xx,
    double & bes, double & besp, const bool calc_der)
{

    double tol = 0.01, prec = 1e-15;
    int imax = 40;

    switch (l)
    {
        case 0:
            if(xx < tol){
                bes = 1.0-std::pow(xx,2)/6.0*(1.0-std::pow(xx,2)/20.0);
                if(calc_der) besp= (-10.0+xx*xx)*xx/30.0;
            }
            else{
                bes = sin(xx)/xx;
                if(calc_der) besp= -(sin(xx)-xx*cos(xx))/std::pow(xx,2);
            }
            break;
        case 1:
            if(xx < tol){
                bes = (10.0-xx*xx)*xx/30.0;
                if(calc_der) besp= (10.0-3.0*xx*xx)/30.0;
            }
            else{
                bes = (sin(xx)-xx*cos(xx))/std::pow(xx,2);
                if(calc_der) besp= ((xx*xx-2.0)*sin(xx)+2.0*xx*cos(xx))/std::pow(xx,3);
            }
            break;
        case 2:
            if(xx < tol){
                bes = xx*xx/15.0-std::pow(xx,4)/210.0;
                if(calc_der) besp= (1.0-xx*xx/7.0)*2.0*xx/15.0;
            }
            else{
                bes = ((3.0-std::pow(xx,2))*sin(xx)-3.0*xx*cos(xx))/std::pow(xx,3);
                if(calc_der) besp= ((4.0*xx*xx-9.0)*sin(xx)+(9.0-xx*xx)*xx*cos(xx))/std::pow(xx,4);
            }
            break;
        case 3:
            if(xx < tol){
                bes=xx*xx*xx/105.0-std::pow(xx,5)/1890.0+std::pow(xx,7)/83160.0;
                if(calc_der) besp= (1.0/35-xx*xx/378.0+std::pow(xx,4)/11880.0)*xx*xx;
            }
            else{
                bes=(15.0*sin(xx)-15.0*xx*cos(xx)-6.0*std::pow(xx,2)*sin(xx)+std::pow(xx,3)*cos(xx))/std::pow(xx,4);
                if(calc_der) besp= ((-60.0+27.0*xx*xx-std::pow(xx,4))*sin(xx)+(60.0*xx-7.0*std::pow(xx,3))*cos(xx))/std::pow(xx,5);
            }
            break;
        default:
            if (std::abs(xx) < prec){
                bes = 0.0;
                besp = 0.0;
            }
            else
            {
                double xxinv = 1.0 / xx;
                if(xx < double(l))
                {
                    double fact = 1.0;
                    for(int il = 1; il < l+1; il ++)
                    {
                        fact=fact*xx/double(2*il+1);
                    }

                    double xx2=0.5*xx*xx;
                    double jn = 1.0, jr = 1.0;
                    int ii = 0;
                    while(std::abs(jr)>=prec && ii<imax)
                    {
                        ii=ii+1;
                        jr=-jr*xx2/double(ii*(2*(l+ii)+1));
                        jn=jn+jr;
                    }
                    bes = jn * fact;
                    if(std::abs(jr)>prec) ModuleBase::WARNING_QUIT("spherical_bessel_function","Bessel function did not converge!");

                    if(calc_der)
                    {
                        double factp=fact*xx/double(2*l+3);
                        double jnp = 1.0;
                        jr = 1.0;
                        ii = 0;
                        while(std::abs(jr)>=prec && ii<imax)
                        {
                            ii=ii+1;
                            jr=-jr*xx2/double(ii*(2*(l+ii)+3));
                            jnp=jnp+jr;
                        }
                        besp=-jnp*factp+jn*fact*xxinv*double(l);
                        if (std::abs(jr)>prec) ModuleBase::WARNING_QUIT("spherical_bessel_function","1st der. of Bessel function did not converge!");
                    }
                }
                else
                {
                    double jn =sin(xx)*xxinv;
                    double jnp=(-cos(xx)+jn)*xxinv;
                    double jr;
                    for(int il=2; il < l+2; il++)
                    {
                        jr=-jn+double(2*il-1)*jnp*xxinv;
                        jn=jnp;
                        jnp=jr;
                    }
                    bes=jn;
                    if(calc_der) besp =-jnp+jn *xxinv*double(l);
                }
            }
    }
}

// Adapted from m_pawrad/pawrad_init and m_pawrad/simp_gen
// as well as other relevant codes from ABINIT
// Only radial grid a*(exp(d*i)-1) for now (will add more if needed later)
// which corresponds to mesh_type=2 in simp_gen
// search for case("r=a*(exp(d*i)-1)") in m_pawpsp.F90, for example
double Paw_Element::simpson_integration(std::vector<double> & f) const
{
    int ir=nr;
    while (std::abs(f[ir]) < 1e-20)
    {
       ir=ir-1;
    }

    //stores factors for carrying out simpson integration
    std::vector<double> simp_fact;
    int simp_int_meshsz;

    int mmax=std::min(ir+1,nr);
    this->prepare_simpson_integration(rr[mmax], simp_int_meshsz, simp_fact);

    double val = 0.0;
    for (int i=0; i<simp_int_meshsz; i++)
    {
        val += f[i] * simp_fact[i];
    }

    return val;
}

void Paw_Element::prepare_simpson_integration(const double r_for_intg, int & meshsz, std::vector<double> & simp_fact) const
{
// generate rad_factor
    std::vector<double> rad_fact;
    rad_fact.resize(nr);

    simp_fact.resize(nr);
    std::fill(simp_fact.begin(), simp_fact.end(), 0.0);

    int isim = 2;
    double stepint = lstep;
   
    rad_fact[0] = rstep;
    for(int ir=1; ir < nr; ir++)
    {
        rad_fact[ir] = rr[ir] + rstep;
    }

// get mesh for integration
    meshsz = nr;
    if (r_for_intg > 0.0)
    {
        int nr_for_intg = int(1e-8+std::log(1.0+r_for_intg/rstep)/lstep)+1;
        int ir = std::min(nr_for_intg, nr);
        if (ir < nr)
        {
            if (std::abs(rr[ir+1]-r_for_intg) < std::abs(rr[ir]-r_for_intg)) ir=ir+1;
        }
        if (ir > 1)
        {
            if (std::abs(rr[ir-1]-r_for_intg) < std::abs(rr[ir]-r_for_intg)) ir=ir-1;
        }
        meshsz = ir;
    }

    //std::cout << "meshsz : " << meshsz << std::endl;

//get simp_fact
    double hh = stepint / 3.0;
    
    simp_fact[meshsz] = hh * rad_fact[meshsz];
    
    int ir_last = 0;
    for(int ir = meshsz; ir > isim - 1; ir-=2)
    {
        simp_fact[ir-1]=4.0*hh*rad_fact[ir-1];
        simp_fact[ir-2]=2.0*hh*rad_fact[ir-2];
        ir_last=ir-2;
    }
    
    simp_fact[ir_last] = 0.5 * simp_fact[ir_last];

    //for(int ir=0;ir<nr;ir++)
    //{
    //    std::cout << simp_fact[ir] << std::endl;
    //}

}

void Paw_Element::spline(const std::vector<double> & r, const std::vector<double> & f,
    std::vector<double> & d2f, const double yp1, const double ypn) const
{
    assert(f.size() > 1);
    assert(f.size() == d2f.size());
    assert(r.size() == f.size());

    int size = f.size();

    std::vector<double> tmp;
    
    tmp.resize(size);

    for(int i = 0; i < size - 1; i++)
    {
        if(r[i] >= r[i+1]) ModuleBase::WARNING_QUIT("Paw_Element","The knots must be strictly increasing!");
    }

    int ibcbeg=1; if(yp1>1e30) ibcbeg=0;
    int ibcend=1; if(ypn>1e30) ibcend=0;

    //Set the first and last equations
    if (ibcbeg==0){
        d2f[0] = 0.0;
        tmp[0] = 0.0;
    }
    else
    {
        d2f[0] = -0.5;
        tmp[0] = (3.0/(r[1]-r[0]))*((f[1]-f[0])/(r[1]-r[0])-yp1);
    }
    
    if (ibcend==0){
        d2f[size-1] = 0.0;
        tmp[size-1] = 0.0;
    }
    else{
        d2f[size-1] = 0.5;
        tmp[size-1] = (3.0/(r[size-1]-r[size-2]))*(ypn-(f[size-1]-f[size-2])/(r[size-1]-r[size-2]));
    }

    //Set the intermediate equations
    for(int i=1;i<size-1;i++)
    {
        double ratio = (r[i]-r[i-1])/(r[i+1]-r[i-1]);
        double pinv = 1.0/(ratio*d2f[i-1] + 2.0);
        d2f[i] = (ratio-1.0)*pinv;
        tmp[i] = (6.0*((f[i+1]-f[i])/(r[i+1]-r[i])-(f[i]-f[i-1]) 
            / (r[i]-r[i-1])) / (r[i+1]-r[i-1])-ratio*tmp[i-1])*pinv;
        if (std::abs(tmp[i])<1e-15) tmp[i]=0.0;
    }

    //Solve the equations
    d2f[size-1] = (tmp[size-1]-d2f[size-1]*tmp[size-2])/(d2f[size-1]*d2f[size-2]+1.0);
    for(int k=size-2;k>-1;k--)
    {
        d2f[k]=d2f[k]*d2f[k+1]+tmp[k];
    }

}

// Note : it is required that r is uniformly spaced
double Paw_Element::splint(const std::vector<double> & r, const std::vector<double> & f,
    const std::vector<double> & d2f, const double x) const
{

    double rmin = r[0];
    double rmax = r[r.size() - 1];

    double de = (rmax - rmin) / double(r.size() - 1);
    
    double de2_dby_six = de* de / 6.0;
    double de_dby_six = de / 6.0;

    // if new argument is beyond end points, take value at end point
    if(x <= rmin) return f[0];
    if(x >= rmax) return f[r.size() - 1];

    int jspl = std::floor((x - rmin)/de);
    double d = x - r[jspl];
    double bb = d / de;
    double aa = 1.0 - bb;
    double cc = aa*(aa*aa - 1.0) * de2_dby_six;
    double dd = bb*(bb*bb - 1.0) * de2_dby_six;
    return aa * f[jspl] + bb*f[jspl+1] + cc*d2f[jspl] + dd*d2f[jspl+1];

}