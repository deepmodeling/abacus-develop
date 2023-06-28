#include "paw_element.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"

// ptilde_l(q) = int_0^{rc} dr r ptilde_l(r) j_l(qr)
void Paw_Element::transform_ptilde()
{
    // Defining the q grid; I will be following Quantum Espresso here
    double dq = 0.01;
    int    nq = int( (std::sqrt(ecutwfc) / dq + 4) * cell_factor );

    ptilde_q.resize(nstates);
    for(int istate = 0; istate < nstates; istate ++)
    {
        ptilde_q[istate].resize(nq);
        int l = lstate[istate]; // the l quantum number

        for (int iq = 0; iq < nq; iq++)
        {
            double q = double(iq) * dq;
            ptilde_q[istate][iq] = this->spherical_bessel_transform(l, ptilde_r[istate], q);
        }
    }
}

// int_0^{rc} dr r fr(r) j_l(qr)
double Paw_Element::spherical_bessel_transform(const int l, std::vector<double> & fr, const double q) const
{

    assert(fr.size() == nr);
    std::vector<double> integrand;
    integrand.resize(nr);

    for(int ir = 0; ir < nr; ir++)
    {
        double x = rr[ir] * q;
        double sph_bes = spherical_bessel_function(l,x);
        integrand[ir] = sph_bes * fr[ir] * rr[ir]; // r j_l(qr) ptilde(r)
    }

    return this->simpson_integration(integrand);
}

//Adapted from m_special_funcs/paw_jbessel_4spline and m_paw_numeric/paw_jbessel of ABINIT
double Paw_Element::spherical_bessel_function(const int l, const double xx)
{

    double bes;
    double tol = 0.001, prec = 1e-56;
    int imax = 40;

    switch (l)
    {
        case 0:
            if(xx < tol) bes = 1.0-std::pow(xx,2)/6.0*(1.0-std::pow(xx,2)/20.0);
            else bes = sin(xx)/xx;
            break;
        case 1:
            if(xx < tol) bes = (10.0-xx*xx)*xx/30.0;
            else bes=(sin(xx)-xx*cos(xx))/std::pow(xx,2);
            break;
        case 2:
            if(xx < tol) bes=xx*xx/15.0-std::pow(xx,4)/210.0;
            else bes=((3.0-std::pow(xx,2))*sin(xx)-3.0*xx*cos(xx))/std::pow(xx,3);
            break;
        case 3:
            if(xx < tol) bes=xx*xx*xx/105.0-std::pow(xx,5)/1890.0+std::pow(xx,7)/83160.0;
            else bes=(15.0*sin(xx)-15.0*xx*cos(xx)-6.0*std::pow(xx,2)*sin(xx)+std::pow(xx,3)*cos(xx))/std::pow(xx,4);
            break;
        default:
            if (std::abs(xx) < prec) bes = 0.0;
            else
            {
                double xxinv = 1.0 / xx;
                if(xx < 1.0)
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
                }
            }
    }

    return bes;

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