#ifndef OCCUPY_H
#define OCCUPY_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/vector3.h"

class Occupy
{

public:

    Occupy();
    ~Occupy();

    static void decision(const std::string &name,const std::string &smearing_method,const double &smearing_sigma);

    static const bool& gauss(void) 
	{
        return use_gaussian_broadening;
    }

    static const bool& fix(void) 
	{
        return fixed_occupations;
    }

public:

    // gaussian_broadening
    static bool use_gaussian_broadening;
    static int gaussian_type;
    static double gaussian_parameter;

    // fixed occupations
    static bool fixed_occupations;

    static void iweights(const int nks,
                         const std::vector<double>& wk,
                         const int nband,
                         const double& nelec,
                         const ModuleBase::matrix& ekb,
                         double& ef,
                         ModuleBase::matrix& wg,
                         const int& is,
                         const std::vector<int>& isk);

    static void gweights(const int nks,
                         const std::vector<double>& wk,
                         const int nband,
                         const double& nelec,
                         const double& smearing_sigma,
                         const int ngauss,
                         const ModuleBase::matrix& ekb,
                         double& ef,
                         double& demet,
                         ModuleBase::matrix& wg,
                         const int& is,
                         const std::vector<int>& isk);

    static void tweights(const int nks,const int nspin,const int nband,
						const double &nelec,const int ntetra,
                         const ModuleBase::matrix &tetra,double **ekb,double &ef,ModuleBase::matrix &wg);

    static double wsweight(const ModuleBase::Vector3<double> &r, ModuleBase::Vector3<double> *rws,const int nrws);

private:
  static void efermig(const ModuleBase::matrix& ekb,
                      const int nbnd,
                      const int nks,
                      const double& nelec,
                      const std::vector<double>& wk,
                      const double& smearing_sigma,
                      const int ngauss,
                      double& ef,
                      const int& is,
                      const std::vector<int>& isk);

  static double sumkg(const ModuleBase::matrix& ekb,
                      const int nband,
                      const int nks,
                      const std::vector<double>& wk,
                      const double& smearing_sigma,
                      const int ngauss,
                      const double& e,
                      const int& is,
                      const std::vector<int>& isk);

  static double wgauss(const double& x, const int n);

  static double w1gauss(const double& x, const int n);

  //============================
  // Needed in tweights
  //============================
  static void efermit(double** ekb,
                      const int nband,
                      const int nks,
                      const double& nelec,
                      const int nspin,
                      const int ntetra,
                      const ModuleBase::matrix& tetra,
                      double& ef);

  static double sumkt(double** ekb,
                      const int nband,
                      const int nks,
                      const int nspin,
                      const int ntetra,
                      const ModuleBase::matrix& tetra,
                      const double& eup);

  static void piksort(const int n, double* a);
};

#endif
