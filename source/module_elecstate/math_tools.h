#include "module_base/blas_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_psi/psi.h"

#ifdef __MPI
void psiMulPsiMpi(const psi::Psi<double>& psi1,
                  const psi::Psi<double>& psi2,
                  psi::Psi<double>& dm_out,
                  const int* desc_psi,
                  const int* desc_dm)
{
    const double one_float = 1.0, zero_float = 0.0;
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = desc_dm[2];
    const int nbands = desc_psi[3];
    pdgemm_(&N_char,
            &T_char,
            &nlocal,
            &nlocal,
            &nbands,
            &one_float,
            psi1.get_pointer(),
            &one_int,
            &one_int,
            desc_psi,
            psi2.get_pointer(),
            &one_int,
            &one_int,
            desc_psi,
            &zero_float,
            dm_out.get_pointer(),
            &one_int,
            &one_int,
            desc_dm);
}

void psiMulPsiMpi(const psi::Psi<std::complex<double>>& psi1,
                  const psi::Psi<std::complex<double>>& psi2,
                  psi::Psi<std::complex<double>>& dm_out,
                  const int* desc_psi,
                  const int* desc_dm)
{
    const complex<double> one_complex = {1.0, 0.0}, zero_complex = {0.0, 0.0};
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = desc_dm[2];
    const int nbands = desc_psi[3];
    pzgemm_(&N_char,
            &T_char,
            &nlocal,
            &nlocal,
            &nbands,
            (const double*)(&one_complex),
            psi1.get_pointer(),
            &one_int,
            &one_int,
            desc_psi,
            psi2.get_pointer(),
            &one_int,
            &one_int,
            desc_psi,
            (const double*)(&zero_complex),
            dm_out.get_pointer(),
            &one_int,
            &one_int,
            desc_dm);
}

#else
void psiMulPsi(psi::Psi<double>& psi1, psi::Psi<double>& psi2, psi::Psi<double>& dm_out)
{
    const double one_float = 1.0, zero_float = 0.0;
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = psi1.get_nbasis();
    const int nbands = psi1.get_nbands();
    dgemm_(&N_char,
           &T_char,
           &nlocal,
           &nlocal,
           &nbands,
           &one_float,
           psi1.get_pointer(),
           &nlocal,
           psi2.get_pointer(),
           &nlocal,
           &zero_float,
           dm_out.get_pointer(),
           &nlocal);
}

void psiMulPsi(psi::Psi<double>& psi1, psi::Psi<double>& psi2, psi::Psi<double>& dm_out)
{
    const int one_int = 1;
    const char N_char = 'N', T_char = 'T';
    const int nlocal = psi1.get_nbasis();
    const int nbands = psi1.get_nbands();
    const complex<double> one_complex = {1.0, 0.0}, zero_complex = {0.0, 0.0};
    zgemm_(&N_char,
           &T_char,
           &nlocal,
           &nlocal,
           &nbands,
           &one_complex,
           psi1.get_pointer(),
           &nlocal,
           psi2.get_pointer(),
           &nlocal,
           &zero_complex,
           dm_out.get_pointer(),
           &nlocal);
}

#endif