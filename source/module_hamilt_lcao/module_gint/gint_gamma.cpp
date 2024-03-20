#include "gint_gamma.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

Gint_Gamma::Gint_Gamma() {

  sender_index_size = 1;
  sender_local_index = nullptr;
  sender_size_process = nullptr;
  sender_displacement_process = nullptr;
  sender_size = 1;
  sender_buffer = nullptr;

  receiver_index_size = 1;
  receiver_global_index = nullptr;
  receiver_size_process = nullptr;
  receiver_displacement_process = nullptr;
  receiver_size = 1;
  receiver_buffer = nullptr;
}

Gint_Gamma::~Gint_Gamma() {

  delete[] sender_local_index;
  delete[] sender_size_process;
  delete[] sender_displacement_process;
  delete[] sender_buffer;

  delete[] receiver_global_index;
  delete[] receiver_size_process;
  delete[] receiver_displacement_process;
  delete[] receiver_buffer;
}