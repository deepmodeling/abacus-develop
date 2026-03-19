#ifndef PSI_PREPARE_BASE_H
#define PSI_PREPARE_BASE_H

namespace psi
{

/**
 * @brief Base class for PSIPrepare without template parameters.
 *
 * This class provides a non-template base class for PSIPrepare<T, Device>,
 * allowing Setup_Psi_pw to store a base class pointer instead of a template pointer.
 * This is part of the gradual refactoring to remove template parameters from Setup_Psi_pw.
 */
class PSIPrepareBase
{
  public:
    PSIPrepareBase() = default;
    virtual ~PSIPrepareBase() = default;
};

} // namespace psi

#endif
