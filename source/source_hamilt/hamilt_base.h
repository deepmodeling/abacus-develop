#ifndef HAMILT_BASE_H
#define HAMILT_BASE_H

#include <string>

namespace hamilt
{

/**
 * @brief Base class for Hamiltonian
 *
 * This is a non-template base class for Hamilt<T, Device>.
 * It provides a common interface for all Hamiltonian types,
 * allowing ESolver to manage Hamiltonian without template parameters.
 */
class HamiltBase
{
  public:
    virtual ~HamiltBase() {}

    /**
     * @brief Update Hamiltonian for a specific k-point
     *
     * @param ik k-point index
     */
    virtual void updateHk(const int ik) { return; }

    /**
     * @brief Refresh the status of Hamiltonian
     *
     * @param yes whether to refresh
     */
    virtual void refresh(bool yes = true) { return; }

    /**
     * @brief Get the class name
     *
     * @return class name
     */
    virtual std::string get_classname() const { return "none"; }

    /**
     * @brief Get the operator chain (as void* to avoid template)
     *
     * @return pointer to operator chain
     */
    virtual void* get_ops() { return nullptr; }
};

} // namespace hamilt

#endif // HAMILT_BASE_H
