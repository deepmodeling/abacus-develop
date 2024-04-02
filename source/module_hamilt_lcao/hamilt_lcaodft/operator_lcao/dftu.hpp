namespace hamilt
{

#ifndef __DFTUTEMPLATE
#define __DFTUTEMPLATE

/// The DFTU class template inherits from class T
/// it is used to calculate the non-local pseudopotential of wavefunction basis
/// Template parameters:
/// - T: base class, it would be OperatorLCAO<TK, TR> or OperatorPW<TK>
template <class T>
class DFTU : public T
{
};

#endif

}