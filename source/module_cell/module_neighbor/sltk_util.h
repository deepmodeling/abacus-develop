#ifndef UTIL_H
#define UTIL_H

#include <string>

/*** Data ***/

static struct { template<typename T> operator T*() const { return static_cast<T*>(0); } } NullPtr;

/*** Function ***/
template<typename exception>
static inline void affirm(const bool b)

{ if (!b) throw exception(); }

template<typename exception>
static inline void affirm(const bool b, const char* const message)

{ if (!b) throw exception(message); }

template<typename exception>
static inline void affirm(const bool b, const std::string& message)

{ if (!b) throw exception(message); }

#endif
