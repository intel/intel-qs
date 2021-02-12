#ifndef IQS_UTILS_HPP
#define IQS_UTILS_HPP

#include <complex>

// Helpful defines, if not already provided.
#define DO_PRAGMA(x) _Pragma(#x)

#ifndef TODO
#define TODO(x) DO_PRAGMA(message("\033[30;43mTODO\033[0m - " #x))
#endif

#ifndef INFO
#define INFO(x) DO_PRAGMA(message("\033[30;46mINFO\033[0m - " #x))
#endif

#define D(x) ((double)(x))
#define UL(x) ((std::size_t)(x))
#define sec() time_in_seconds()
#define xstr(s) __str__(s)
#define __str__(s) #s

/////////////////////////////////////////////////////////////////////////////////////////

using ComplexSP = std::complex<float>;
using ComplexDP = std::complex<double>;

template<typename T>
struct extract_value_type //lets call it extract_value_type
{
    typedef T value_type;
};

template<template<typename> class X, typename T>
struct extract_value_type<X<T>>   //specialization
{
    typedef T value_type;
};

/////////////////////////////////////////////////////////////////////////////////////////

double time_in_seconds(void);

/////////////////////////////////////////////////////////////////////////////////////////

namespace qhipster {

/// Utility method to inform on the currently set compiler flags.
void WhatCompileDefinitions();

}

/////////////////////////////////////////////////////////////////////////////////////////

#endif	// header guard IQS_UTILS_HPP
