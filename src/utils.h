#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#ifdef _MSC_VER
#include "xmmintrin.h"
#include <nmmintrin.h>
#include <immintrin.h>
#endif

template <typename T>
void _my_prefetch(T ptr)
{
#ifdef _WIN32
	_mm_prefetch((const char*) ptr, _MM_HINT_T0);
#else
	__builtin_prefetch((const char*) ptr);
#endif
}

#if defined(_MSC_VER)  /* Visual Studio */
#define REFRESH_FORCE_INLINE __forceinline
#define REFRESH_NO_INLINE __declspec(noinline)
#elif defined(__GNUC__)
#define REFRESH_FORCE_INLINE __inline__ __attribute__((always_inline, unused))
#define REFRESH_NO_INLINE __attribute__((noinline))
#else
#define REFRESH_FORCE_INLINE
#define REFRESH_NO_INLINE
#endif
