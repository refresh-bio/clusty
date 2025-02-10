// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2025, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#include "conversion.h"

char Conversions::digits[100000 * 5];
uint64_t Conversions::powers10[15];
double Conversions::neg_powers10[15];
Conversions::_si Conversions::_init;