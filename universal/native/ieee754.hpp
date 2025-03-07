#pragma once
// ieee754.hpp: manipulation functions for IEEE-754 native types
//
// Copyright (C) 2017-2023 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <sstream>
#include <iomanip>
#include <cmath>    // for frexpf/frexp/frexpl  float/double/long double fraction/exponent extraction
#include <limits>
#include <tuple>

// configure the low level compiler interface to deal with floating-point bit manipulation
#include "../utility/architecture.hpp"
#include "../utility/compiler.hpp"
#include "../utility/bit_cast.hpp"
#include "../utility/long_double.hpp"

// set up the database of compiler/architecture specific floating-point parameters
#include "../native/ieee754_parameter.hpp"
#include "../native/ieee754_decoder.hpp"
#include "../native/ieee754_type_tag.hpp"

////////////////////////////////////////////////////////////////////////////////////////
// enable throwing specific exceptions long double special handling of 128bit fields
#if !defined(BITBLOCK_THROW_ARITHMETIC_EXCEPTION)
// default is to use std::cerr for signalling an error
#define BITBLOCK_THROW_ARITHMETIC_EXCEPTION 0
#endif

// if the compiler environment allows, set up
// constexpr compatible bit casts, otherwise
// fallback to nonconstexpr bit casts.
#include "extract_fields.hpp"

// functions that do not need to be constexpr
#include "nonconst_bitcast.hpp"
#include "ieee754_float.hpp"
#include "ieee754_double.hpp"
#include "ieee754_longdouble.hpp"
// above includes are a refactoring of this old include
//#include <universal/native/nonconstexpr754.hpp>

// support functions
#include "integers.hpp"
#include "ieee754_parameter_ostream.hpp"
#include "manipulators.hpp"
#include "attributes.hpp"
#include "../traits/arithmetic_traits.hpp"

// numeric helpers
#include "ieee754_numeric.hpp"
