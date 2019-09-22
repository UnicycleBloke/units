/////////////////////////////////////////////////////////////////////////////////////////////
// The MIT License (MIT)
//
// Copyright (c) 2019 Alan Chambers (unicycle.bloke@gmail.com)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
/////////////////////////////////////////////////////////////////////////////////////////////
#pragma once
#include "units.h"
#include "quantities.h"


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Compile time tests
template <typename T1, typename T2, typename = decltype(std::declval<T1>() + std::declval<T2>())>
constexpr bool has_add(T1 a, T2 b) { return true; }
constexpr bool has_add(...)        { return false; }

template <typename T1, typename T2, typename = decltype(std::declval<T1>() - std::declval<T2>())>
constexpr bool has_sub(T1 a, T2 b) { return true; }
constexpr bool has_sub(...)        { return false; }

template <typename T1, typename T2, typename = decltype(std::declval<T1>() * std::declval<T2>())>
constexpr bool has_mul(T1 a, T2 b) { return true; }
constexpr bool has_mul(...)        { return false; }

template <typename T1, typename T2, typename = decltype(std::declval<T1>() / std::declval<T2>())>
constexpr bool has_div(T1 a, T2 b) { return true; }
constexpr bool has_div(...)        { return false; }

template <typename T, typename = decltype(si::exp(std::declval<T>()))>
constexpr bool has_exp(T)   { return true; }
constexpr bool has_exp(...) { return false; }


// Equivalence of construction methods.
static_assert(meter_t{100} == 100_m);
static_assert(100_m == meter_t{100});
static_assert(meter_t{100} == 100.0_m);
static_assert(100.0_m == meter_t{100});

static_assert(meter_t(100) == 100_m);
static_assert(100_m == meter_t(100));
static_assert(meter_t(100) == 100.0_m);
static_assert(100.0_m == meter_t(100));

// Combining dimension
static_assert(200_mps == 4000_m / 20_s);
static_assert(200_mps == 4_km / 20000_ms);
static_assert(1000 / 50_s == 20_Hz);
static_assert(100 / 50_s == 2_Hz);
static_assert(1000 / 500_s == 2_Hz);
static_assert(1000 / 50_ms == 20000_Hz);
static_assert(1000 / 50_us == 20_MHz);
static_assert(60 * 1000_ms == minute_t{1});
static_assert(1000_ms * 60 == minute_t{1});
static_assert(minute_t{1} == 60_ms * 1000);

static_assert(1_kg == 50 * 20_g);
static_assert(1_kg == 50_g * 20);
static_assert(1_kg == 5000_g / 5);
static_assert(50_Hz == 1 / 20_ms);

// Comparison
static_assert(10_m == 1000_cm);
static_assert(100_cm == 1000_mm);
static_assert(10_m > 999_cm);
static_assert(100_cm > 999_mm);
static_assert(10_m >= 1000_cm);
static_assert(100_cm >= 1000_mm);
static_assert(10_m <= 1000_cm);
static_assert(100_cm <= 1000_mm);
static_assert(10_m < 1001_cm);
static_assert(100_cm < 1001_mm);
static_assert(10_m != 1001_cm);
static_assert(100_cm != 1001_mm);

// Conversion to bool
// Explicit
static_assert(bool(meter_t{100}) == true);
static_assert(bool(meter_t{0}) == false);
static_assert(false != bool(meter_t{100}));
static_assert(true != bool(meter_t{0}));
// Double NOT
static_assert(!!meter_t{100} == true);
static_assert(!!meter_t{0} == false);
static_assert(false != !!meter_t{100});
static_assert(true != !!meter_t{0});
// NOT
static_assert(!meter_t{100} == false);
static_assert(!meter_t{0} == true);
static_assert(!100_m == false);
static_assert(!0_m == true);

// Conversion of ratios
static_assert(meter_t{100} == 100_m);
static_assert(kilometer_t{10} == 10000_m);
static_assert(kilometer_t{10} == 10000000_mm);
static_assert(12345_mm == 1234.5_cm);
static_assert(12_m == 12000_mm);

constexpr meter_t m1 = meter_t(1.23);
constexpr centimeter_t cm1 = m1;
static_assert(m1.value() == static_cast<base_type>(1.23)); 
static_assert(cm1.value() == 123);

// Copy construct type with different ratio 
static_assert(centimeter_t(m1) == 123_cm);

// Multiplication and division, including with scalars 
static_assert(12_mps == 48_m / 4_s);
static_assert(12_mps * 4_s == 48_m);
static_assert(12_mps * 4 == 48_mps); 
static_assert(4 * 12_mps == 48_mps); 
static_assert(12_mps / 4 == 3_mps); 
static_assert(12 / 4_s == 3_Hz); 

// Addition and subtraction
static_assert(12_s + 3_ms == 12.003_s);
static_assert(12_us + 3_ms == 3.012_ms);
static_assert(12_s - 3_ms == 11.997_s);
static_assert(1_ms - 3_ms == -2.0_ms);
static_assert(3_ms - 1_s == -0.997_s);
static_assert(3_ms - 1_s == -997_ms);

static_assert(has_add(1_m, 1_m));
static_assert(!has_add(1_m, 1_s));
static_assert(has_add(1 / 1_Hz , 1_s));
static_assert(!has_exp(1_m));
static_assert(has_exp(1_m / 1_m));

// abs(const quantity<T, D, R, Tag>& q) noexcept
static_assert(si::abs(-123_m) == 123_m);
static_assert(si::abs(-123_m) != -123_m);
// fmod(const quantity<T, D, R, Tag>& q) noexcept
// remainder(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
// remquo(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2, int* quo) noexcept 
// fmax(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
static_assert(si::fmax(12_m, 13_m) == 13_m);
static_assert(si::fmax(-12_m, -13_m) == -12_m);
// fmin(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
static_assert(si::fmin(12_m, 13_m) == 12_m);
static_assert(si::fmin(-12_m, -13_m) == -13_m);
// fdim(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
static_assert(si::fdim(12_m, 13_m) == 0_m);
static_assert(si::fmin(13_m, 12_m) == 1_m);
static_assert(std::abs(si::fdim(13_m, 12_m).value() - std::fdim(13, 12)) < std::numeric_limits<base_type>::epsilon());

// // Exponential functions - can take only scalar quantities, and the ratio is normalised.
// exp(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::exp(3_m / 1_m).value() - std::exp(3)) < std::numeric_limits<base_type>::epsilon());
// exp2(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::exp2(3_m / 1_m).value() - std::exp2(3)) < std::numeric_limits<base_type>::epsilon());
// expm1(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::expm1(3_m / 1_m).value() - std::expm1(3)) < std::numeric_limits<base_type>::epsilon());
// log(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::log(3_m / 1_m).value() - std::log(3)) < std::numeric_limits<base_type>::epsilon());
// log10(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::log10(3_m / 1_m).value() - std::log10(3)) < std::numeric_limits<base_type>::epsilon());
// log2(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::log2(3_m / 1_m).value() - std::log2(3)) < std::numeric_limits<base_type>::epsilon());
// log1p(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::log1p(3_m / 1_m).value() - std::log1p(3)) < std::numeric_limits<base_type>::epsilon());

// // Power functions - TODO pow, sqrt and cbrt
// hypot(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 

// // Trigonometric functions - can take only scalar quantities, and the ratio is normalised.
// sin(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::sin(3_m / 1_m).value() - std::sin(3)) < std::numeric_limits<base_type>::epsilon());
// cos(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::cos(3_m / 1_m).value() - std::cos(3)) < std::numeric_limits<base_type>::epsilon());
// tan(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::tan(3_m / 1_m).value() - std::tan(3)) < std::numeric_limits<base_type>::epsilon());
// asin(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::asin(1_m / 3_m).value() - std::asin(1./3)) < std::numeric_limits<base_type>::epsilon());
// acos(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::acos(1_m / 3_m).value() - std::acos(1./3)) < std::numeric_limits<base_type>::epsilon());
// atan(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::atan(1_m / 3_m).value() - std::atan(1./3)) < std::numeric_limits<base_type>::epsilon());
// atan2(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
static_assert(std::abs(si::atan2(1_m, 3_m).value() - std::atan(1./3)) < std::numeric_limits<base_type>::epsilon());

// // Hyperbolic functions - can take only scalar quantities, and the ratio is normalised.
// sinh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::sinh(3_m / 1_m).value() - std::sinh(3)) < 2*std::numeric_limits<base_type>::epsilon());
// cosh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::cosh(3_m / 1_m).value() - std::cosh(3)) < 3*std::numeric_limits<base_type>::epsilon());
// tanh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::tanh(3_m / 1_m).value() - std::tanh(3)) < std::numeric_limits<base_type>::epsilon());
// asinh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::asinh(1_m / 3_m).value() - std::asinh(1./3)) < std::numeric_limits<base_type>::epsilon());
// acosh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::acosh(9_m / 3_m).value() - std::acosh(3)) < std::numeric_limits<base_type>::epsilon());
// atanh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::atanh(1_m / 3_m).value() - std::atanh(1./3)) < std::numeric_limits<base_type>::epsilon());

// //Error and gamma functions 
// erf(const quantity<T, scalar_d, R, Tag>& q) noexcept
// erfc(const quantity<T, scalar_d, R, Tag>& q) noexcept
// tgamma(const quantity<T, scalar_d, R, Tag>& q) noexcept
// lgamma(const quantity<T, scalar_d, R, Tag>& q) noexcept

// // Nearest integer floating point operations 
// ceil(const quantity<T, D, R, Tag>& q) noexcept 
static_assert(si::ceil(123.1_m) == 124_m);
static_assert(si::ceil(-123.1_m) == -123_m);
// floor(const quantity<T, D, R, Tag>& q) noexcept 
static_assert(si::floor(123.9_m) == 123_m);
static_assert(si::floor(-123.9_m) == -124_m);
// trunc(const quantity<T, D, R, Tag>& q) noexcept 
static_assert(si::trunc(123.1_m) == 123_m);
static_assert(si::trunc(-123.1_m) == -123_m);
static_assert(si::trunc(123.9_m) == 123_m);
static_assert(si::trunc(-123.9_m) == -123_m);
// round(const quantity<T, D, R, Tag>& q) noexcept 
static_assert(si::round(123.1_m) == 123_m);
static_assert(si::round(-123.1_m) == -123_m);
static_assert(si::round(123.9_m) == 124_m);
static_assert(si::round(-123.9_m) == -124_m);


base_type test_double(base_type m, base_type s)
{ 
    return m * m * m / s; 
}


auto test_units(meter_t m, second_t s)
{
    return m * m * m / s;
}
