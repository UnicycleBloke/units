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
#include "units.h"
#include "quantities.h"


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Compile time tests

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

// Make sure that dividing like quantities creates a scalar, but no of the underlying representation.
// We rely on this for some of the tests.
static_assert(std::is_same_v<std::remove_reference_t<decltype(1_m / 1_m)>, scalar_t>);
static_assert(std::is_same_v<std::remove_reference_t<decltype(1_m *= 1)>, decltype(1_m)>);
static_assert(std::is_same_v<std::remove_reference_t<decltype(1_m /= 1)>, decltype(1_m)>);
static_assert(std::is_same_v<std::remove_reference_t<decltype(1_m += 1_m)>, decltype(1_m)>);
static_assert(std::is_same_v<std::remove_reference_t<decltype(1_m -= 1_m)>, decltype(1_m)>);

#define UNITS_EPSILON std::numeric_limits<base_type>::epsilon()
static_assert((1_mm / 1_m).value() == 1);
static_assert(std::abs((1_mm / 1_m).scaled_value() - 0.001) < UNITS_EPSILON);
static_assert(std::is_same_v<decltype(1_mm / 1_m)::ratio, std::milli>);

// We can only add-assign or add values with the same dimension.
template <typename T1, typename T2, typename = decltype(std::declval<T1>() += std::declval<T2>())>
constexpr bool has_add_assign(T1 a, T2 b) { return true; }
constexpr bool has_add_assign(...)        { return false; }

static_assert(has_add_assign(1_m, 1_m) == true);
static_assert(has_add_assign(1_m, 1_s) == false);
static_assert(has_add_assign(1_m, 1) == false);
static_assert(has_add_assign(1, 1_m) == false);

template <typename T1, typename T2, typename = decltype(std::declval<T1>() + std::declval<T2>())>
constexpr bool has_add(T1 a, T2 b) { return true; }
constexpr bool has_add(...)        { return false; }

static_assert(has_add(1_m, 1_m) == true);
static_assert(has_add(1_m, 1_s) == false);
static_assert(has_add(1_m, 1) == false);
static_assert(has_add(1, 1_m) == false);

// We can only sub-assign or sub values with the same dimension.
template <typename T1, typename T2, typename = decltype(std::declval<T1>() -= std::declval<T2>())>
constexpr bool has_sub_assign(T1 a, T2 b) { return true; }
constexpr bool has_sub_assign(...)        { return false; }

static_assert(has_sub_assign(1_m, 1_m) == true);
static_assert(has_sub_assign(1_m, 1_s) == false);
static_assert(has_sub_assign(1_m, 1) == false);
static_assert(has_sub_assign(1, 1_m) == false);

template <typename T1, typename T2, typename = decltype(std::declval<T1>() - std::declval<T2>())>
constexpr bool has_sub(T1 a, T2 b) { return true; }
constexpr bool has_sub(...)        { return false; }

static_assert(has_sub(1_m, 1_m) == true);
static_assert(has_sub(1_m, 1_s) == false);
static_assert(has_sub(1_m, 1) == false);
static_assert(has_sub(1, 1_m) == false);

// We can only div-assign scalars. Any types are valid for non-assigning div.
template <typename T1, typename T2, typename = decltype(std::declval<T1>() /= std::declval<T2>())>
constexpr bool has_div_assign(T1 a, T2 b) { return true; }
constexpr bool has_div_assign(...)        { return false; }

static_assert(has_div_assign(1_m, 1_m) == false);
static_assert(has_div_assign(1_m, 1_s) == false);
static_assert(has_div_assign(1_m, 1_m / 1_m) == true);
static_assert(has_div_assign(1_m, 1) == true);
static_assert(has_div_assign(1, 1_m) == false);

template <typename T1, typename T2, typename = decltype(std::declval<T1>() / std::declval<T2>())>
constexpr bool has_div(T1 a, T2 b) { return true; }
constexpr bool has_div(...)        { return false; }

static_assert(has_div(1_m, 1_m) == true);
static_assert(has_div(1_m, 1_s) == true);
static_assert(has_div(1_m, 1) == true);
static_assert(has_div(1, 1_m) == true);

// We can only mul-assign scalars. Any types are valid for non-assigning mul.
template <typename T1, typename T2, typename = decltype(std::declval<T1>() *= std::declval<T2>())>
constexpr bool has_mul_assign(T1 a, T2 b) { return true; }
constexpr bool has_mul_assign(...)        { return false; }

static_assert(has_mul_assign(1_m, 1_m) == false);
static_assert(has_mul_assign(1_m, 1_s) == false);
static_assert(has_mul_assign(1_m, 1_m / 1_m) == true);
static_assert(has_mul_assign(1_m, 1) == true);
static_assert(has_mul_assign(1, 1_m) == false);

template <typename T1, typename T2, typename = decltype(std::declval<T1>() * std::declval<T2>())>
constexpr bool has_mul(T1 a, T2 b) { return true; }
constexpr bool has_mul(...)        { return false; }

static_assert(has_mul(1_m, 1_m) == true);
static_assert(has_mul(1_m, 1_s) == true);
static_assert(has_mul(1_m, 1) == true);
static_assert(has_mul(1, 1_m) == true);

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
static_assert(si::fdim(13_m, 12_m) == 1_m);
static_assert(std::abs(si::fdim(13_m, 12_m).value() - std::fdim(13, 12)) < UNITS_EPSILON);

// // Exponential functions - can take only scalar quantities, and the ratio is normalised.
// exp(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::exp(3_m / 1_m) - std::exp(3)) < UNITS_EPSILON);
// exp2(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::exp2(3_m / 1_m) - std::exp2(3)) < UNITS_EPSILON);
// expm1(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::expm1(3_m / 1_m) - std::expm1(3)) < UNITS_EPSILON);
// log(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::log(3_m / 1_m) - std::log(3)) < UNITS_EPSILON);
// log10(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::log10(3_m / 1_m) - std::log10(3)) < UNITS_EPSILON);
// log2(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::log2(3_m / 1_m) - std::log2(3)) < UNITS_EPSILON);
// log1p(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::log1p(3_m / 1_m) - std::log1p(3)) < UNITS_EPSILON);

// // Power functions - TODO pow, sqrt and cbrt
// hypot(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 

// // Trigonometric functions - can take only scalar quantities, and the ratio is normalised.
// sin(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::sin(3_m / 1_m) - std::sin(3)) < UNITS_EPSILON);
// cos(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::cos(3_m / 1_m) - std::cos(3)) < UNITS_EPSILON);
// tan(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::tan(3_m / 1_m) - std::tan(3)) < UNITS_EPSILON);
// asin(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::asin(1_m / 3_m) - std::asin(1./3)) < UNITS_EPSILON);
// acos(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::acos(1_m / 3_m) - std::acos(1./3)) < UNITS_EPSILON);
// atan(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::atan(1_m / 3_m) - std::atan(1./3)) < UNITS_EPSILON);
// atan2(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
static_assert(std::abs(si::atan2(1_m, 3_m) - std::atan(1./3)) < UNITS_EPSILON);

// // Hyperbolic functions - can take only scalar quantities, and the ratio is normalised.
// sinh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::sinh(3_m / 1_m) - std::sinh(3)) < 2*UNITS_EPSILON);
// cosh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::cosh(3_m / 1_m) - std::cosh(3)) < 3*UNITS_EPSILON);
// tanh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::tanh(3_m / 1_m) - std::tanh(3)) < UNITS_EPSILON);
// asinh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::asinh(1_m / 3_m) - std::asinh(1./3)) < UNITS_EPSILON);
// acosh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::acosh(9_m / 3_m) - std::acosh(3)) < UNITS_EPSILON);
// atanh(const quantity<T, scalar_d, R, Tag>& q) noexcept
static_assert(std::abs(si::atanh(1_m / 3_m) - std::atanh(1./3)) < UNITS_EPSILON);

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

// Tagged quantities
using diameter_t = si::quantity<base_type, meter_d, std::ratio<1>, struct circle_tag>;
using radius_t   = si::quantity_scale<diameter_t, std::ratio<2>>;

static_assert(radius_t{10} == diameter_t{20});
static_assert(has_add(radius_t{1}, diameter_t{1}) == true);
static_assert(has_add(radius_t{1}, meter_t{1}) == false);
static_assert(has_add(radius_t{1}.convert_tag<void>(), meter_t{1}) == true);


int main()
{
    auto add = 12_m;
    assert((add += 3_m) == 15_m);
    assert((add += 3_cm) == 15.03_m);

    auto sub = 12_m;
    assert((sub -= 3_m) == 9_m);
    assert((sub -= 3_cm) == 8.97_m);

    auto div = 12_m;
    assert((div /= 3) == 4_m);

    auto mul = 12_m;
    assert((mul *= 3) == 36_m);
}