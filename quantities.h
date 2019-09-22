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


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// It would be nice if integral types worked, too, but not so important. 
// Also there would be issues of truncation with divides and ratios.
//using base_type = std::complex<float>; // Lots of operations are not constexpr, or there are other issues.
//using base_type = int; // Doesn't work dues to rounding issues mainly.
using base_type = float;

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
using kilogram_d     = si::make_dimension<si::kilogram_u<1>>;
using meter_d        = si::make_dimension<si::meter_u<1>>;
using second_d       = si::make_dimension<si::second_u<1>>;
using ampere_d       = si::make_dimension<si::ampere_u<1>>;
using kelvin_d       = si::make_dimension<si::kelvin_u<1>>;
using mole_d         = si::make_dimension<si::mole_u<1>>;
using candela_d      = si::make_dimension<si::candela_u<1>>;

using meter_sq_d     = si::make_dimension<si::meter_u<2>>;
using meter_per_s_d  = si::dimension_divide<meter_d, second_d>;
using meter_per_ss_d = si::dimension_divide<meter_per_s_d, second_d>;

using radian_d    = si::scalar_d; // rad
using steradian_d = si::scalar_d; // sr
using hertz_d     = si::make_dimension<si::second_u<-1>>; // Hz
using newton_d    = si::make_dimension<si::kilogram_u<1>, si::meter_u<1>, si::second_u<-2>>; // N
using pascal_d    = si::make_dimension<si::kilogram_u<1>, si::meter_u<-1>, si::second_u<-2>>; // Pa
using joule_d     = si::dimension_multiply<newton_d, meter_d>; // J
using watt_d      = si::dimension_divide<joule_d, second_d>; // W
using coulomb_d   = si::make_dimension<si::ampere_u<1>, si::second_u<1>>; // C
using volt_d      = si::dimension_divide<joule_d, coulomb_d>; // V
using farad_d     = si::dimension_divide<coulomb_d, volt_d>; // F
using ohm_d       = si::dimension_divide<volt_d, ampere_d>; // Ohm
using siemens_d   = si::dimension_divide<ampere_d, volt_d>; // S
using weber_d     = si::dimension_multiply<volt_d, second_d>; // Wb
using tesla_d     = si::dimension_divide<weber_d, meter_sq_d>; // T
using henry_d     = si::dimension_divide<weber_d, ampere_d>; // H

//#define CREATE_QUANTITY(name, dimension, suffix) \
//using name = si::quantity<base_type, dimension>; \
//constexpr name operator""_##suffix(long double value) noexcept { return name(value); } \
//constexpr name operator""_##suffix(unsigned long long int value) noexcept { return name(value); }


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Mass units
// The kilogram is the base unit, but this messes with the standard scale factors milli, kilo, ...
// so we create the gram first.
using gram_t = si::quantity<base_type, kilogram_d>;
constexpr gram_t operator""_g(long double value) noexcept { return gram_t(value); }
constexpr gram_t operator""_g(unsigned long long int value) noexcept { return gram_t(value); }

using milligram_t = si::milli<gram_t>;
constexpr milligram_t operator""_mg(long double value) noexcept { return milligram_t(value); }
constexpr milligram_t operator""_mg(unsigned long long int value) noexcept { return milligram_t(value); }

using kilogram_t = si::kilo<gram_t>;
constexpr kilogram_t operator""_kg(long double value) noexcept { return kilogram_t(value); }
constexpr kilogram_t operator""_kg(unsigned long long int value) noexcept { return kilogram_t(value); }

using pound_t = si::quantity_scale<gram_t, si::make_ratio<454>>;
constexpr pound_t operator""_lb(long double value) noexcept { return pound_t(value); }
constexpr pound_t operator""_lb(unsigned long long int value) noexcept { return pound_t(value); }

// Distance units
using meter_t = si::quantity<base_type, meter_d>;
constexpr meter_t operator""_m(long double value) noexcept { return meter_t(value); }
constexpr meter_t operator""_m(unsigned long long int value) noexcept { return meter_t(value); }

using kilometer_t = si::kilo<meter_t>;
constexpr kilometer_t operator""_km(long double value) noexcept { return kilometer_t(value); }
constexpr kilometer_t operator""_km(unsigned long long int value) noexcept { return kilometer_t(value); }

using centimeter_t = si::centi<meter_t>;
constexpr centimeter_t operator""_cm(long double value) noexcept { return centimeter_t(value); }
constexpr centimeter_t operator""_cm(unsigned long long int value) noexcept { return centimeter_t(value); }

using millimeter_t = si::milli<meter_t>;
constexpr millimeter_t operator""_mm(long double value) noexcept { return millimeter_t(value); }
constexpr millimeter_t operator""_mm(unsigned long long int value) noexcept { return millimeter_t(value); }


using diameter_t = si::quantity<base_type, meter_d, std::ratio<1>, struct circle_tag>;
using radius_t   = si::quantity_scale<diameter_t, std::ratio<2>>;

static_assert(radius_t{10} == diameter_t{20});
// Fails because the tags are different.
//static_assert(radius_t{10} == meter_t{10});

// Time units
using second_t = si::quantity<base_type, second_d>;
constexpr second_t operator""_s(long double value) noexcept { return second_t(value); }
constexpr second_t operator""_s(unsigned long long int value) noexcept { return second_t(value); }

using minute_t = si::quantity_scale<second_t, si::make_ratio<60>>;
using hour_t   = si::quantity_scale<minute_t, si::make_ratio<60>>;
using day_t    = si::quantity_scale<hour_t,   si::make_ratio<24>>;

using millisecond_t = si::milli<second_t>;
constexpr millisecond_t operator""_ms(long double value) noexcept { return millisecond_t(value); }
constexpr millisecond_t operator""_ms(unsigned long long int value) noexcept { return millisecond_t(value); }

using microsecond_t = si::micro<second_t>;
constexpr microsecond_t operator""_us(long double value) noexcept { return microsecond_t(value); }
constexpr microsecond_t operator""_us(unsigned long long int value) noexcept { return microsecond_t(value); }

// Frequency units
using hertz_t = si::quantity<base_type, hertz_d>;
constexpr hertz_t operator""_Hz(long double value) noexcept { return hertz_t(value); }
constexpr hertz_t operator""_Hz(unsigned long long int value) noexcept { return hertz_t(value); }

using megahertz_t = si::mega<hertz_t>;
constexpr megahertz_t operator""_MHz(long double value) noexcept { return megahertz_t(value); }
constexpr megahertz_t operator""_MHz(unsigned long long int value) noexcept { return megahertz_t(value); }

// Velocity units
using meter_per_s_t = si::quantity_divide<meter_t, second_t>;
constexpr meter_per_s_t operator""_mps(long double value) noexcept { return meter_per_s_t(value); }
constexpr meter_per_s_t operator""_mps(unsigned long long int value) noexcept { return meter_per_s_t(value); }

// Acceleration units
using meter_per_ss_t = si::quantity_divide<meter_per_s_t, second_t>;
constexpr meter_per_ss_t operator""_mpss(long double value) noexcept { return meter_per_ss_t(value); }
constexpr meter_per_ss_t operator""_mpss(unsigned long long int value) noexcept { return meter_per_ss_t(value); }

