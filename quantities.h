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

using radian_d       = si::scalar_d; // rad
using steradian_d    = si::scalar_d; // sr
using hertz_d        = si::make_dimension<si::second_u<-1>>; // Hz
using newton_d       = si::make_dimension<si::kilogram_u<1>, si::meter_u<1>, si::second_u<-2>>; // N
using pascal_d       = si::make_dimension<si::kilogram_u<1>, si::meter_u<-1>, si::second_u<-2>>; // Pa
using joule_d        = si::dimension_multiply<newton_d, meter_d>; // J
using watt_d         = si::dimension_divide<joule_d, second_d>; // W
using coulomb_d      = si::make_dimension<si::ampere_u<1>, si::second_u<1>>; // C
using volt_d         = si::dimension_divide<joule_d, coulomb_d>; // V
using farad_d        = si::dimension_divide<coulomb_d, volt_d>; // F
using ohm_d          = si::dimension_divide<volt_d, ampere_d>; // Ohm
using siemens_d      = si::dimension_divide<ampere_d, volt_d>; // S
using weber_d        = si::dimension_multiply<volt_d, second_d>; // Wb
using tesla_d        = si::dimension_divide<weber_d, meter_sq_d>; // T
using henry_d        = si::dimension_divide<weber_d, ampere_d>; // H


#define CREATE_QUANTITY(name, dimension, scale, suffix) \
using name = si::quantity<base_type, dimension, scale>; 

#define CREATE_SUFFIX(name, suffix) \
constexpr name operator""_##suffix(long double value) noexcept { return name(value); } \
constexpr name operator""_##suffix(unsigned long long int value) noexcept { return name(value); }

#define CREATE_QUANTITY_AND_SUFFIX(name, dimension, scale, suffix) \
CREATE_QUANTITY(name, dimension, scale, suffix) \
CREATE_SUFFIX(name, suffix)

// The kilogram is the base unit, but this messes with the standard scale 
// factors milli, kilo, ... so we create the gram first.
CREATE_QUANTITY_AND_SUFFIX(gram_t,      kilogram_d,  std::ratio<1>, g)   
CREATE_QUANTITY_AND_SUFFIX(meter_t,     meter_d,     std::ratio<1>, m)   
CREATE_QUANTITY_AND_SUFFIX(second_t,    second_d,    std::ratio<1>, s)   
CREATE_QUANTITY_AND_SUFFIX(ampere_t,    ampere_d,    std::ratio<1>, A)   
CREATE_QUANTITY_AND_SUFFIX(kelvin_t,    kelvin_d,    std::ratio<1>, K)   
CREATE_QUANTITY_AND_SUFFIX(mole_t,      mole_d,      std::ratio<1>, mol)   
CREATE_QUANTITY_AND_SUFFIX(candela_t,   candela_d,   std::ratio<1>, cd)   

CREATE_QUANTITY_AND_SUFFIX(radian_t,    radian_d,    std::ratio<1>, rad)   
CREATE_QUANTITY_AND_SUFFIX(steradian_t, steradian_d, std::ratio<1>, sr)   
CREATE_QUANTITY_AND_SUFFIX(hertz_t,     hertz_d,     std::ratio<1>, Hz)   
CREATE_QUANTITY_AND_SUFFIX(newton_t,    newton_d,    std::ratio<1>, N)   
CREATE_QUANTITY_AND_SUFFIX(pascal_t,    pascal_d,    std::ratio<1>, Pa)   
CREATE_QUANTITY_AND_SUFFIX(joule_t,     joule_d,     std::ratio<1>, J)   
CREATE_QUANTITY_AND_SUFFIX(watt_t,      watt_d,      std::ratio<1>, W)   
CREATE_QUANTITY_AND_SUFFIX(coulomb_t,   coulomb_d,   std::ratio<1>, C)   
CREATE_QUANTITY_AND_SUFFIX(volt_t,      volt_d,      std::ratio<1>, V)   
CREATE_QUANTITY_AND_SUFFIX(farad_t,     farad_d,     std::ratio<1>, F)   
CREATE_QUANTITY_AND_SUFFIX(ohm_t,       ohm_d,       std::ratio<1>, Ohm)   
CREATE_QUANTITY_AND_SUFFIX(siemens_t,   siemens_d,   std::ratio<1>, S)   
CREATE_QUANTITY_AND_SUFFIX(weber_t,     weber_d,     std::ratio<1>, Wb)   
CREATE_QUANTITY_AND_SUFFIX(tesla_t,     tesla_d,     std::ratio<1>, T)   
CREATE_QUANTITY_AND_SUFFIX(henry_t,     henry_d,     std::ratio<1>, H)   

// Just for testing really. Maybe not actually required at all.
using scalar_t = si::quantity<base_type, si::scalar_d>;
constexpr scalar_t operator""_scalar(long double value) noexcept { return scalar_t(value); }
constexpr scalar_t operator""_scalar(unsigned long long int value) noexcept { return scalar_t(value); }


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Mass units
using milligram_t = si::milli<gram_t>;
CREATE_SUFFIX(milligram_t, mg)

using kilogram_t = si::kilo<gram_t>;
CREATE_SUFFIX(kilogram_t, kg)

using pound_t = si::quantity_scale<gram_t, si::make_ratio<454>>;
CREATE_SUFFIX(pound_t, lb)

// Distance units
using kilometer_t = si::kilo<meter_t>;
CREATE_SUFFIX(kilometer_t, km)

using centimeter_t = si::centi<meter_t>;
CREATE_SUFFIX(centimeter_t, cm)

using millimeter_t = si::milli<meter_t>;
CREATE_SUFFIX(millimeter_t, mm)

// Time units
using minute_t = si::quantity_scale<second_t, si::make_ratio<60>>;
using hour_t   = si::quantity_scale<minute_t, si::make_ratio<60>>;
using day_t    = si::quantity_scale<hour_t,   si::make_ratio<24>>;

using millisecond_t = si::milli<second_t>;
CREATE_SUFFIX(millisecond_t, ms)

using microsecond_t = si::micro<second_t>;
CREATE_SUFFIX(microsecond_t, us)

// Frequency units
using megahertz_t = si::mega<hertz_t>;
CREATE_SUFFIX(megahertz_t, MHz)

// Velocity units
using meter_per_s_t = si::quantity_divide<meter_t, second_t>;
CREATE_SUFFIX(meter_per_s_t, mps)

// Acceleration units
using meter_per_ss_t = si::quantity_divide<meter_per_s_t, second_t>;
CREATE_SUFFIX(meter_per_ss_t, mpss)

