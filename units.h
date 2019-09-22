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
#include <type_traits>
#include <ratio>
#include <cassert>
#include <cmath>
#include <complex>


namespace si {


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Representation for powers of SI base units. The base units effectively form the basis of 
// a vector space (R7) for dimensions.
template <int EXP> struct kilogram_u { static constexpr int exp = EXP; };
template <int EXP> struct meter_u { static constexpr int exp = EXP; };
template <int EXP> struct second_u { static constexpr int exp = EXP; };
template <int EXP> struct ampere_u { static constexpr int exp = EXP; };
template <int EXP> struct kelvin_u { static constexpr int exp = EXP; };
template <int EXP> struct mole_u { static constexpr int exp = EXP; };
template <int EXP> struct candela_u { static constexpr int exp = EXP; };


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// A dimension for a quantity is basically a point in the vector space of which the units 
// form the basis. 
template <
    typename KG_T  = kilogram_u<0>, 
    typename M_T   = meter_u<0>, 
    typename S_T   = second_u<0>, 
    typename A_T   = ampere_u<0>, 
    typename K_T   = kelvin_u<0>, 
    typename MOL_T = mole_u<0>, 
    typename CD_T  = candela_u<0>>
struct dimension
{
    // There's probably a nicer way to do this. The idea is that a dimension is created 
    // only from the set of seven SI base unit types above. Templates being what they are,
    // alternative types could have been used so long as they had the correct interface.
    // Maybe this restriction is unnecessary, but we want to be sure that the types of 
    // two dimensions for the same point in the space of SI units are the same.
    static_assert(std::is_same_v<KG_T, kilogram_u<KG_T::exp>>);
    static_assert(std::is_same_v<M_T, meter_u<M_T::exp>>);
    static_assert(std::is_same_v<S_T, second_u<S_T::exp>>);
    static_assert(std::is_same_v<A_T, ampere_u<A_T::exp>>);
    static_assert(std::is_same_v<K_T, kelvin_u<K_T::exp>>);
    static_assert(std::is_same_v<MOL_T, mole_u<MOL_T::exp>>);
    static_assert(std::is_same_v<CD_T, candela_u<CD_T::exp>>);

    // Provide access to the exponents for each base unit.
    using kg  = KG_T;
    using m   = M_T;
    using s   = S_T;
    using A   = A_T;
    using K   = K_T;
    using mol = MOL_T;
    using cd  = CD_T;
};


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// This is used to invert dimensions below.
using scalar_d = dimension<>;


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Multiplying or dividing dimensions adds or subtracts the corresponding vectors.
// Could possibly add powers and roots. Roots would appear to require
// changing the base unit exponents to be rational rather than integral powers.
// Then we could use rational powers, too...
template <typename D1, typename D2>
using dimension_multiply = dimension<
    kilogram_u<D1::kg::exp + D2::kg::exp>, 
    meter_u<D1::m::exp + D2::m::exp>, 
    second_u<D1::s::exp + D2::s::exp>,
    ampere_u<D1::A::exp + D2::A::exp>,
    kelvin_u<D1::K::exp + D2::K::exp>,
    mole_u<D1::mol::exp + D2::mol::exp>,
    candela_u<D1::cd::exp + D2::cd::exp>>;

template <typename D1, typename D2>
using dimension_divide = dimension<
    kilogram_u<D1::kg::exp - D2::kg::exp>, 
    meter_u<D1::m::exp - D2::m::exp>, 
    second_u<D1::s::exp - D2::s::exp>,
    ampere_u<D1::A::exp - D2::A::exp>,
    kelvin_u<D1::K::exp - D2::K::exp>,
    mole_u<D1::mol::exp - D2::mol::exp>,
    candela_u<D1::cd::exp - D2::cd::exp>>;


namespace detail {


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Type function to generate a single component of a dimension based on an argument 
// from a type list. The argument represents a single SI base unit, to some power).
template <typename U>
struct get_dimension_impl;

// Convenience function.
template <typename U>
using get_dimension = typename get_dimension_impl<U>::dim;

template <int N>
struct get_dimension_impl<kilogram_u<N>>
{
    using dim = dimension<kilogram_u<N>>;
};

template <int N>
struct get_dimension_impl<meter_u<N>>
{
    using dim = dimension<kilogram_u<0>, meter_u<N>>;
};

template <int N>
struct get_dimension_impl<second_u<N>>
{
    using dim = dimension<kilogram_u<0>, meter_u<0>, second_u<N>>;
};

template <int N>
struct get_dimension_impl<ampere_u<N>>
{
    using dim = dimension<kilogram_u<0>, meter_u<0>, second_u<0>, ampere_u<N>>;
};

template <int N>
struct get_dimension_impl<kelvin_u<N>>
{
    using dim = dimension<kilogram_u<0>, meter_u<0>, second_u<0>, ampere_u<0>, kelvin_u<N>>;
};

template <int N>
struct get_dimension_impl<mole_u<N>>
{
    using dim = dimension<kilogram_u<0>, meter_u<0>, second_u<0>, ampere_u<0>, kelvin_u<0>, mole_u<N>>;
};

template <int N>
struct get_dimension_impl<candela_u<N>>
{
    using dim = dimension<kilogram_u<0>, meter_u<0>, second_u<0>, ampere_u<0>, kelvin_u<0>, mole_u<0>, candela_u<N>>;
};


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Type function to convert a type list of base units (with powers) in any order into a 
// dimension. A dimension is some point in the vector space of possible combinations of 
// SI base units.
template <typename... Ds>
struct make_dimension_impl
{
    // Scalar type for when there are no template arguments.
    using dim = dimension<>;
};

template <typename D, typename... Ds>
struct make_dimension_impl<D, Ds...>
{
    // Recursively combine the base units. This is all compile time magic. 
    using dim = dimension_multiply<get_dimension<D>, typename make_dimension_impl<Ds...>::dim>;
};


} // namespace detail {


// Convenience function.
template <typename... Ds>
using make_dimension = typename detail::make_dimension_impl<Ds...>::dim;


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Convenience type function to ensure that rationals are reduced properly. Unreduced 
// rationals are allowed with std::ratio, but these are distinct types from the reduced 
// versions. We don't want that. 
template <intmax_t N, intmax_t D = 1>
using make_ratio = typename std::ratio<N, D>::type;


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Represents a physical quantity with some particular dimension and ratio. It contains a 
// single arithmetic data member, typically of a floating point type. A whole family of 
// types with the same dimension can be created: e.g. m, km, mm, cm, ...
template <typename T, typename D, typename R = std::ratio<1>, typename Tag = void>
class quantity
{
public:
    // Simple restriction to keep life easier. Could extend to user defined types that have 
    // the relevant properties. Would be nice to extend to integers and complex values, if 
    // nothing else. 
    static_assert(std::is_floating_point_v<T>);
    // Ensure that the ratio is reduced and makes sense.
    static_assert(std::is_same_v<R, make_ratio<R::num, R::den>>);
    static_assert(R::num * R::den > 0);
    // Only si::dimension types are allowed.
    static_assert(std::is_same_v<D, dimension< 
        kilogram_u<D::kg::exp>, meter_u<D::m::exp>, second_u<D::s::exp>,
        ampere_u<D::A::exp>, kelvin_u<D::K::exp>, mole_u<D::mol::exp>,
        candela_u<D::cd::exp>>>);

    using type  = T;
    using dim   = D;
    using ratio = R;
    using tag   = Tag;

    explicit constexpr quantity(T value) noexcept : m_value{value} {}

    template <typename R2>
    explicit constexpr quantity(const quantity<T, D, R2, Tag>& other) noexcept
    : m_value{other.value() * R2::num * R::den / R2::den / R::num} 
    {
    }

    // Note that this ignores the ratio. The result is 1.23 for 1.23_kg and 1.23_mg.
    constexpr T value() const noexcept { return m_value; }
    constexpr T scaled_value() const noexcept { return m_value * R::num / R::den; }

    // Conversion to quantities with different ratios.
    template <typename R2>
    constexpr operator quantity<T, D, R2, Tag>() const noexcept
    {
        using Q2 = quantity<T, D, R2, Tag>;
        return Q2(m_value * R::num * R2::den / R::den / R2::num);
    }

    // This is problematic unless explicit, because implicit conversion to bool
    // will allow comparisons to succeed that ought to fail.
    explicit constexpr operator bool() const noexcept
    {
        return m_value != 0;
    }

    constexpr bool operator!() const noexcept
    {
        return m_value == 0;
    }

    // It only makes sense to multiply-assign with a scalar. 
    template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
    quantity& operator*=(U scalar) noexcept 
    {
        //static_assert(std::is_arithmetic<U>::value);
        m_value *= scalar;
        return *this;
    }

    template <typename R2>
    quantity& operator*=(quantity<T, scalar_d, R2, Tag> scalar) noexcept 
    {
        m_value *= (scalar.value() * R2::num / R2::den);
        return *this;
    }

    // It only makes sense to divide-assign with a scalar. 
    template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
    quantity& operator/=(U scalar) noexcept 
    {
        //static_assert(std::is_arithmetic<U>::value);
        m_value /= scalar;
        return *this;
    }

    template <typename R2>
    quantity& operator/=(quantity<T, scalar_d, R2, Tag> scalar) noexcept 
    {
        m_value /= (scalar.value() * R2::num / R2::den);
        return *this;
    }

    // It only makes sense to add-assign with quantities having the same dimension.
    template <typename R2>
    quantity& operator+=(const quantity<T, D, R2, Tag>& q2) noexcept
    {
        m_value += (q2.value() * R2::num * R::den / R2::den / R::num);
        return *this;
    }

    // It only makes sense to subtract-assign with quantities having the same dimension.
    template <typename R2>
    quantity& operator-=(const quantity<T, D, R2, Tag>& q2) noexcept
    {
        m_value -= (q2.value() * R2::num * R::den / R2::den / R::num);
        return *this;
    }

private:
    T m_value; 
};


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Simple type function to create a famile of types with the same dimension and tag, but 
// different ratios.
template <typename Q, typename R2>
using quantity_scale = quantity<typename Q::type, typename Q::dim, 
    std::ratio_multiply<typename Q::ratio, R2>, typename Q::tag>;

template <typename Q1, typename Q2>
using quantity_multiply = decltype(std::declval<Q1>() * std::declval<Q2>());

template <typename Q1, typename Q2>
using quantity_divide = decltype(std::declval<Q1>() / std::declval<Q2>());

// yotta, zetta, zepto and yocto are only supported on some platforms.
//template <typename Q> using yotta = quantity_scale<Q, std::yotta>;
//template <typename Q> using zetta = quantity_scale<Q, std::zetta>;
template <typename Q> using exa   = quantity_scale<Q, std::exa>;
template <typename Q> using peta  = quantity_scale<Q, std::peta>;
template <typename Q> using tera  = quantity_scale<Q, std::tera>;
template <typename Q> using giga  = quantity_scale<Q, std::giga>;
template <typename Q> using mega  = quantity_scale<Q, std::mega>;
template <typename Q> using kilo  = quantity_scale<Q, std::kilo>;
template <typename Q> using hecto = quantity_scale<Q, std::hecto>;
template <typename Q> using deca  = quantity_scale<Q, std::deca>;
template <typename Q> using deci  = quantity_scale<Q, std::deci>;
template <typename Q> using centi = quantity_scale<Q, std::centi>;
template <typename Q> using milli = quantity_scale<Q, std::milli>;
template <typename Q> using micro = quantity_scale<Q, std::micro>;
template <typename Q> using nano  = quantity_scale<Q, std::nano>;
template <typename Q> using pico  = quantity_scale<Q, std::pico>;
template <typename Q> using femto = quantity_scale<Q, std::femto>;
template <typename Q> using atto  = quantity_scale<Q, std::atto>;
//template <typename Q> using zepto = quantity_scale<Q, std::zepto>;
//template <typename Q> using yocto = quantity_scale<Q, std::yocto>;


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Multiply two quantities to create a new quantity with the relevant dimension and ratio.
template <typename T, typename Tag, typename D1, typename D2, typename R1, typename R2>
constexpr auto operator*(const quantity<T, D1, R1, Tag>& q1, const quantity<T, D2, R2, Tag>& q2) noexcept 
// Is it really necessary to duplicated return type? auto works with deduction, so no, but
// I am not certain that it will return Q rather than Q& or whatever.
    -> quantity<T, dimension_multiply<D1, D2>, std::ratio_multiply<R1, R2>, Tag>
{
    using Q = quantity<T, dimension_multiply<D1, D2>, std::ratio_multiply<R1, R2>, Tag>;
    return Q(q1.value() * q2.value());
}

// Multiply by a scalar.
template <typename T, typename Tag, typename D, typename R, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
constexpr quantity<T, D, R, Tag> operator*(const quantity<T, D, R, Tag>& q, U scalar) noexcept 
{
    static_assert(std::is_arithmetic<U>::value);
    return quantity<T, D, R, Tag>(q.value() * scalar);
}

// Multiply by a scalar.
template <typename T, typename Tag, typename D, typename R, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
constexpr quantity<T, D, R, Tag> operator*(U scalar, const quantity<T, D, R, Tag>& q) noexcept 
{
    static_assert(std::is_arithmetic<U>::value);
    return quantity<T, D, R, Tag>(q.value() * scalar);
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Divide two quantities to create a new quantity with the relevant dimension and ratio.
template <typename T, typename Tag, typename D1, typename D2, typename R1, typename R2>
constexpr auto operator/(const quantity<T, D1, R1, Tag>& q1, const quantity<T, D2, R2, Tag>& q2) noexcept
    -> quantity<T, dimension_divide<D1, D2>, std::ratio_divide<R1, R2>, Tag>
{
    using Q = quantity<T, dimension_divide<D1, D2>, std::ratio_divide<R1, R2>, Tag>;
    return Q(q1.value() / q2.value());
}

// Divide by a scalar.
template <typename T, typename Tag, typename D, typename R, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
constexpr quantity<T, D, R, Tag> operator/(const quantity<T, D, R, Tag>& q, U scalar) noexcept 
{
    static_assert(std::is_arithmetic<U>::value);
    return quantity<T, D, R, Tag>(q.value() / scalar);
}

// Divide a scalar.
template <typename T, typename Tag, typename D, typename R, typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
constexpr auto operator/(U scalar, const quantity<T, D, R, Tag>& q) noexcept 
    -> quantity<T, dimension_divide<scalar_d, D>, std::ratio_divide<std::ratio<1>, R>, Tag>
{
    // The dimension and the ratio need to be inverted in this case.
    static_assert(std::is_arithmetic<U>::value);
    using Q = quantity<T, dimension_divide<scalar_d, D>, std::ratio_divide<std::ratio<1>, R>, Tag>;
    return Q(scalar / q.value());
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Add two quantities which have the same dimension but may have different ratios.
template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, D, R1, Tag> operator+(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return quantity<T, D, R1, Tag>(q1.value() + q2.value() * R2::num * R1::den / R2::den / R1::num);
}

// Subtract two quantities which have the same dimension but may have different ratios.
template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, D, R1, Tag> operator-(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return quantity<T, D, R1, Tag>(q1.value() - q2.value() * R2::num * R1::den / R2::den / R1::num);
}

template <typename T, typename Tag, typename D, typename R>
constexpr quantity<T, D, R, Tag> operator-(const quantity<T, D, R, Tag>& q) noexcept 
{
    return quantity<T, D, R, Tag>(-q.value());
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Relational operators
template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr bool operator==(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept
{
    return (q1.value() * R1::num * R2::den) == (q2.value() * R2::num * R1::den);
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr bool operator!=(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept
{
    return (q1.value() * R1::num * R2::den) != (q2.value() * R2::num * R1::den);
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr bool operator<(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept
{
    return (q1.value() * R1::num * R2::den) < (q2.value() * R2::num * R1::den);
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr bool operator<=(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept
{
    return (q1.value() * R1::num * R2::den) <= (q2.value() * R2::num * R1::den);
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr bool operator>(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return (q1.value() * R1::num * R2::den) > (q2.value() * R2::num * R1::den);
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr bool operator>=(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return (q1.value() * R1::num * R2::den) >= (q2.value() * R2::num * R1::den);
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Common mathematical functions
// These all basically pass through the underlying value to the matching function in cmath.

// Basic operations 
template <typename T, typename Tag, typename D, typename R>
constexpr quantity<T, D, R, Tag> abs(const quantity<T, D, R, Tag>& q) noexcept
{
    return quantity<T, D, R, Tag>(std::abs(q.value()));
}

template <typename T, typename Tag, typename D, typename R>
constexpr quantity<T, D, R, Tag> fmod(const quantity<T, D, R, Tag>& q) noexcept
{
    return quantity<T, D, R, Tag>(std::fmod(q.value()));
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, D, R1, Tag> remainder(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return quantity<T, D, R1, Tag>(std::remainder(q1.value(), q2.value() * R2::num * R1::den / R2::den / R1::num));
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, D, R1, Tag> remquo(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2, int* quo) noexcept 
{
    return quantity<T, D, R1, Tag>(std::remquo(q1.value(), q2.value() * R2::num * R1::den / R2::den / R1::num, quo));
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, D, R1, Tag> fmax(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return quantity<T, D, R1, Tag>(std::fmax(q1.value(), q2.value() * R2::num * R1::den / R2::den / R1::num));
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, D, R1, Tag> fmin(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return quantity<T, D, R1, Tag>(std::fmin(q1.value(), q2.value() * R2::num * R1::den / R2::den / R1::num));
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, D, R1, Tag> fdim(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return quantity<T, D, R1, Tag>(std::fdim(q1.value(), q2.value() * R2::num * R1::den / R2::den / R1::num));
}


// Exponential functions - can take only scalar quantities, and the ratio is normalised.
template <typename T, typename Tag, typename R>
constexpr T exp(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::exp(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T exp2(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::exp2(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T expm1(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::expm1(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T log(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::log(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T log10(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::log10(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T log2(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::log2(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T log1p(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::log1p(q.value() * R::num / R::den);
}


// Power functions - TODO pow, sqrt and cbrt
template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, D, R1, Tag> hypot(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return quantity<T, D, R1, Tag>(std::hypot(q1.value(), q2.value() * R2::num * R1::den / R2::den / R1::num));
}


// Trigonometric functions - can take only scalar quantities, and the ratio is normalised.
template <typename T, typename Tag, typename R>
constexpr T sin(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::sin(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T cos(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::cos(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T tan(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::tan(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T asin(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::asin(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T acos(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::acos(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T atan(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::atan(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr T atan2(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return std::atan2(q1.value() * R1::num / R1::den, q2.value() * R2::num / R2::den);
}


// Hyperbolic functions - can take only scalar quantities, and the ratio is normalised.
template <typename T, typename Tag, typename R>
constexpr T sinh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::sinh(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T cosh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::cosh(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T tanh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::tanh(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T asinh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::asinh(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T acosh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::acosh(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T atanh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::atanh(q.value() * R::num / R::den);
}


//Error and gamma functions 
template <typename T, typename Tag, typename R>
constexpr T erf(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::erf(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T erfc(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::erfc(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T tgamma(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::tgamma(q.value() * R::num / R::den);
}

template <typename T, typename Tag, typename R>
constexpr T lgamma(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return std::lgamma(q.value() * R::num / R::den);
}


// Nearest integer floating point operations 
template <typename T, typename Tag, typename D, typename R>
constexpr quantity<T, D, R, Tag> ceil(const quantity<T, D, R, Tag>& q) noexcept 
{
    return quantity<T, D, R, Tag>(std::ceil(q.value()));
}

template <typename T, typename Tag, typename D, typename R>
constexpr quantity<T, D, R, Tag> floor(const quantity<T, D, R, Tag>& q) noexcept 
{
    return quantity<T, D, R, Tag>(std::floor(q.value()));
}

template <typename T, typename Tag, typename D, typename R>
constexpr quantity<T, D, R, Tag> trunc(const quantity<T, D, R, Tag>& q) noexcept 
{
    return quantity<T, D, R, Tag>(std::trunc(q.value()));
}

template <typename T, typename Tag, typename D, typename R>
constexpr quantity<T, D, R, Tag> round(const quantity<T, D, R, Tag>& q) noexcept 
{
    return quantity<T, D, R, Tag>(std::round(q.value()));
}


} // namespace si {


