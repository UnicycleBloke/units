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

// template <typename D, int N>
// using dimension_power = dimension<
//     kilogram_u<D::kg::exp * N>, 
//     meter_u<D::m::exp * N>, 
//     second_u<D::s::exp * N>,
//     ampere_u<D::A::exp * N>,
//     kelvin_u<D::K::exp * N>,
//     mole_u<D::mol::exp * N>,
//     candela_u<D::cd::exp * N>>;


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Type function to generate a single component of a dimension based on an argument 
// from a type list. The argument represents a single SI base unit, to some power).
namespace detail {

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
// This is used to invert dimensions below.
using scalar_d = make_dimension<>;


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
    using tag   = Tag;
    using dim   = D;
    using ratio = R;

    explicit constexpr quantity(T value) noexcept : m_value{value} {}

    template <typename R2>
    explicit constexpr quantity(const quantity<T, D, R2, Tag>& other) noexcept
    : m_value{other.value() * R2::num * R::den / R2::den / R::num} 
    {
    }

    constexpr T value() const noexcept { return m_value; }

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
    template <typename U>
    quantity& operator*=(U scalar) noexcept 
    {
        static_assert(std::is_arithmetic<U>::value);
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
    template <typename U>
    quantity& operator/=(U scalar) noexcept 
    {
        static_assert(std::is_arithmetic<U>::value);
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

template <typename Q1, typename Q2>
using quantity_multiply = decltype(std::declval<Q1>() * std::declval<Q2>());

template <typename Q1, typename Q2>
using quantity_divide = decltype(std::declval<Q1>() / std::declval<Q2>());


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
template <typename T, typename Tag, typename D, typename R, typename U>
constexpr quantity<T, D, R, Tag> operator*(const quantity<T, D, R, Tag>& q, U scalar) noexcept 
{
    static_assert(std::is_arithmetic<U>::value);
    return quantity<T, D, R, Tag>(q.value() * scalar);
}

// Multiply by a scalar.
template <typename T, typename Tag, typename D, typename R, typename U>
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
template <typename T, typename Tag, typename D, typename R, typename U>
constexpr quantity<T, D, R, Tag> operator/(const quantity<T, D, R, Tag>& q, U scalar) noexcept 
{
    static_assert(std::is_arithmetic<U>::value);
    return quantity<T, D, R, Tag>(q.value() / scalar);
}

// Divide a scalar.
template <typename T, typename Tag, typename D, typename R, typename U>
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
template <typename T, typename Tag, typename D1, typename D2, typename R1, typename R2>
constexpr quantity<T, D1, R1, Tag> operator+(const quantity<T, D1, R1, Tag>& q1, const quantity<T, D2, R2, Tag>& q2) noexcept 
{
    return quantity<T, D1, R1, Tag>(q1.value() + q2.value() * R2::num * R1::den / R2::den / R1::num);
}

// Subtract two quantities which have the same dimension but may have different ratios.
template <typename T, typename Tag, typename D1, typename D2, typename R1, typename R2>
constexpr quantity<T, D1, R1, Tag> operator-(const quantity<T, D1, R1, Tag>& q1, const quantity<T, D2, R2, Tag>& q2) noexcept 
{
    return quantity<T, D1, R1, Tag>(q1.value() - q2.value() * R2::num * R1::den / R2::den / R1::num);
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
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> exp(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::exp(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> exp2(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::exp2(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> expm1(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::expm1(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> log(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::log(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> log10(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::log10(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> log2(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::log2(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> log1p(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::log1p(q.value() * R::num / R::den));
}


// Power functions - TODO pow, sqrt and cbrt
template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, D, R1, Tag> hypot(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return quantity<T, D, R1, Tag>(std::hypot(q1.value(), q2.value() * R2::num * R1::den / R2::den / R1::num));
}


// Trigonometric functions - can take only scalar quantities, and the ratio is normalised.
template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> sin(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::sin(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> cos(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::cos(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> tan(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::tan(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> asin(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::asin(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> acos(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::acos(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> atan(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::atan(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename D, typename R1, typename R2>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> atan2(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::atan2(q1.value() * R1::num / R1::den, q2.value() * R2::num / R2::den));
}


// Hyperbolic functions - can take only scalar quantities, and the ratio is normalised.
template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> sinh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::sinh(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> cosh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::cosh(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> tanh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::tanh(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> asinh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::asinh(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> acosh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::acosh(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> atanh(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::atanh(q.value() * R::num / R::den));
}


//Error and gamma functions 
template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> erf(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::erf(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> erfc(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::erfc(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> tgamma(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::tgamma(q.value() * R::num / R::den));
}

template <typename T, typename Tag, typename R>
constexpr quantity<T, scalar_d, std::ratio<1>, Tag> lgamma(const quantity<T, scalar_d, R, Tag>& q) noexcept
{
    return quantity<T, scalar_d, std::ratio<1>, Tag>(std::lgamma(q.value() * R::num / R::den));
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


//using namespace si;


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
using hertz_d        = si::make_dimension<si::second_u<-1>>;
using meter_per_s_d  = si::dimension_divide<meter_d, second_d>;
using meter_per_ss_d = si::dimension_divide<meter_per_s_d, second_d>;


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Mass units
// The kilogram is the base unit, but this messes with the standard scale factors milli, kilo, ...
// so we create the gram first.
using gram_t = si::quantity<base_type, kilogram_d>;
constexpr gram_t operator""_g(long double value) noexcept { return gram_t(value); }
constexpr gram_t operator""_g(unsigned long long int value) noexcept { return gram_t(value); }

//using milligram_t = si::quantity_scale<gram_t, std::milli>;
using milligram_t = si::milli<gram_t>;
constexpr milligram_t operator""_mg(long double value) noexcept { return milligram_t(value); }
constexpr milligram_t operator""_mg(unsigned long long int value) noexcept { return milligram_t(value); }

using kilogram_t = si::kilo<gram_t>;
constexpr kilogram_t operator""_kg(long double value) noexcept { return kilogram_t(value); }
constexpr kilogram_t operator""_kg(unsigned long long int value) noexcept { return kilogram_t(value); }

using pound_t = si::quantity_scale<gram_t, si::make_ratio<454>>;
constexpr pound_t operator""_lb(long double value) noexcept { return pound_t(value); }
constexpr pound_t operator""_lb(unsigned long long int value) noexcept { return pound_t(value); }

static_assert(1_kg == 1000_g);
static_assert(1_lb == 454_g);

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


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Compile time tests

// abs(const quantity<T, D, R, Tag>& q) noexcept
static_assert(si::abs(-123_m) == 123_m);
static_assert(si::abs(-123_m) != -123_m);
// fmod(const quantity<T, D, R, Tag>& q) noexcept
// remainder(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
// remquo(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2, int* quo) noexcept 
// fmax(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
static_assert(si::fmax(12_m, 13_m) == 13_m);
// fmin(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 
static_assert(si::fmin(12_m, 13_m) == 12_m);
// fdim(const quantity<T, D, R1, Tag>& q1, const quantity<T, D, R2, Tag>& q2) noexcept 

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


base_type test_double(base_type m, base_type s)
{ 
    return m * m * m / s; 
}


auto test_units(meter_t m, second_t s)
{
    return m * m * m / s;
}
