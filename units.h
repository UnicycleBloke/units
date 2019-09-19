#include <type_traits>
#include <ratio>
#include <cassert>


namespace si {


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Representation for powers of SI base units. The base units effective form the basis of 
// a vector space for dimensions.
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
typename KG_T = kilogram_u<0>, 
typename M_T = meter_u<0>, 
typename S_T = second_u<0>, 
typename A_T = ampere_u<0>, 
typename K_T = kelvin_u<0>, 
typename MOL_T = mole_u<0>, 
typename CD_T = candela_u<0>>
struct dimension
{
    // There's probably a nicer way to do this. The idea is that a dimension is created 
    // only from the set of seven SI base unit types above. Templates being what they are,
    // alternative types could have been used so long as they had the correct interface.
    // Maybe this restriction is unnecessary, but we want to be sure that the types of 
    // two dimensions in the space of SI units are the same.
    static_assert(std::is_same<KG_T, kilogram_u<KG_T::exp>>::value);
    static_assert(std::is_same<M_T, meter_u<M_T::exp>>::value);
    static_assert(std::is_same<S_T, second_u<S_T::exp>>::value);
    static_assert(std::is_same<A_T, ampere_u<A_T::exp>>::value);
    static_assert(std::is_same<K_T, kelvin_u<K_T::exp>>::value);
    static_assert(std::is_same<MOL_T, mole_u<MOL_T::exp>>::value);
    static_assert(std::is_same<CD_T, candela_u<CD_T::exp>>::value);

    // Provide access to the exponents for each base unit.
    using kg = KG_T;
    using m = M_T;
    using s = S_T;
    using A = A_T;
    using K = K_T;
    using mol = MOL_T;
    using cd = CD_T;
};


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Multiplying or dividing dimensions adds or subtracts the base vectors.
// Could possibly add powers and roots. Roots would appear to require
// changing the base units to be rational rather than integral powers.
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
// kilogram_u<D::kg::exp * N>, 
// meter_u<D::m::exp * N>, 
// second_u<D::s::exp * N>,
// ampere_u<D::A::exp * N>,
// kelvin_u<D::K::exp * N>,
// mole_u<D::mol::exp * N>,
// candela_u<D::cd::exp * N>>;


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

// Convenience function.
template <typename... Ds>
using make_dimension = typename make_dimension_impl<Ds...>::dim;


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
template <typename T, typename Dim, typename Ratio = std::ratio<1>>
class quantity
{
public:
    static_assert(std::is_arithmetic<T>::value);
    static_assert(std::is_same<Ratio, make_ratio<Ratio::num, Ratio::den>>::value);
    static_assert(std::is_same<Dim, dimension< 
        kilogram_u<Dim::kg::exp>, meter_u<Dim::m::exp>, second_u<Dim::s::exp>,
        ampere_u<Dim::A::exp>, kelvin_u<Dim::K::exp>, mole_u<Dim::mol::exp>,
        candela_u<Dim::cd::exp>>>::value);

    using type = T;
    using dim = Dim;
    using ratio = Ratio;

    explicit constexpr quantity(T value) noexcept : m_value{value} {}
    constexpr T value() const noexcept { return m_value; }

    // Copy constructions from quantities with different ratios. 
    template <typename Ratio2>
    explicit constexpr quantity(const quantity<T, Dim, Ratio2>& other) noexcept
    : m_value{other.value() * Ratio2::num * Ratio::den / Ratio2::den / Ratio::num} 
    {
    }

    // Conversion to quantities with different ratios.
    template <typename Ratio2>
    constexpr operator quantity<T, Dim, Ratio2>() const noexcept
    {
        using Q2 = quantity<T, Dim, Ratio2>;
        return Q2(m_value * Ratio::num * Ratio2::den / Ratio::den / Ratio2::num);
    }

    // Non-explicit form leads to operator overload ambiguity - and errors - with the scalar 
    // multiply and divide because the bool implicitly converts to int. Maybe just use !!.
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

    quantity& operator*=(quantity<T, dimension<>> scalar) noexcept 
    {
        m_value *= scalar.value();
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

    quantity& operator/=(quantity<T, dimension<>> scalar) noexcept 
    {
        m_value /= scalar.value();
        return *this;
    }

    // It only makes sense to add-assign with quantities having the same dimension.
    template <typename Ratio2>
    quantity& operator+=(const quantity<T, Dim, Ratio2>& q2) noexcept
    {
        m_value += (q2.value() * Ratio2::num * Ratio::den / Ratio2::den / Ratio::num);
        return *this;
    }

    // It only makes sense to subtract-assign with quantities having the same dimension.
    template <typename Ratio2>
    quantity& operator-=(const quantity<T, Dim, Ratio2>& q2) noexcept
    {
        m_value -= (q2.value() * Ratio2::num * Ratio::den / Ratio2::den / Ratio::num);
        return *this;
    }

private:
    T m_value; 
};


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Multiply two quantities to create a new quantity with the relevant dimension and ratio.
template <typename T, typename D1, typename D2, typename R1, typename R2>
constexpr auto operator*(const quantity<T, D1, R1>& q1, const quantity<T, D2, R2>& q2) noexcept 
// Is it really necessary to duplicated return type? auto works with deduction, so no, but
// I am not certain that it will return Q rather than Q& or whatever.
    -> quantity<T, dimension_multiply<D1, D2>, std::ratio_multiply<R1, R2>>
{
    using Q = quantity<T, dimension_multiply<D1, D2>, std::ratio_multiply<R1, R2>>;
    return Q(q1.value() * q2.value());
}

// Multiply by a scalar.
template <typename T, typename D, typename R, typename U>
constexpr quantity<T, D, R> operator*(const quantity<T, D, R>& q, U scalar) noexcept 
{
    static_assert(std::is_arithmetic<U>::value);
    return quantity<T, D, R>(q.value() * scalar);
}

// Multiply by a scalar.
template <typename T, typename D, typename R, typename U>
constexpr quantity<T, D, R> operator*(U scalar, const quantity<T, D, R>& q) noexcept 
{
    static_assert(std::is_arithmetic<U>::value);
    return quantity<T, D, R>(q.value() * scalar);
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Divide two quantities to create a new quantity with the relevant dimension and ratio.
template <typename T, typename D1, typename D2, typename R1, typename R2>
constexpr auto operator/(const quantity<T, D1, R1>& q1, const quantity<T, D2, R2>& q2) noexcept
    -> quantity<T, dimension_divide<D1, D2>, std::ratio_divide<R1, R2>>
{
    using Q = quantity<T, dimension_divide<D1, D2>, std::ratio_divide<R1, R2>>;
    return Q(q1.value() / q2.value());
}


// Divide by a scalar.
template <typename T, typename D, typename R, typename U>
constexpr quantity<T, D, R> operator/(const quantity<T, D, R>& q, U scalar) noexcept 
{
    static_assert(std::is_arithmetic<U>::value);
    return quantity<T, D, R>(q.value() / scalar);
}

// Divide a scalar.
template <typename T, typename D, typename R, typename U>
constexpr auto operator/(U scalar, const quantity<T, D, R>& q) noexcept 
    -> quantity<T, dimension_divide<scalar_d, D>, std::ratio_divide<std::ratio<1>, R>>
{
    // The dimension and the ratio need to be inverted in this case.
    static_assert(std::is_arithmetic<U>::value);
    return quantity<T, dimension_divide<scalar_d, D>, std::ratio_divide<std::ratio<1>, R>>(scalar / q.value());
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// Add two quantities which have the same dimension but may have different ratios.
template <typename T, typename D1, typename D2, typename R1, typename R2>
constexpr quantity<T, D1, R1> operator+(const quantity<T, D1, R1>& q1, const quantity<T, D2, R2>& q2) noexcept 
{
    return quantity<T, D1, R1>(q1.value() + q2.value() * R2::num * R1::den / R2::den / R1::num);
}

// Subtract two quantities which have the same dimension but may have different ratios.
template <typename T, typename D1, typename D2, typename R1, typename R2>
constexpr quantity<T, D1, R1> operator-(const quantity<T, D1, R1>& q1, const quantity<T, D2, R2>& q2) noexcept 
{
    return quantity<T, D1, R1>(q1.value() - q2.value() * R2::num * R1::den / R2::den / R1::num);
}

template <typename T, typename D, typename R>
constexpr quantity<T, D, R> operator-(const quantity<T, D, R>& q) noexcept 
{
    return quantity<T, D, R>(-q.value());
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename D, typename R1, typename R2>
constexpr bool operator==(const quantity<T, D, R1>& q1, const quantity<T, D, R2>& q2) noexcept
{
    return (q1.value() * R1::num * R2::den) == (q2.value() * R2::num * R1::den);
}

template <typename T, typename D, typename R1, typename R2>
constexpr bool operator!=(const quantity<T, D, R1>& q1, const quantity<T, D, R2>& q2) noexcept
{
    return (q1.value() * R1::num * R2::den) != (q2.value() * R2::num * R1::den);
}

template <typename T, typename D, typename R1, typename R2>
constexpr bool operator<(const quantity<T, D, R1>& q1, const quantity<T, D, R2>& q2) noexcept
{
    return (q1.value() * R1::num * R2::den) < (q2.value() * R2::num * R1::den);
}

template <typename T, typename D, typename R1, typename R2>
constexpr bool operator<=(const quantity<T, D, R1>& q1, const quantity<T, D, R2>& q2) noexcept
{
    return (q1.value() * R1::num * R2::den) <= (q2.value() * R2::num * R1::den);
}

template <typename T, typename D, typename R1, typename R2>
constexpr bool operator>(const quantity<T, D, R1>& q1, const quantity<T, D, R2>& q2) noexcept 
{
    return (q1.value() * R1::num * R2::den) > (q2.value() * R2::num * R1::den);
}

template <typename T, typename D, typename R1, typename R2>
constexpr bool operator>=(const quantity<T, D, R1>& q1, const quantity<T, D, R2>& q2) noexcept 
{
    return (q1.value() * R1::num * R2::den) >= (q2.value() * R2::num * R1::den);
}


} // namespace si {


//using namespace si;


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// It would be nice if integral types worked, too, but not so important. 
// Also there would be issues of truncation with divides and ratios.
using base_type = float; // float;


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//#define MAKE_DIMENSION(dim, kg, m, s, A, K, mol, cd) \
//using dim = make_dimension< \
// kilogram_u<kg>, meter_u<m>, second_u<s>, ampere_u<A>, kelvin_u<K>, mole_u<mol>, candela_u<cd>>;
using kilogram_d     = si::make_dimension<si::kilogram_u<1>>;
using meter_d        = si::make_dimension<si::meter_u<1>>;
using second_d       = si::make_dimension<si::second_u<1>>;
using hertz_d        = si::make_dimension<si::second_u<-1>>;
using meter_per_s_d  = si::dimension_divide<meter_d, second_d>;
using meter_per_ss_d = si::dimension_divide<meter_per_s_d, second_d>;


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//#define MAKE_QUANTITY(name, dim, ratio, suffix) \
//using name = quantity<double, dim, ratio>; \
//constexpr name operator""_ ## suffix(long double value) noexcept { return name(value); } \
//constexpr name operator""_ ## suffix(unsigned long long int value) noexcept { return name(value); }

// Mass units
using kilogram_t = si::quantity<base_type, kilogram_d>;
constexpr kilogram_t operator""_kg(long double value) noexcept { return kilogram_t(value); }
constexpr kilogram_t operator""_kg(unsigned long long int value) noexcept { return kilogram_t(value); }

using gram_t = si::quantity<base_type, kilogram_d, std::milli>;
constexpr gram_t operator""_g(long double value) noexcept { return gram_t(value); }
constexpr gram_t operator""_g(unsigned long long int value) noexcept { return gram_t(value); }

using pound_t = si::quantity<base_type, kilogram_d, si::make_ratio<454, 1000>>;
constexpr pound_t operator""_lb(long double value) noexcept { return pound_t(value); }
constexpr pound_t operator""_lb(unsigned long long int value) noexcept { return pound_t(value); }

// Distance units
using meter_t = si::quantity<base_type, meter_d>;
constexpr meter_t operator""_m(long double value) noexcept { return meter_t(value); }
constexpr meter_t operator""_m(unsigned long long int value) noexcept { return meter_t(value); }

using kilometer_t = si::quantity<base_type, meter_d, std::kilo>;
constexpr kilometer_t operator""_km(long double value) noexcept { return kilometer_t(value); }
constexpr kilometer_t operator""_km(unsigned long long int value) noexcept { return kilometer_t(value); }

using centimeter_t = si::quantity<base_type, meter_d, std::centi>;
constexpr centimeter_t operator""_cm(long double value) noexcept { return centimeter_t(value); }
constexpr centimeter_t operator""_cm(unsigned long long int value) noexcept { return centimeter_t(value); }

using millimeter_t = si::quantity<base_type, meter_d, std::milli>;
constexpr millimeter_t operator""_mm(long double value) noexcept { return millimeter_t(value); }
constexpr millimeter_t operator""_mm(unsigned long long int value) noexcept { return millimeter_t(value); }

// Time units
using second_t = si::quantity<base_type, second_d>;
constexpr second_t operator""_s(long double value) noexcept { return second_t(value); }
constexpr second_t operator""_s(unsigned long long int value) noexcept { return second_t(value); }

using minute_t = si::quantity<base_type, second_d, si::make_ratio<60>>;
using hour_t = si::quantity<base_type, second_d, si::make_ratio<3600>>;
using day_t = si::quantity<base_type, second_d, si::make_ratio<86400>>;

using millisecond_t = si::quantity<base_type, second_d, std::milli>;
constexpr millisecond_t operator""_ms(long double value) noexcept { return millisecond_t(value); }
constexpr millisecond_t operator""_ms(unsigned long long int value) noexcept { return millisecond_t(value); }

using microsecond_t = si::quantity<base_type, second_d, std::micro>;
constexpr microsecond_t operator""_us(long double value) noexcept { return microsecond_t(value); }
constexpr microsecond_t operator""_us(unsigned long long int value) noexcept { return microsecond_t(value); }

// Frequency units
using hertz_t = si::quantity<base_type, hertz_d>;
constexpr hertz_t operator""_Hz(long double value) noexcept { return hertz_t(value); }
constexpr hertz_t operator""_Hz(unsigned long long int value) noexcept { return hertz_t(value); }

using megahertz_t = si::quantity<base_type, hertz_d, std::mega>;
constexpr megahertz_t operator""_MHz(long double value) noexcept { return megahertz_t(value); }
constexpr megahertz_t operator""_MHz(unsigned long long int value) noexcept { return megahertz_t(value); }

// Velocity units
using meter_per_s_t = si::quantity<base_type, meter_per_s_d>;
constexpr meter_per_s_t operator""_mps(long double value) noexcept { return meter_per_s_t(value); }
constexpr meter_per_s_t operator""_mps(unsigned long long int value) noexcept { return meter_per_s_t(value); }

// Acceleration units
using meter_per_ss_t = si::quantity<base_type, meter_per_ss_d>;
constexpr meter_per_ss_t operator""_mpss(long double value) noexcept { return meter_per_ss_t(value); }
constexpr meter_per_ss_t operator""_mpss(unsigned long long int value) noexcept { return meter_per_ss_t(value); }



template <typename Q, typename QTag>
struct tagged_quantity { Q value; };





using radius_t = tagged_quantity<meter_t, struct radius_tag>;

void circle(radius_t r)
{

}

void compile_time_tests()
{
    radius_t r {23_m};
    circle(r);

    meter_t r2{23};
    circle(r2);

    radius_t r3{r2};
    circle(r3);

    meter_t m3{100};
    m3 *= si::quantity<base_type, si::scalar_d>(3);
    assert(m3.value() == 300);
    m3 /= si::quantity<base_type, si::scalar_d>(3);
    assert(m3.value() == 100);

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
    static_assert(bool(meter_t{100}) == true);
    static_assert(bool(meter_t{0}) == false);
    static_assert(false != bool(meter_t{100}));
    static_assert(true != bool(meter_t{0}));
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
    static_assert(12_mps * 4 == 48_mps); // Ambiguous apparently
    static_assert(4 * 12_mps == 48_mps); // Ambiguous apparently
    static_assert(12_mps / 4 == 3_mps); // Ambiguous apparently
    static_assert(12 / 4_s == 3_Hz); // Ambiguous apparently

    // Addition and subtraction
    static_assert(12_s + 3_ms == 12.003_s);
    static_assert(12_us + 3_ms == 3.012_ms);
    static_assert(12_s - 3_ms == 11.997_s);
    static_assert(1_ms - 3_ms == -2.0_ms);
    static_assert(3_ms - 1_s == -0.997_s);
    static_assert(3_ms - 1_s == -997_ms);

    // Assignment with arithmetic - these can't be constexpr, I think, as they modify the target object.
    meter_t m2(123);
    assert(m2 == 123_m);
    m2 *= 10;
    assert(m2 == 1230_m);
}


base_type test_double(base_type m, base_type s)
{ 
    return m * m * m / s; 
}


auto test_units(meter_t m, second_t s)
{
    return m * m * m / s;
}
