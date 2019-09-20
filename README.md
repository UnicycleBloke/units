# ub-units - A simple units library for C++17 

**ub-units** is a simple compile-time C++ library which provides dimensional analysis and automatic scaling of physical units. The intention is to support strong types in order to reduce errors in mathematical expressions and when passing arguments to functions. For example, consider the following two functions:

```c++
double speed_native(double distance, double time)
{
    return distance / time;
}

meter_per_sec_t speed_units(meter_t distance, second_t time) 
{ 
    return distance / time; 
}
```

In the case of `speed_native()`, it would be simple to pass the parameters in the wrong order, leading to runtime errors that might be hard to diagnose. There is also an issue of scale: is the distance in feet, fathoms, miles, meters, or parsecs? You could always add comments to inform the user, or give meaningful names to the parameters, but these could be ignored, forgotten or misunderstood. Maybe a two-argument function isn't really a big problem, but suppose it was a ten-argument function...

In the case of `speed_units()`, `meter_t` and `second_t` are not simply typedefs of `double`, but distinct types. This makes it impossible to mix up the order of parameters: the code would not compile. 

These two functions were compiled by Compiler Explorer using ARM GCC with `-std=c++17 -O` enabled, and yielded identical assembler when the representation was `float`, and near identical assembler when it was `double`. The additional type safety is free.

For the issue of scale, `speed_units()` expects a distance in meters, but other types with the same dimension are automatically converted by the relevant multiplier. So we can do this:

```c++
kilometer_t km{75}; 
hour_t h{3};
// Returns meter_per_sec_t{(75 * 1000) / (3 * 3600)}
meter_per_sec_t mps = speed_units(km, h); 

// Literals for the various user-defined types are also permitted. 
// The results of expressions are converted to match the assigned 
// variable's type.
lightyear_per_sec_t = speed_units(42_ly, 3_ms); // 14'000 ly/s
```

Every type is associated at compile time with a dimension and a rational scale factor. Arithmetic expressions return values whose types have appropriate dimensions and scale factors.

## Background

I was inspired to write **ub-units** by [Jonathan Boccara](https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/)'s blog articles on strong types for strong interfaces. I was particularly interested in applying this idea to embedded software, where I have sometimes fallen foul of the errors strong types are intended to address. I wanted simple-to-use types which compile to little more than operations on the underlying representations: typically floats or doubles.

Along the way, I have found a couple of alternative libraries which do very similar things: one from [Nic Holthaus](https://github.com/nholthaus/units), and another from [Mateusz Pusz](https://github.com/mpusz/units). I haven't used these or studied them in any detail. Pusz's library appears to be suggested as a candidate for standardisation and includes all kinds of C++20 features that are beyond my current understanding. He presented a great [talk](https://www.youtube.com/watch?v=wKchCktZPHU) about it at C++Now 2019, which is very interesting, if somewhat challenging for bears of very little brain like myself.

I am a relative beginner in template metaprogamming, and have sought to keep my code as simple as possible so that (a) I understand how it works and, (b) other beginners might understand how it works. To be honest, keeping it simple is my approach to all software development, so this is no surprise. 

I feel certain that there are many TMP bells and whistles which I could add to improve the genericity,readability of compilation errors, type constraints, and so on. I'd like to learn more about that stuff but the goal here was to create a library that is simple to use, and which is entirely made up of compile time constraints so that the compiler can optimise all of it away to basic arithmetic operations on the underlying representation. 

I would like to thank [Matt Godbolt](https://github.com/mattgodbolt) for the amazing Compiler Explorer, which made developing this library a lot easier due to its WYSIWYG feedback on whether my templatey goodness would compile or not.

## Supported operations

All the basic operations for arithmetic types are supported. Bitwise operations and logical operations are not supported.

- `+`, `-`, `+=`, `-=`: for quantities with the same dimension.
- `*`, `/`, `*=`, `/=`: for quantities of any dimension, and scalars.
- `!` returns true if the value of the representation is zero. Otherwise false. Operator `bool()` is explicit to avoid ambiguities with conversion to int.
- `==`, `!=`, `<`, `<=`, `>`, `>=`: for quantities with the same dimension, returning bool
- Most `<cmath>` functions are supported: exponential, trigonometric, hyperbolic and other special functions only accept scalars (all exponents zero).

Almost all operations are `constexpr`.

# Design

The library design is based on the principle that each **dimension** for some physical quantity - such as `m/s^2` for acceleration - can be represented by a unique point in the vector space (**R7**) of which the seven fundamental [SI base units](https://www.bipm.org/en/measurement-units/base-units.html) form the basis. 
- Each class of quantity (e.g. distance, time, speed, ...) is associated at compile time with a dimension. 
- Each such class can be represented by any number of distinct types which have relative scaling ratios (e.g. for mass we have kilograms, grams, milligrams, pounds, daltons, ...). 
- Quantities with the same dimension but different ratios are automatically converted as necessary. 
- Quantities with different dimensions cannot be automatically converted, and result in compile-time errors if used where a different dimension is expected. 
- Arithmetic expressions with quantities results in new quantities which have the appropriate dimensions and ratios, or compilation errors if the operation makes no sense (such as adding distance to time). 

## SI base units

The SI base units are represented by very simple structures which capture the exponent of the unit as a template parameter. For example in the dimension `m/s^2` the `meter_u` unit has an exponent of `1`, and the `second_u` unit has an exponent of `-2`. The others are all zero.

```c++
template <int EXP> struct kilogram_u { static constexpr int exp = EXP; };
template <int EXP> struct meter_u { static constexpr int exp = EXP; };
template <int EXP> struct second_u { static constexpr int exp = EXP; };
template <int EXP> struct ampere_u { static constexpr int exp = EXP; };
template <int EXP> struct kelvin_u { static constexpr int exp = EXP; };
template <int EXP> struct mole_u { static constexpr int exp = EXP; };
template <int EXP> struct candela_u { static constexpr int exp = EXP; };
```

## Creating dimensions from base units

A `dimension` is just a struct with seven associated types which capture the exponents of the seven SI base units. This can be imagined as a unique vector in the space of dimensions.

```c++
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
    using kg  = KG_T;
    using m   = M_T;
    using s   = S_T;
    using A   = A_T;
    using K   = K_T;
    using mol = MOL_T;
    using cd  = CD_T;
};
```

There are a number of simple helper functions to create dimension types more easily.

```c++
using scalar_d        = dimension<>;
// 'make_dimension' is a variadic template which
// recursively multiplies its parameters together.
using meter_d         = make_dimension<meter_u<1>>;
using second_d        = make_dimension<second_u<1>>;
using meter_per_sec_d = make_dimension<meter_u<1>, second_u<-1>>;;
// OR... 'dimension_divide' adds the exponents (i.e. add the vectors)
using meter_per_sec_d = dimension_divide<meter_d, second_d>;
// No idea what this is for... 'dimension_multiply' subtracts the exponents
using kgsq_m_per_s_d  = dimension_multiply<meter_per_sec_d, kilogram_u<2>>;
```

## Creating custom quantities

Quantities are types which associate three elements:

- The underlying type of the representation (typically a floating point type)
- The `dimension` of the quantity
- The scale factor of the quantity (represented by `std::ratio`). The default ratio is `1`.

We can also easily define suffixes for user-defined literals for our quantities. Though I am reluctant to use macros, this code looks as if it would benefit from them, especially if many quantities are going to be defined.

```c++
using meter_t = quantity<double, meter_d>;
// This version is for floating point inputs (e.g. 123.4_m, 123._m)
constexpr meter_t operator""_m(long double value) noexcept 
    { return meter_t(value); }
// This version is for integral inputs: we can omit the trailing 
// period (e.g. 123_m)
constexpr meter_t operator""_m(unsigned long long int value) noexcept 
    { return meter_t(value); }

// Kilometers have the same representation and dimension, but a 
// different multiplier. 
using kilometer_t = quantity<double, meter_d, std::kilo>;
constexpr kilometer_t operator""_km(long double value) noexcept 
    { return kilometer_t(value); }
constexpr kilometer_t operator""_km(unsigned long long int value) noexcept 
    { return kilometer_t(value); }

static_assert(5_km == 5'000_m);    
```

As with dimensions, there are a number of helper functions to make creating quantities easier.

```c++
using gram_t          = quantity_scale<kilogram_t, std::milli>;
using pound_t         = quantity_scale<kilogram_t, make_ratio<454, 1000>>;
using meter_per_s_t   = quantity_divide<meter_t, second_t>;
using meter_per_ss_t  = quantity_divide<meter_per_s_t, second_t>;
using kilgram_meter_t = quantity_multiply<kilogram_t, meter_t>;
```

## Creating tagged quantities

In addition to the representation, dimension and ratio, an additional template type parameter for quantities is used to specify a **tag**. By default this is `void`. The purpose of the tag is to allow us to create type which have the same dimension, but which do not automatically convert into each other. The reason for this is to help disambiguate sitations where the dimension alone is not sufficient. For example, consider a function which calculates the area of a circle:

```c++
meter_sq_t circle_area(meter_t radius) 
{
    constexpr meter_t::type PI = 3.1415...;
    return PI * radius * radius;
}
```

Though the dimensions are checked, it would easy to accidentally pass the diameter or some other value. Any type with length dimension would be automatically converted to meters and we could have subtle errors in the software.

The solution is to create a distinct family of quantities which represent radii. We do this by passing a unique tag type to identify such quantities.

```c++
using diameter_t = quantity<double, meter_d, std::ratio<1>, struct circle_tag>;
using radius_t   = quantity_scale<diameter_t, std::ratio<2>>;

static_assert(radius_t{10} == diameter_t{20});

// Fails to compile because the tags are different.
static_assert(diameter_t{10} == meter_t{10});

// TODO requires work here to convert to meter_sq_t...
meter_sq_t circle_area(radius_t radius); 

radius_t r{3};
circle_area(r); // Calculates PI * r * r;

diameter_t d{6};
circle_area(d); // Calculates PI * (d/2) * (d/2);
```

This is not an ideal solution - the user code construct the radius with an incorrect value. There is only so much we can do...


