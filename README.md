# ub-units - A simple-ish units library for C++17 

**ub-units** is a simple compile-time modern C++ library which provides dimensional analysis and automatic scaling of physical units. The intention is to support strong types in order to reduce errors in mathematical expressions and when passing arguments to functions. For example, consider the following two functions:

```c++
double speed(double distance, double time)
{
    return distance / time;
}

meter_per_sec_t speed(meter_t distance, second_t time) 
{ 
    return distance / time; 
}
```

In the first case, it would be simple to pass the parameters in the wrong order, leading to runtime errors that might be hard to diagnose. There is also an issue of scale: is the distance in feet, miles, meters, or parsecs? You could add comments to inform the user, or give meaningful names to the parameters, but these could be ignored or misunderstood.

In the second case, `meter_t` and `second_t` are not simply typedefs of double, but distinct types. This makes it impossible to mix up the order of parameters: the code would not compile. For the issue of scale, the function expects a distance in meters, but other types with the same dimension are automatically converted as necessary. So we can do this:

```c++
kilometer_t km{75};
hour_t h{3};
meter_per_sec_t mps = speed(km, h); // Returns meter_per_sec_t{(75 * 1000) / (3 * 3600)}
```

## User defined literals

In addition to associating types with dimensions and scale factors, we can create user-defined literals suffixes for any type, so that we can do the following:

```c++
static_assert(kilometer_t{75} == 75_km);
static_assert(kilometer_t{75} == 75'000_m);
```

## Tagged types

There is still a possible problem in some cases, in which knowing the dimension and scale factor of some quantity is not sufficient. For example, consider calculating the area of a circle:

```c++
meter_sq_t circle_area(meter_t radius);
```

What is to prevent us accidentally passing the diameter instead of the radius? We can created simple wrappers around types with dimensions. The wrapper includes a tag which marks the type as not convertible to or from other types with the same dimension, plus an additional scale factor to allow conversions between types with the same tag.

```c++
meter_sq_t circle_area(radius_t radius);

radius_t r{3_m};
circle_area(r); // Calculates pi * r * r;

diameter_t d{6_m};
circle_area(d); // Calculates pi * (d/2) * (d/2);
```

# Design

The library design is based on the idea that each dimension (such as `m/s^2` for acceleration) is represented by a unique point in the vector space (R7) of which the seven fundamental SI base units form the basis. Each type of quantity is associated at compile time with a dimension and a scale factor (using `std::ratio`). Quantities with the same dimension but different ratios are automatically converted as necessary. Quantities with different dimensions cannot be automatically converted, and result in compile-time errors if used where a different dimension is expected. Arithmetic with quantities results in new quantities which have the appropriate dimensions and ratios, or compilation errors if the operation makes no sense (such as adding distance to time). 

## Creating custom quantities

Dimensions and quantities are intended to be simple to create. 

```c++
// Create new dimensions
using meter_d         = si::make_dimension<si::meter_u<1>>;
using second_d        = si::make_dimension<si::second_u<1>>;
using meter_per_sec_d = si::dimension_divide<meter_d, second_d>;

// Create new quantities, and literal suffixes (optional)
using meter_t = si::quantity<double, meter_d>;
constexpr meter_t operator""_m(long double value) noexcept { return meter_t(value); }
constexpr meter_t operator""_m(unsigned long long int value) noexcept { return meter_t(value); }

using kilometer_t = si::quantity<double, meter_d, std::kilo>;
constexpr kilometer_t operator""_km(long double value) noexcept { return kilometer_t(value); }
constexpr kilometer_t operator""_km(unsigned long long int value) noexcept { return kilometer_t(value); }

...

meter_per_sec_t speed(meter_t d, second_t s) { return d / s; }

auto v1 = speed(100_km, 1_h); // returns 100,000 / 3,600 m/s
auto v2 = speed(100_km, 3_m); // does not compile

micrometer um = 1234_mm; // 1234000um
millisecond ms{1000};
auto v3 = speed(um, ms); // returns 1.234 m/s
auto v4 = speed(ms, um); // does not compile

```

## Creating tagged quantities

TODO Additional tagging of types so that, for example radius_t and diameter_t, having the same dimension (m), are distinct types. This property could be orthogonal to the dimensional analysis, or not. Not really sure how arithmetic with types having different tags would work out. Does the tag come with a second ratio so that, for example, a radius can convert to a 
diameter. Or is the tag an additional pseudo-base in the dimension vector... 

The library code itself is entry level template metaprogramming that hopefully other mere mortals can understand. I wanted to avoid a lot of black magic that I didn't really understand sufficiently well. This probably means that the code is not as tightly constrained as it should be for type deductions and the like, but it's a start. 

The following code demonstrates some of the capabilities:

```c++
    // Equivalence of construction methods.
    static_assert(meter_t{100} == 100_m);

    // Combining dimension
    static_assert(200_mps == 4_km / 20000_ms);
    static_assert(1000 / 50_s == 20_Hz);
    static_assert(1000 / 50_us == 20_MHz);
    static_assert(60 * 1000_ms == minute_t{1});

    // Comparison
    static_assert(100_cm == 1000_mm);
    static_assert(100_cm > 999_mm);

    // Conversion of ratios
    static_assert(kilometer_t{10} == 10000000_mm);

    // Multiplication and division, including with scalars 
    static_assert(12_mps == 48_m / 4_s);
    static_assert(12_mps * 4_s == 48_m);
    static_assert(12_mps * 4 == 48_mps); 
    static_assert(4 * 12_mps == 48_mps); 

    // Addition and subtraction
    static_assert(12_s + 3_ms == 12.003_s);
    static_assert(3_ms - 1_s == -997_ms);

    auto speed(meter_t distance, second_t time)
    {
        return distance / time;
    }
    
```