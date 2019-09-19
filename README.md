# ub-units - A simple-ish units library for C++17 

This is a compile time modern C++ library which provides dimensional analysis and automatic scaling of physical units. The design is based around the idea that each dimension (such as m/s^2 for acceleration) is represented by a unique point in the vector space of which the seven fundamental SI base units form the basis. Each type of quantity is associated at compile time with a dimension and a scale factor (from std::ratio). Dimensions are simple to create, and are automatically calculated for arithmetic expressions.

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
```