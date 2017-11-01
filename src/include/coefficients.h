#pragma once

#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace constants {
    using namespace boost::multiprecision;

    typedef cpp_dec_float_50 bmp_real;
    typedef std::vector<cpp_dec_float_50> mp_coefficients;

    // Значение pi
    const bmp_real MP_PI = boost::math::constants::pi<bmp_real>();

    // Значение I1(0)
    const bmp_real I_1_0 = MP_PI * MP_PI / 12.0;

    // Значение I3(0)
    const bmp_real I_3_0 = 7 * MP_PI * MP_PI * MP_PI * MP_PI / 120;

} // namespace constants