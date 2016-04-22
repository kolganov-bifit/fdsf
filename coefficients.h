#pragma once

#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace coefficients {
    using namespace boost::multiprecision;

    typedef std::vector<cpp_dec_float_50> coefficient;

    const coefficient a_k1 = { 0.2715113138214362780964488,
                               0.0562661238060587633169245,
                               0.0067420740469345689743873,
                               0.0005169505155333205859985,
                               0.0000194771836765773190602 };

    const coefficient b_k1 = { 0.0215113138214352840651584,
                               0.0231105175729721417901084,
                               0.0003669081577365413477999,
                               0.0000610424408732720110769 };

    const coefficient a_k2 = {};

    const coefficient b_k2 = {};

    const coefficient a_k3 = {};

    const coefficient b_k3 = {};

    const coefficient a_k4 = {};

    const coefficient b_k4 = {};

} // namespace coefficients