#pragma once

#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>

namespace fdsf{
    using namespace boost::multiprecision;

    typedef cpp_dec_float_50 bmp_real;

    namespace integer {
        // ��������� �������� ������� �� ������� k=1 � ����� x
        bmp_real fd_one(bmp_real x);

        // ��������� �������� ������� �� ������� k=2 � ����� x
        bmp_real fd_two(bmp_real x);

        // ��������� �������� ������� �� ������� k=3 � ����� x
        bmp_real fd_3(bmp_real x);

        // ��������� �������� ������� �� ������� k=4 � ����� x
        bmp_real fd_4(bmp_real x);
    } // namespace integer

} // namespace fdsf