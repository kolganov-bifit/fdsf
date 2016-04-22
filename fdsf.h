#pragma once

#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>

namespace fdsf{
    using namespace boost::multiprecision;

    // ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x ��� ��������
    // �������� t. ������� ������������ � ���� t^k /(exp(x)+exp(t))
    cpp_dec_float_50 FermiDirakFunction(cpp_dec_float_50 t, cpp_dec_float_50 x, cpp_dec_float_50 k);

    namespace integer {
        // ��������� �������� ������� �� ������� k=1 � ����� x
        cpp_dec_float_50 fd_one(cpp_dec_float_50 x);

        // ��������� �������� ������� �� ������� k=2 � ����� x
        cpp_dec_float_50 fd_two(cpp_dec_float_50 x);

        // ��������� �������� ������� �� ������� k=3 � ����� x
        cpp_dec_float_50 fd_3(cpp_dec_float_50 x);

        // ��������� �������� ������� �� ������� k=4 � ����� x
        cpp_dec_float_50 fd_4(cpp_dec_float_50 x);
    } // namespace integer

} // namespace fdsf