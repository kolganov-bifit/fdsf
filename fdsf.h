#pragma once

#include <vector>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/number.hpp>

namespace fdsf{
    using namespace boost::multiprecision;

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x при заданном
    // значении t. Функция представлена в виде t^k /(exp(x)+exp(t))
    cpp_dec_float_50 FermiDirakFunction(cpp_dec_float_50 t, cpp_dec_float_50 x, cpp_dec_float_50 k);

    namespace integer {
        // Вычисляет значение функции ФД индекса k=1 в точке x
        cpp_dec_float_50 fd_one(cpp_dec_float_50 x);

        // Вычисляет значение функции ФД индекса k=2 в точке x
        cpp_dec_float_50 fd_two(cpp_dec_float_50 x);

        // Вычисляет значение функции ФД индекса k=3 в точке x
        cpp_dec_float_50 fd_3(cpp_dec_float_50 x);

        // Вычисляет значение функции ФД индекса k=4 в точке x
        cpp_dec_float_50 fd_4(cpp_dec_float_50 x);
    } // namespace integer

} // namespace fdsf