#pragma once

#include "coefficients.h"
#include <vector>


namespace fdsf{
    using namespace constants;

    // Вычисляет значение функции ФД индекса k=1 в точке x
    bmp_real fermi_dirak_k1(bmp_real x);

    // Вычисляет значение функции ФД индекса k=2 в точке x
    bmp_real fermi_dirak_k2(bmp_real x);

    // Вычисляет значение функции ФД индекса k=3 в точке x
    bmp_real fermi_dirak_k3(bmp_real x);

    // Вычисляет значение функции ФД индекса k=4 в точке x
    bmp_real fermi_dirak_k4(bmp_real x);

    // Вычисляет значение функции ФД индекса k = -3/2 в точке x
    bmp_real fermi_dirak_m_3half(bmp_real x);

    // Вычисляет значение функции ФД индекса k = -1/2 в точке x
    bmp_real fermi_dirak_m_half(bmp_real x);

    // Вычисляет значение функции ФД индекса k = 1/2 в точке x
    bmp_real fermi_dirak_half(bmp_real x);

    // Вычисляет значение функции ФД индекса k = 3/2 в точке x
    bmp_real fermi_dirak_3half(bmp_real x);

    // Вычисляет значение функции ФД индекса k = 5/2 в точке x
    bmp_real fermi_dirak_5half(bmp_real x);

    // Вычисляет значение функции ФД индекса k = 7/2 в точке x
    bmp_real fermi_dirak_7half(bmp_real x);

} // namespace fdsf