#include "fdsf.h"
#include "coefficients.h"
#include <boost/math/special_functions/gamma.hpp>
#include <iomanip>
#include <limits>

using namespace fdsf;

// Значение pi
const cpp_dec_float_50 PI = boost::math::constants::pi<cpp_dec_float_50>();

// Значение Ik(0) для индекса k = 0,1,2,3
//const static cpp_dec_float_50 I_k_0[] = { log(cpp_dec_float_50(2))
#if 0
cpp_dec_float_50 fdsf::get_T_max(cpp_dec_float_50 X, int k)
{
    cpp_dec_float_50 a1;
    if (k != 0) {
        a1 = pow(fdsf::factorial(k + 1), -1 / k);
    }
    int i = 1;
    cpp_dec_float_50 y = log(1 + exp(X));
    cpp_dec_float_50 T = X - log(epsilon), I_approximate;
    while (true)
    {
        // 3-х итераций вполне достаточно для определения Tmax
        if (i > 3) {
            break;
        }
        if (k == 0) {
            I_approximate = y;
        }
        else {
            I_approximate = fdsf::factorial(k)*y*pow(1 + a1*y, k);
        }
        // Итерационное вычисление Tmax
        T = X - log(epsilon) + log(pow(T + 1, k) / I_approximate);
        i++;
    }

    return sqrt(T);
}
#endif

// Новая схема Горнера для прецизионного вычисления y
static cpp_dec_float_50 GornerSchemeForPrecesionY(int N, cpp_dec_float_50 x)
{
    //cpp_dec_float_50 alpha = x < epsilon ? exp(x) : exp(-x);
    const cpp_dec_float_50 alpha = exp(x);
    const cpp_dec_float_50 z = alpha / (2 + alpha);
    cpp_dec_float_50 sum = 1.0 / (2 * N + 1);

    for (int i = N - 1; i >= 0; i--) {
        sum = 1.0 / (2 * i + 1.0) + z*z*sum;
    }

    //return x < epsilon ? 2*z*sum : x + 2*z*sum;
    return 2 * z*sum;
}

cpp_dec_float_50 fdsf::FermiDirakFunction(cpp_dec_float_50 t, cpp_dec_float_50 x, cpp_dec_float_50 k)
{
#if 0
    // new , зачин для полуцелых
    return 2 * pow(t, 2 * k + 1) / (1 + exp(t*t - x));
#endif
    return pow(t, k) / (boost::math::tgamma(cpp_dec_float_50(k))*(exp(x) + exp(t)));
}

// хочется сделать одну функцию, которая будет подтягивать коэффиуциенты в зависимости от индекса
static cpp_dec_float_50 get_negative_value(const cpp_dec_float_50 x, const cpp_dec_float_50 k)
{
    using namespace coefficients;
    const size_t N_base = 4;
    cpp_dec_float_50 S1 = 0, S2 = 0;
    cpp_dec_float_50 I_negative;
    cpp_dec_float_50 y = log(1 + exp(x));

    for (int n = 0; n < N_base + 1; n++) {
        S1 = S1 + a_k1[n]*pow(y, n + 1);
    }
    for (int m = 0; m < N_base; m++) {
        S2 = S2 + b_k1[m]*pow(y, m + 1);
    }

    I_negative = (1 + S1) / (1 + S2);
    //delta_base.at(j) = (F_base.at(j) / z.at(j) - 1);
    return I_negative;
}

cpp_dec_float_50 fdsf::integer::fd_one(cpp_dec_float_50 x)
{
    const cpp_dec_float_50 I_1_0 = PI*PI / 12;
    const int k = 1;
  
    cpp_dec_float_50 I_1_minus_x = get_negative_value(-x, k);

    return x*x / 2 + 2*I_1_0 - I_1_minus_x;
}

cpp_dec_float_50 fdsf::integer::fd_two(cpp_dec_float_50 x)
{
    const cpp_dec_float_50 I_1_0 = PI*PI / 12.0;
    const int k = 2;

    cpp_dec_float_50 I_2_minus_x = get_negative_value(-x, k);

    return x*x*x / 3 + 4*x*I_1_0 + I_2_minus_x;
}

cpp_dec_float_50 fdsf::integer::fd_3(cpp_dec_float_50 x)
{
    const cpp_dec_float_50 I_1_0 = PI*PI / 12.0;
    const cpp_dec_float_50 I_3_0 = 7 * PI*PI*PI*PI / 120;
    const int k = 3;

    cpp_dec_float_50 I_3_minus_x = get_negative_value(-x, k);

    return x*x*x*x / 4 + 6*x*x*I_1_0 + 2*I_3_0 - I_3_minus_x;
}

cpp_dec_float_50 fdsf::integer::fd_4(const cpp_dec_float_50 x)
{
    const cpp_dec_float_50 I_1_0 = PI*PI / 12.0;
    const cpp_dec_float_50 I_3_0 = 7 * PI*PI*PI*PI / 120;
    const int k = 4;

    cpp_dec_float_50 I_4_minus_x = get_negative_value(-x, k);
    return x*x*x*x*x/5 + 8*x*x*x*I_1_0 + 8*x*I_3_0 + I_4_minus_x;
}
