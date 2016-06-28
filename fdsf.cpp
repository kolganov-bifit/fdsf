#include "fdsf.h"
#include "coefficients.h"
#include <iomanip>
#include <limits>

namespace fdsf {

    // Значение Ik(0) для индекса k = 0,1,2,3
    //const static bmp_real I_k_0[] = { log(bmp_real(2))
#if 0
    bmp_real get_T_max(bmp_real X, int k)
    {
        bmp_real a1;
        if (k != 0) {
            a1 = pow(fdsf::factorial(k + 1), -1 / k);
        }
        int i = 1;
        bmp_real y = log(1 + exp(X));
        bmp_real T = X - log(epsilon), I_approximate;
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
    static bmp_real GornerSchemeForPrecesionY(int N, bmp_real x)
    {
        //bmp_real alpha = x < epsilon ? exp(x) : exp(-x);
        const bmp_real alpha = exp(x);
        const bmp_real z = alpha / (2 + alpha);
        bmp_real sum = 1.0 / (2 * N + 1);

        for (int i = N - 1; i >= 0; i--) {
            sum = 1.0 / (2 * i + 1.0) + z*z*sum;
        }

        //return x < epsilon ? 2*z*sum : x + 2*z*sum;
        return 2 * z*sum;
    }

    // Вычисляет значение функции ФД индекса k = 1, 2, 3 в точке x при заданном
    // значении t. Функция представлена в виде t^k /(exp(x)+exp(t))
    static bmp_real FermiDirakFunction(bmp_real t, bmp_real x, bmp_real k)
    {
#if 0
        // new , зачин для полуцелых
        return 2 * pow(t, 2 * k + 1) / (1 + exp(t*t - x));
#endif
        return pow(t, k) / (boost::math::tgamma(bmp_real(k))*(exp(x) + exp(t)));
    }

    // хочется сделать одну функцию, которая будет подтягивать коэффиуциенты в зависимости от индекса
    static bmp_real get_negative_value(const bmp_real x, const bmp_real k,
        const constants::mp_coefficients& a,
        const constants::mp_coefficients& b)
    {

        const size_t N_base = 4;
        bmp_real S1 = 0, S2 = 0;
        bmp_real I_negative;
        bmp_real y = log(1 + exp(x));

        using namespace constants;

        for (int n = 0; n < N_base + 1; n++) {
            S1 = S1 + a[n] * pow(y, n + 1);
        }
        for (int m = 0; m < N_base; m++) {
            S2 = S2 + b[m] * pow(y, m + 1);
        }

        I_negative = (1 + S1) / (1 + S2);
        // возвращаем значение аппроксимации согласно формуле (14) в статье
        return boost::math::tgamma(k + 1)*y*pow(I_negative, k);
    }

    bmp_real integer::fd_one(const bmp_real x)
    {
        using namespace constants;
        const bmp_real I_1_0 = mp_pi*mp_pi / 12;
        const int k = 1;

        bmp_real I_1_minus_x = get_negative_value(-x, k, constants::a_k1, constants::b_k1);

        return x*x / 2 + 2 * I_1_0 - I_1_minus_x;
    }

    bmp_real integer::fd_two(const bmp_real x)
    {
        using namespace constants;
        const bmp_real I_1_0 = mp_pi*mp_pi / 12.0;
        const int k = 2;

        bmp_real I_2_minus_x = get_negative_value(-x, k, a_k2, b_k2);

        return x*x*x / 3 + 4 * x*I_1_0 + I_2_minus_x;
    }

    bmp_real integer::fd_3(const bmp_real x)
    {
        using namespace constants;
        const bmp_real I_1_0 = mp_pi*mp_pi / 12.0;
        const bmp_real I_3_0 = 7 * mp_pi*mp_pi*mp_pi*mp_pi / 120;
        const int k = 3;

        bmp_real I_3_minus_x = get_negative_value(-x, k, a_k3, b_k3);

        return x*x*x*x / 4 + 6 * x*x*I_1_0 + 2 * I_3_0 - I_3_minus_x;
    }

    bmp_real integer::fd_4(const bmp_real x)
    {
        using namespace constants;
        const bmp_real I_1_0 = mp_pi*mp_pi / 12.0;
        const bmp_real I_3_0 = 7 * mp_pi*mp_pi*mp_pi*mp_pi / 120;
        const int k = 4;

        bmp_real I_4_minus_x = get_negative_value(-x, k, a_k4, b_k4);
        return x*x*x*x*x / 5 + 8 * x*x*x*I_1_0 + 8 * x*I_3_0 + I_4_minus_x;
    }
}