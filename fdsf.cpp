#include "fdsf.h"
#include <iomanip>
#include <limits>

namespace fdsf {

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

    bmp_real integer::fermi_dirak_k1(const bmp_real x)
    {
        using namespace constants;
        const int k = 1;
        const mp_coefficients a = { 0.2715113138214362780964488,
                                    0.0562661238060587633169245,
                                    0.0067420740469345689743873,
                                    0.0005169505155333205859985,
                                    0.0000194771836765773190602 };

        const mp_coefficients b = { 0.0215113138214352840651584,
                                    0.0231105175729721417901084,
                                    0.0003669081577365413477999,
                                    0.0000610424408732720110769 };

        bmp_real I_1_minus_x = get_negative_value(-x, k, a, b);

        return x*x / 2 + 2 * I_1_0 - I_1_minus_x;
    }

    bmp_real integer::fermi_dirak_k2(const bmp_real x)
    {
        using namespace constants;
        const int k = 2;
        const mp_coefficients a = { 0.2263816364340698560028783,
                                    0.0533684335574798857246766,
                                    0.0062904756340795211604491,
                                    0.0005023228274452983506998,
                                    0.0000189379675088061004880 };

        const mp_coefficients b = { 0.0388816364340691133155655,
                                    0.0243043998742774445085992,
                                    0.0006290985326433190105734,
                                    0.0000657018161945458806177 };

        bmp_real I_2_minus_x = get_negative_value(-x, k, a, b);

        return x*x*x / 3 + 4 * x*I_1_0 + I_2_minus_x;
    }

    bmp_real integer::fermi_dirak_k3(const bmp_real x)
    {
        using namespace constants;

        const int k = 3;
        const mp_coefficients a = { 0.1583482145380455955096383,
                                    0.0460645149909308107878344,
                                    0.0048861379108841469134267,
                                    0.0004336733305971515517559,
                                    0.0000173435613795895152436 };

        const mp_coefficients b = { 0.0125148812047107612191739,
                                    0.0266693407000929631393759,
                                    0.0003285431094547362504004,
                                    0.0000820910787890062715299 };

        bmp_real I_3_minus_x = get_negative_value(-x, k, a, b);

        return x*x*x*x / 4 + 6 * x*x*I_1_0 + 2 * I_3_0 - I_3_minus_x;
    }

    bmp_real integer::fermi_dirak_k4(const bmp_real x)
    {
        using namespace constants;

        const int k = 4;

        const mp_coefficients a = { 0.0560148791230902149024568,
                                    0.0351117957891800867706741,
                                    0.0021834386943672331415760,
                                    0.0002464861525522946634693,
                                    0.0000092228177886669241259 };

        const mp_coefficients b = { -0.0611726208769112866900252,
                                     0.0279968542816146833953639,
                                    -0.0007512148294307540141223,
                                     0.0000860680747142919882956 };

        bmp_real I_4_minus_x = get_negative_value(-x, k, a, b);
        return x*x*x*x*x / 5 + 8 * x*x*x*I_1_0 + 8 * x*I_3_0 + I_4_minus_x;
    }
}