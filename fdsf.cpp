#include "fdsf.h"
#include "coefficients.h"
#include <boost/math/special_functions/gamma.hpp>
#include <iomanip>
#include <limits>

using namespace fdsf;

// �������� pi
const bmp_real PI = boost::math::constants::pi<bmp_real>();

// �������� Ik(0) ��� ������� k = 0,1,2,3
//const static bmp_real I_k_0[] = { log(bmp_real(2))
#if 0
bmp_real fdsf::get_T_max(bmp_real X, int k)
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
        // 3-� �������� ������ ���������� ��� ����������� Tmax
        if (i > 3) {
            break;
        }
        if (k == 0) {
            I_approximate = y;
        }
        else {
            I_approximate = fdsf::factorial(k)*y*pow(1 + a1*y, k);
        }
        // ������������ ���������� Tmax
        T = X - log(epsilon) + log(pow(T + 1, k) / I_approximate);
        i++;
    }

    return sqrt(T);
}
#endif

// ����� ����� ������� ��� ������������� ���������� y
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

// ��������� �������� ������� �� ������� k = 1, 2, 3 � ����� x ��� ��������
// �������� t. ������� ������������ � ���� t^k /(exp(x)+exp(t))
static bmp_real FermiDirakFunction(bmp_real t, bmp_real x, bmp_real k)
{
#if 0
    // new , ����� ��� ���������
    return 2 * pow(t, 2 * k + 1) / (1 + exp(t*t - x));
#endif
    return pow(t, k) / (boost::math::tgamma(bmp_real(k))*(exp(x) + exp(t)));
}

// ������� ������� ���� �������, ������� ����� ����������� ������������� � ����������� �� �������
static bmp_real get_negative_value(const bmp_real x, const bmp_real k)
{
    using namespace coefficients;
    const size_t N_base = 4;
    bmp_real S1 = 0, S2 = 0;
    bmp_real I_negative;
    bmp_real y = log(1 + exp(x));

    for (int n = 0; n < N_base + 1; n++) {
        S1 = S1 + a_k1[n]*pow(y, n + 1);
    }
    for (int m = 0; m < N_base; m++) {
        S2 = S2 + b_k1[m]*pow(y, m + 1);
    }

    I_negative = (1 + S1) / (1 + S2);
    return I_negative;
}

bmp_real fdsf::integer::fd_one(bmp_real x)
{
    const bmp_real I_1_0 = PI*PI / 12;
    const int k = 1;
  
    bmp_real I_1_minus_x = get_negative_value(-x, k);

    return x*x / 2 + 2*I_1_0 - I_1_minus_x;
}

bmp_real fdsf::integer::fd_two(bmp_real x)
{
    const bmp_real I_1_0 = PI*PI / 12.0;
    const int k = 2;

    bmp_real I_2_minus_x = get_negative_value(-x, k);

    return x*x*x / 3 + 4*x*I_1_0 + I_2_minus_x;
}

bmp_real fdsf::integer::fd_3(bmp_real x)
{
    const bmp_real I_1_0 = PI*PI / 12.0;
    const bmp_real I_3_0 = 7 * PI*PI*PI*PI / 120;
    const int k = 3;

    bmp_real I_3_minus_x = get_negative_value(-x, k);

    return x*x*x*x / 4 + 6*x*x*I_1_0 + 2*I_3_0 - I_3_minus_x;
}

bmp_real fdsf::integer::fd_4(const bmp_real x)
{
    const bmp_real I_1_0 = PI*PI / 12.0;
    const bmp_real I_3_0 = 7 * PI*PI*PI*PI / 120;
    const int k = 4;

    bmp_real I_4_minus_x = get_negative_value(-x, k);
    return x*x*x*x*x/5 + 8*x*x*x*I_1_0 + 8*x*I_3_0 + I_4_minus_x;
}
