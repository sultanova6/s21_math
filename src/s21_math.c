#include "s21_math.h"

int s21_abs(int x) {
    return x < 0 ? x * -1: x;
}

long double s21_fabs(double x) {
    return x < 0 ? x * -1: x;
}

long double s21_exp(double x) {
    long double r = 1, t = 1,  i = 1;
    int f = x < 0 ? 1: 0;
    x = f == 1 ? x * -1: x;
    while (s21_fabs(r) > s21_EPS) {
        r *= x / i;
        i += 1;
        t += r;
        if (t > s21_MAX_double) {
            t = s21_INF;
            break;
        }
    }
    t = f == 1 ? t > s21_MAX_double ? 0 : 1. / t : t;
    return t > s21_MAX_double ? s21_INF : t;
}

long double s21_fmod(double x, double y) {
    return (x && y) || y != 0 ? (x / y - (int)(x / y)) * y : s21_NAN;
}

double re_mod(double x) {
    return (x > 2 * s21_PI || x < -2 * s21_PI) ?  s21_fmod(x, 2 * s21_PI): x;
}

long double s21_sin(double x) {
    x = re_mod(x);
    long double r = 0, p = x, f = 1;
    for (int i = 0; s21_fabs(p / f) > s21_EPS; i++) {
        r += p / f;
        p *= -1 * x * x;
        f *= (2 * (i + 1)) * (2 * (i + 1) + 1);
    }
    return r;
}

long double s21_cos(double x) {
    x = re_mod(x);
    return s21_sin(s21_PI / 2 + x);
}

long double s21_tan(double x) {
    return s21_sin(x) / s21_cos(x);
}

long double s21_asin(double x) {
    double t = x, r = x;
    if (x < -1 || x > 1) {
        r = s21_NAN;
    } else if (x == -1 || x == 1) {
        r = s21_PI / (2 * x);
    } else {
        for (long double c = 1; c < 10000000; c++) {
            t *= ((x * x) * (2 * c - 1) * (2 * c - 1)) /
                   ((2 * c) * (2 * c + 1));
            r += t;
        }
    }
    return r;
}

long double s21_acos(double x) {
    return (x < -1 || x > 1) ? s21_NAN :
    (x == -1 || x == 1) ? s21_PI * (1 - x) / 2 :
    s21_PI / 2. - s21_asin(x);
}

long double p_atan(double x) {
    long double r = x, t = x, i = 1;
    while (s21_fabs(r) > s21_EPS) {
        r = -1 * r * x * x * (2 * i - 1) / (2 * i + 1);
        i += 1;
        t += r;
    }
    return t;
}

long double s21_atan(double x) {
    return (x < 1 && x > -1 && x != 0) ? p_atan(x):
    x == 0 ? 0:
    x == 1 ? s21_PI / 4:
    x == -1 ? -s21_PI / 4:
    x > 1 ? s21_PI / 2 - p_atan(1 / x):
    x < -1 ? -s21_PI / 2 - p_atan(1 / x):
    0;
}

void tr_e(double x, opt *param) {
    long double i = 1;
    int e = 0;
    param->s_n = x < 0 ? -1 : 1;
    x *= param->s_n;
    while (x >= 10) {
        x /= 10.;
        i *= 10;
        e++;
    }
    param->m_a = x;
    param->pow = i;
    param->E = e;
}

long double p_log(double x) {
    x--;
    long double r = x, t = x;
    long double i = 2;
    while (s21_fabs(r) > s21_EPS) {
        r *= -x * (i - 1) / i;
        i += 1;
        t += r;
    }
    return t;
}

long double e_log(double x) {
    opt param;
    tr_e(x, &param);
    x = param.m_a * param.s_n / 10;
    return p_log(x) + (param.E + 1) * s21_ln10;
}

long double s21_log(double x) {
    return (x > 0 && x < 2) ? p_log(x) :
    (x < 0) ? s21_NAN :
    (x == 0) ? -s21_INF :
    e_log(x);
}

long double s21_pow(double b, double e) {
    long double num;
    if (b < 0) {
        if ((long int)e == e) {
            if (e > 0) {
                num = b;
                for (long int i = 0; i < (long int)e - 1; i++) {
                    num *= b;
                }
            } else if (e == 0) {
                num = 1;
            } else {
                num = 1 / b;
                for (long int i = 0; i < (long int)e * (-1) - 1; i++) {
                    num /= b;
                }
            }
        } else {
            if (e == -s21_INF || e == s21_INF) {
                if (b * (-1) < 1) {
                    num = 0;
                } else if (b * (-1) == 1) {
                    num = 1;
                } else {
                    if (e == -s21_INF) {
                        num = 0;
                    } else {
                        num = s21_INF;
                    }
                }
            } else {
                num = -s21_NAN;
            }
        }
    } else if (b == 0) {
        num = e == 0 ? 1: 0;
    } else if (b == 1) {
        num = 1;
    } else {
        if ((long int)e == e) {
            if (e > 0) {
                num = b;
                for (long int i = 0; i < (long int)e - 1; i++) {
                    num *= b;
                }
            } else if (e == 0) {
                num = 1;
            } else {
                num = 1 / b;
                for (long int i = 0; i < (long int)e * (-1) - 1; i++) {
                    num /= b;
                }
            }
        } else {
            num = s21_exp(e * (double)s21_log(b));
        }
    }
    return num;
}

long double w_sqrt(long double r, long double t, double x) {
    while (s21_fabs(r - t) > s21_EPS) {
        t = r;
        r = (t + x / t) / 2;
    }
    return r;
}

long double s21_sqrt(double x) {
    long double r = 4, t = 0;
    return (x < 0) ? -s21_NAN : w_sqrt(r, t, x);
}

long double s21_floor(double x) {
    if (x < 0) {
        x -= 0.9999999;
    }
    int answer = (int)x;
    return answer;
}

long double s21_ceil(double x) {
    int xi = (int)x;
    double temp = x - xi;
    if (temp > 0.0) x += 0.9999999;
    return (int)x;
}
