#include <iostream>
#include <cmath>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace std;

const double eps = 0.001;

double f1(double x) { return pow(2, x) + 1; }

double f1_d(double x) { return log(2) * pow(2, x); }

double f1_dd(double x) { return log(2) * log(2) * pow(2, x); }

double f2(double x) { return pow(x, 5); }

double f2_d(double x) { return 5 * pow(x, 4); }

double f2_dd(double x) { return 20 * pow(x, 3); }

double f3(double x) { return (1 - x) / 3; }

double f3_d(double x) { return -x / 3; }

double f3_dd(double x) { return -1.0 / 3; }

double dihotomy_root(double a, double b, double(f)(double), double(g)(double))
{
    double fga = f(a) - g(a), fgb = f(b) - g(b), x;
    if (fga * fgb > 0)
        throw runtime_error("(f - g)(a) * (f - g)(b) > 0");
    while (abs(a - b) > eps)
    {
        x = (a + b) / 2;
        if ((f(x) - g(x)) * fga > 0)
            a = x;
        else
            b = x;
    }
    return (a + b) / 2;
}

double combined_root(double a, double b, double(f)(double), double(g)(double), double(f_d)(double), double(g_d)(double), double(f_dd)(double), double(g_dd)(double))
{
    double fga = f(a) - g(a), fgb = f(b) - g(b), fga_d, fgb_d, fga_dd, fgb_dd;
    if (fga * fgb > 0)
        throw runtime_error("(f - g)(a) * (f - g)(b) > 0");
    while (abs(a - b) > eps)
    {
        fga = f(a) - g(a);
        fgb = f(b) - g(b);
        fga_d = f_d(a) - g_d(a);
        fgb_d = f_d(b) - g_d(b);
        fga_dd = f_dd(a) - g_dd(a);
        fgb_dd = f_dd(b) - g_dd(b);

        if (fga_d * fga_dd < 0)
            a -= fga * (a - b) / (fga - fgb);
        else
            a -= fga / fga_d;
        if (fgb_d * fgb_dd < 0)
            b -= fgb * (b - a) / (fgb - fga);
        else
            b -= fgb / fgb_d;
    }
    return (a + b) / 2;
}

double l_rect_integral(double a, double b, double f(double))
{
    double h = (b - a) / 1000, s = 0;
    for (int i = 0; i < 1000; i++)
        s += f(a + i * h) * h;
    return s;
}

double trap_integral(double a, double b, double(f)(double))
{
    double h = (b - a) / 1000, s = 0, x1 = a, x2 = a;
    for (int i = 1; i <= 1000; i++)
    {
        x2 += h;
        s += f((x1 + x2) / 2) * h;
        x1 = x2;
    }
    return s;
}

void plot_the_whole_thing(double f1f2, double f1f3, double f2f3)
{
    plt::figure_size(1200, 780);
    vector<double> x, y1, y2, y3, crv;
    double mn = min(min(f1f2, f2f3), f1f3), mx = max(max(f1f2, f1f3), f2f3);
    for (double i = mn - 0.7; i <= mx + 0.7; i += 0.01)
    {
        x.push_back(i); y1.push_back(f1(i)); y2.push_back(f2(i)); y3.push_back(f3(i));
        if ((i <= f1f2 + f2f3 + f1f3 - mn - mx) && (i >= mn))
            crv.push_back(f3(i));
        else if ((i > f1f2 + f2f3 + f1f3 - mn - mx) && (i <= mx))
            crv.push_back(f2(i));
        else crv.push_back(f1(i));
    }

    double f1r = f1(f1f2), f2r = f2(f2f3), f3r = f3(f1f3);
    plt::named_plot("f1", x, y1);
    plt::named_plot("f2", x, y2);
    plt::named_plot("f3", x, y3);
    plt::ylim(min(min(f1r, f2r), f3r) - 0.7, max(max(f1r, f3r), f2r) + 0.7);
    std::map<std::string, std::string> keywds;
    keywds["facecolor"] = "lightsalmon";
    keywds["edgecolor"] = "none";
    plt::fill_between(x, y1, crv, keywds);
    plt::title("Area within three curves");
    plt::legend();
    plt::save("./plot.png");
    plt::show();
}

int main()
{
    double f1f2 = dihotomy_root(-10, 10, f1, f2),
        f1f3 = dihotomy_root(-10, 10, f1, f3),
        f2f3 = dihotomy_root(-10, 10, f2, f3),
        f1_int, f2_int, f3_int;

    std::cout << "\n\n"
        << "The intersection point of f1 and f2 is "
        << f1f2 << "\n"
        << "The intersection point of f1 and f3 is "
        << f1f3 << "\n"
        << "The intersection point of f2 and f3 is "
        << f2f3 << "\n";

    f1_int = (f1f2 < f1f3) ? trap_integral(f1f2, f1f3, f1) : trap_integral(f1f3, f1f2, f1);
    f2_int = (f1f2 < f2f3) ? trap_integral(f1f2, f2f3, f2) : trap_integral(f2f3, f1f2, f2);
    f3_int = (f1f3 < f2f3) ? trap_integral(f1f3, f2f3, f3) : trap_integral(f2f3, f1f3, f3);

    double area = 2 * max(max(f1_int, f2_int), f3_int) - (f1_int + f2_int + f3_int);
    std::cout << "\n"
        << "The area within the three curves equals to "
        << area << "\n";
    
    plot_the_whole_thing(f1f2, f1f3, f2f3);

    return 0;
}