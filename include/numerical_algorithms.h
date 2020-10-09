#ifndef NUM_ALG
#define NUM_ALG

#include <array>
#include <functional>
#include <vector>

constexpr int MAX_STEPS = 100;

double trapezoid_rule(std::function<double(double)> func, double a, double b, int n);
std::array<double, 2> polynomial_interpolation(const std::vector<double>& xa, const std::vector<double>& ya, double x);
double romberg(std::function<double(double)> func, double a, double b);

#endif
