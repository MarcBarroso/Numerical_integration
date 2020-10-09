#ifndef NUM_ALG
#define NUM_ALG

#include <array>
#include <functional>
#include <vector>

double trapezoid_rule(const std::function<double(double)>& func, double a, double b, int n);
std::array<double, 2> polynomial_interpolation(const std::vector<double>& xa, const std::vector<double>& ya, double x);
double romberg(const std::function<double(double)>& func, double a, double b);
std::array<std::vector<double>, 2> gauss_legendre_quad(double a, double b, int n);
double gauss_legendre_integration(const std::function<double(double)>& func, double a, double b, int n);

#endif
