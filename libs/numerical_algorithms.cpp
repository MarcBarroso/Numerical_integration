#include "numerical_algorithms.h"

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

using namespace std;

constexpr double EPS = 1e-11;
constexpr int MAX_STEPS = 100;

// Implementation of the 'extended trapezoidal rule'.
// Sec 4.1 in 'Numerical Recipies'.
// Error: ~1/n^2.
// INPUT: the function to be integrated in the interval (a, b), and the degree of the rule.
// OUTPUT: value of the integral.
double trapezoid_rule(const function<double(double)>& func, double a, double b, int n)
{
    if(n==1)
        return 0.5*(b-a)*(func(a)+func(b));

    double it = 1 << (n-1);
    double h = (b-a)/it;
    double sum = 0.0;
    sum+=0.5*func(a);
    for(int i=1; i < it; i++){
        sum += func(a+i*h);
    }
    sum+=0.5*func(b);

    return h*sum;
}

// Implementation of the polynomial interpolation of degree N-1 through N points.
// Sec 3.1 in 'Numerical Recipies'.
// INPUT: points and value of function at those points, and the position we want to know the value of the function at.
// OUTPUT: y: the value of the polynomial at x, and dy: an aproximation to the error of the approximation.
array<double, 2> polynomial_interpolation(const vector<double>& xa, const vector<double>& ya, double x)
{
    if(xa.size() != ya.size()){ cout << "ERROR in polynomial integration: size of vector xa != size of vector ya" << endl; }
    int n = xa.size(), ns=0;
    double y, dy;
    double c[xa.size()], d[xa.size()];
    double den;

    double dif = abs(x-xa[0]);
    double t_dif;

    for(int i=0;i<n;i++){
        t_dif = abs(x-xa[i]);
        if(t_dif < dif){
            ns = i;
            dif = t_dif;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }

    y = ya[ns];

    for(int m=1;m<=n;m++){
        for(int i=0;i<n-m;i++){
            double ho = xa[i]-x;
            double hp = xa[i+m]-x;
            double w = c[i+1]-d[i];
            den = ho-hp;
            if(den == 0.0){ cout << "ERROR in polynomial integration: two points with same x. " << endl; }
            den = w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        if(2*ns < (n-m)){
            dy = c[ns];
            y += dy;
        }else{
            if(ns > 0){
                ns--;
                dy = d[ns];
                y += dy;
            }
        }
    }

    array<double, 2> arr_sol = {y, dy};
    return arr_sol;
}

// Romberg's algorithm for computing the value of the integral numerically. It relies on successive evaluations
// of the trapezoidal rule with diferent number of points, and the extrapolate to get the result for h=0 (or points-> infty).
// The precision is controlled by the value of EPS. One can also change the parameters K and MAX_STEPS.
// Sec 4.3 in 'Numerical Recipies'.
// INTPUT: function we want to integrate in the interval (a, b).
// OUTPUT: value of the integral.
double romberg(const function<double(double)>& func, double a, double b)
{
    constexpr int K = 5;

    vector<double> s, h;
    array<double, 2> sol;

    h.push_back(1.0);
    for(int j=1;j<=MAX_STEPS;j++){
        s.push_back(trapezoid_rule(func, a, b, j));
        if(j >= K){
            sol = polynomial_interpolation(h, s, 0.0);
            if(abs(sol[1]) <= EPS*abs(sol[0])) 
                return sol[0];
        }
        h.push_back(0.25*h.back());
    }

    return -1;
}

// Finds the points and the weights for the Gauss-Legendre quadrature in the region (a, b) with n points.
// Sec 4.5 in 'Numerical Recipies'.
// INTPUT: function we want to integrate in the interval (a, b).
// OUTPUT: points of absicisses and weight needed for the G-L quad.
std::array<std::vector<double>, 2> gauss_legendre_quad(double a, double b, int n)
{
    vector<double> x, w;
    x.resize(n);
    w.resize(n);
    int m = (n+1)/2;

    double xm = 0.5*(a+b);
    double xl = 0.5*(b-a);
    double z;

    for(int i=0;i<m;i++){
        z = cos(M_PI*(i+1.0-0.25)/(n+0.5));
        double pp;
        double z1;
        do{
            double p0 = 1.0;
            double pn = 1.0;
            double po = 0.0;
            for(int j=0;j<n;j++){
                po = p0;
                p0 = pn;
                pn = ((2.0*j+1.0)*z*p0-j*po)/(j+1);
            }
            pp = n*(pn-p0)/(z*z-1.0);
            z1 = z;
            z=z1-pn/pp;
        }while(abs(z-z1) > EPS);

        x[i] = xm-xl*z;
        x[n-i-1] = xm+xl*z;
        w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
        w[n-i-1] = w[i];
    }

    std::array<std::vector<double>, 2> tmp = {x, w};
    return tmp;
}

// Computation of the integral by G-L quadrature. Uses the gauss_legendre_quad to obtain the points where the 
// function needs to be evaluated and multiplicative weights.
// TODO: implement an error function to estimate how well we do.
// INTPUT: function we want to integrate in the interval (a, b).
// OUTPUT: value of the integral.
double gauss_legendre_integration(const function<double(double)>& func, double a, double b, int n)
{
    auto quad = gauss_legendre_quad(a, b, n);
    double sum=0.0;

    for(int i=0; i<quad[0].size(); i++){
        sum+=func(quad[0][i])*quad[1][i];
    }

    return sum;
}
