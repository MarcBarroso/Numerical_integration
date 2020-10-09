#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <array>
#include <vector>

using namespace std;

double func1(double x){
    return exp(x)*sin(x);
    //return pow(x,4)*log(x+sqrt(x*x+1));
}

double func2(double x, double a){
    return sin(x)*cos(x);
}

double trapezoid_rule(function<double(double)> func, double a, double b, int n)
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

double romberg(function<double(double)> func, double a, double b)
{
    constexpr int K = 5;
    constexpr double eps = 1e-20;

    vector<double> s, h;
    array<double, 2> sol;

    h.push_back(1.0);
    for(int j=1;j<=100;j++){
        s.push_back(trapezoid_rule(func, a, b, j));
        if(j >= K){
            sol = polynomial_interpolation(h, s, 0.0);
            //cout << "s, j, sol[0], sol[1]: " << s.back() << " " << j << " " << sol[0] << " " << sol[1] << endl;
            if(abs(sol[1]) <= eps*abs(sol[0])) 
                return sol[0];
        }
        h.push_back(0.25*h.back());
    }

    return -1;
}

int main(){
    using namespace std::placeholders;

    auto start = chrono::steady_clock::now();

    auto integral = romberg(std::bind(func1, _1), 0, 2);
    double real_val = 5.39689100903380;
    cout << "Numerical value : " << setprecision(20) << integral  << endl;
    cout << "Real value      : " << setprecision(20) << real_val << endl;
    cout << "Error           : " <<  abs(integral-real_val) << endl;

    //vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    //vector<double> y = {func1(1), func1(2), func1(3), func1(4), func1(5), func1(6), func1(7), func1(8), func1(9), func1(10)};

    //array<double, 2> sol = polynomial_interpolation(x, y, 5.5);
    //cout << "The interpolation yields: " << sol[0] << endl;
    //cout << "The real value is: " << func1(5.5) << endl;
    //cout << "Real error: " << abs(sol[0]-func1(5.5)) << endl;
    //cout << "Estimated error: " << abs(sol[1]) << endl;

    chrono::duration<double> elapsed_seconds = chrono::steady_clock::now()-start;
    cout << "------------------------------------" << endl;
    cout << "Elapsed time: " << elapsed_seconds.count() << "s" << endl;
    return 0;
}
