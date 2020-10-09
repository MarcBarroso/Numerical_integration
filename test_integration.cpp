#include "numerical_algorithms.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double func1(double x){
    return exp(x)*sin(x);
    //return pow(x,4)*log(x+sqrt(x*x+1));
}

double func2(double x, double a){
    return sin(x)*cos(x);
}

class damper{
    function<double(double)> damp;
    double value;

    public:
        damper(function<double(double)> damp_) : damp(damp_){
            value = romberg(damp, 0, 1);
            cout << "VALUE: " << setprecision(20) << value << endl;
        }
};

int main(){
    auto start = chrono::steady_clock::now();
    using namespace std::placeholders;

    damper damp1(std::bind(func2, _1, 1.0));

    //auto integral = romberg(std::bind(func1, _1), 0, 2);
    //double real_val = 5.39689100903380;
    //cout << "Numerical value : " << setprecision(20) << integral  << endl;
    //cout << "Real value      : " << setprecision(20) << real_val << endl;
    //cout << "Error           : " <<  abs(integral-real_val) << endl;

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
