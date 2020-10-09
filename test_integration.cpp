#include "numerical_algorithms.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double func1(double x){
    return exp(-x*x);
}

double func2(double x, double a){
    return cos(x/a)*cos(x/M_PI);
}

double func3(double x, double a, double b, double c){
    return (x-a)/((x-b)+c);
}

class damper{
    function<double(double)> damp;
    double value;

    public:
        damper(function<double(double)> damp_) : damp(damp_){
            value = romberg(damp, 0, 1);
            cout << "Damping function ready! Value: " << setprecision(20) << value << endl;
        }

        void print_dump(double x){
            cout << "Value of the damping function @ x = " << x << " : " << damp(x) << endl;
        }


};

int main(){
    auto start = chrono::steady_clock::now();
    using namespace std::placeholders;

    cout << "----------------------------" << endl;

    damper damp1(std::bind(func2, _1, 2.0));
    damp1.print_dump(0.5);

    cout << "----------------------------" << endl;

    damper damp2(func1);
    damp2.print_dump(0.5);

    cout << "----------------------------" << endl;

    cout << "Comparison of Romberg's method and Gauss-Legendre quadrature" << endl;
    cout << "Romberg: " << setprecision(20) << romberg(std::bind(func3, _1, 5, 1, 2), 0, 1) << endl;
    cout << "Gauss-Legendre (15 points): " << setprecision(20) << gauss_legendre_integration(std::bind(func3, _1, 5, 1, 2), 0, 1, 15) << endl;

    chrono::duration<double> elapsed_seconds = chrono::steady_clock::now()-start;
    cout << "----------------------------" << endl;
    cout << "Elapsed time: " << setprecision(5) << elapsed_seconds.count() << "s" << endl;
    return 0;
}
