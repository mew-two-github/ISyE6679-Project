#include<bits/stdc++.h>
#include<codi.hpp>
using namespace std;


template<typename T>
// Function supports any data type T with which * is defined
T func(const T& x) {
  T t = x * x * x * x * x * x * x;
  return t * 3.0;
}

using t1s = codi::RealForwardGen<double>;
using t2s = codi::RealForwardGen<t1s>;
 
int main() {
 
  {
    t2s aFor = 2.0;
    // set all first order directions in order to get the 2. order derivative
    aFor.value().gradient() = 1.0;
    aFor.gradient().value() = 1.0;
 
    t2s cFor = func(aFor);
 
    std::cout << "t0s:   " << cFor.value().value() << std::endl;
    std::cout << "t1_1s: " << cFor.value().gradient() << std::endl;
    std::cout << "t1_2s: " << cFor.gradient().value() << std::endl;
    std::cout << "t2s:   " << cFor.gradient().gradient() << std::endl;

    // Expected Output 
    // f(2) = 3*2^7 = 384
    // f'(2) = 21*2^6 = 1344
    // f''(2) = 126*2^5 = 4032 
  }
  return 0;
}