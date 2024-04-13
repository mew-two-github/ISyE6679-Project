#include<bits/stdc++.h>
#include<codi.hpp>
using namespace std;


using t1s = codi::RealForwardGen<double>;
using t2s = codi::RealForwardGen<t1s>;
// Function supports any data type T with which * is defined
t2s* func(const t2s*  x) {
  t2s t[] = { x[0] * x[0] * x[0], x[0]*x[1]*30, 5*x[1]*x[1] };
  return t ;
}

 
int main() {
 
  {
    t2s* aFor[] = {2.0,3.0};
    // set all first order directions in order to get the 2. order derivative
  //     t2s cFor[3] = func(aFor);
    for(int i=0;i<2;++i)
    {
        aFor[i].value().gradient() = 1.0;
        aFor[i].gradient().value() = 1.0;
        t2s cFor[3] = func(aFor);
        for(int j = 0; j < 3; ++j)
        {
          cFor[j] = 
        }
    }
 
    
 
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