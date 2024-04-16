#include<bits/stdc++.h>
#include<codi.hpp>
using namespace std;
#include "outputHelpers.hpp"
 
using EH = codi::EvaluationHelper;
// using Real = codi::RealForward;
template<typename Real>
void func(Real const*  x, Real* t) {
  // t = Real();
  t[0] = 4 * x[1] * x[0] + 2*x[1]*x[2] + x[2];
}
template<typename Real>
void funcwrap(std::vector<Real> const &x, std::vector<Real> &y) {
  size_t n = x.size() / 2;
  func(&x[0], &y[0]);
}
auto hes = EH::createHessian(3, 3);
int main()
{
  vector<double> x(3), y(1);
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;
  funcwrap(x,y);
  cout<<y[0];
  cout<<endl;
  EH::evalHessian(funcwrap<EH::HessianComputationType>, x, 3, hes);
  for(size_t j = 0; j < hes.getN(); j += 1) 
  {
    std::cout << "  ";
    for(size_t k = 0; k < hes.getN(); k += 1) 
    {
      if(k != 0) 
      {
        std::cout << ", ";
      }
      std::cout << hes(0, j, k);
    }
    std::cout << "\n";
  }
  double val = hes(0, 1, 0);
  cout<<val<<endl;
  return 0;
}