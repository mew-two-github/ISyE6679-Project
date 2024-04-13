#include<bits/stdc++.h>
#include<codi.hpp>
using namespace std;

using EH = codi::EvaluationHelper;

template<typename Real>
void func(const Real*  x, Real* t) {
  t[0] = x[0] * x[0] * x[0];
  t[1] = x[0]*x[1]*30;
  t[2] = 5*x[1]*x[1];
}

auto hes = EH::createHessian(3, 2);
int main()
{
  double x[] = {1,2}, y[3];
  func(x,y);
  for(int i = 0; i <3; ++i)
  {
    cout<<y[i]<<" ";
  }
  cout<<endl;
  return 0;
}