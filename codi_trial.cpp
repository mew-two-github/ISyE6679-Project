#include <iostream>
#include <vector>
#include <cmath>
#include <codi.hpp>
#include <iostream>
 
// template<typename Real>;
 
using Real = codi::RealReverse;

using namespace std;

void H(const Real* X, double* RealPowerMax, double* RealPowerMin, double* ReactPowerMax, double* ReactPowerMin, 
                        double* VoltMax, double* VoltMin, double* FromBus, double* ToBus, double* VoltAngMax, 
                        double* VoltAngMin, int n, int m, Real* result) 
    {
		// Variable Definitions
		Real temp = 0;
		int Xplace = 0, counter = 0;
		
		// Function Operations
		
		for (int i=0; i < n; i++) {
			// Formula 4
			//Max
			temp = X[Xplace] - RealPowerMax[i];
			result[counter++] = temp;
			//Min
			temp = -1*X[Xplace] + RealPowerMin[i];
			result[counter++] = temp;
			Xplace++;

		}
		for (int i=0; i < n; i++) {
		
			// Formula 5
			//Max
			temp = X[Xplace] - ReactPowerMax[i];
			result[counter++] = temp;
			//Min
			temp = -1*X[Xplace] + ReactPowerMin[i];
			result[counter++] = temp;
			Xplace++;
		}

	}


void reverseModeJacobianComputation() {
 
  using Tape = typename Real::Tape;
  int n = 2; // just power

  Real x[4] = {1,2,3,4};
  Real y[8]; // 2 constraints per term

  double RealPowerMax[2]= {5,50}, RealPowerMin[2] = {0,2}, ReactPowerMax[2] = {2,20}, ReactPowerMin[2] = {-1, -5}, VoltMax[2] ={0}, 
  VoltMin[2] ={0},FromBus[2]={0}, ToBus[2]= {0}, VoltAngMax[1] = {0}, VoltAngMin[1] = {0};

  int m=1;

  cout << endl;
 
  codi::Jacobian<double> jacobian(8,4);
 
  Tape& tape = Real::getTape();
  tape.setActive();
 
  // Step 1: Record the tape
  for(size_t i = 0; i < 4; ++i) {
    tape.registerInput(x[i]);
  }
 
  H( x, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin, n, m, y);
 
  tape.registerOutput(y[0]);
  tape.registerOutput(y[1]);
 
  tape.setPassive();
 
  // Step 2: Iterate over the output dimension
  for(size_t i = 0; i < 8; ++i) {
    y[i].gradient() = 1.0; // Step 3: Set the seeding for the i-th output variable
 
    tape.evaluate();
 
    // Step 4: Get the gradients from the inputs.
    for(size_t j = 0; j < 4; ++j) {
      jacobian(i, j) = x[j].getGradient();
    }
 
    tape.clearAdjoints(); // Step 5: Clear the adjoints
  }
 
  cout << "Reverse mode Jacobian:" << endl;
  cout << "f(1 .. 4) = (" << y[0] << ", " << y[1] << ")" << endl;
  cout << "df/dx (1 .. 4) = \n" << jacobian << endl;
 
  tape.reset(false);
}


int main() 
{
  reverseModeJacobianComputation();
  // int n = 2; // just power

  // const Real x[4] = {1,2,3,4};
  // Real y[8]; // 2 constraints per term

  // double RealPowerMax[2]= {5,50}, RealPowerMin[2] = {0,2}, ReactPowerMax[2] = {2,20}, ReactPowerMin[2] = {-1, -5}, VoltMax[2] ={0}, 
  // VoltMin[2] ={0},FromBus[2]={0}, ToBus[2]= {0}, VoltAngMax[1] = {0}, VoltAngMin[1] = {0};

  // int m=1;

  // cout << endl;
  // H( x, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin, y, n, m);
  // for(const Real &val : y)
  //   cout << "Output " << val << endl;
  // reverseModeJacobianComputation();
 
  return 0;
}