#ifndef DEFINE_H2_
#define DEFINE_H2_

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <CoDiPack-master/include/codi.hpp>

using Real = codi::RealForward;

using namespace std;

//function declarations

void f(const Real* X, vector<double> a, vector<double> b, vector<double> c, Real* result) {

	// Function Operations (reminder that Real Power are the first values of the X vector)
	// Formula 1
	for (int i=0; i < a.size(); i++) {
		//temp += a[i]*pow(X[i], 2) + b[i]*X[i] + c[i];
		//printf("a[%d]*X[%d]**2 + b[%d]*X[%d] + c[%d] = %f\n", i, i, i, i, i, temp);
		//result.push_back(temp);
		result[0] += a[i]*pow(X[i], 2) + b[i]*X[i] + c[i];
		//std::cout << result << ' ';
	}
}

/*
void G(const Real* X, vector<double> BusID_Gen, vector<double> BusID_REG, vector<int> FromBus, vector<int> ToBus,
				 double* PowerReg, vector<double> RealPowerDemand, vector<double> ReactPowerDemand, vector<double> G,
				 vector<double> B, vector<double> Real_P_G, vector<double> React_P_G, vector<double> Real_P_REG,
				 vector<double> React_P_REG, Real* result) {
	// Variable Definitions
	Real temp = 0;
	int Xplace = 0, counter = 0;
	//vector<double> result;
	vector<vector<double> > G_map, B_map, theta;
	int Nb = RealPowerDemand.size();
	
	// Fill theta, G and B values
	for( int i = 0; i < Nb; i++ )
	{
		vector<double> v(Nb, 0.0);
		G_map.push_back(v);
		B_map.push_back(v);
		theta.push_back(v);
	}
	int Nbranches = FromBus.size();
	for (int i=0; i < Nbranches; i++) 
	{
		theta[FromBus[i]-1][ToBus[i]-1] = X[3*Nb+FromBus[i]-1] - X[3*Nb+ToBus[i]-1];
		G_map[FromBus[i]-1][ToBus[i]-1] = G[i];
		B_map[FromBus[i]-1][ToBus[i]-1] = B[i];
	}
	
	// Fill Power terms
	vector<double> P_G(Nb, 0.0), P_R(Nb, 0.0), Q_G(Nb, 0.0), Q_R(Nb, 0.0);
	int NG = BusID_Gen.size(), NR = BusID_REG.size();
	for( int i = 0; i < Nb; i++ )
	{
		for( int j = 0; j < NG; j++ )
		{
			if( i + 1 == BusID_Gen[j] )
			{	
				P_G[i] = X[j];
				Q_G[i] = X[j + NG];
				continue;
			}
		}
		for( int j = 0; j < NR; j++ )
		{
			if( i + 1 == BusID_REG[j] )
			{	
				P_R[i] = Real_P_REG[j];
				Q_R[i] = React_P_REG[j];
				continue;
			}
		}
	}

	double temp = 0, temp2 = 0;
	
	// Real Power Balances
	for(int i = 0; i < Nb; i++ )
	{
		temp = P_G[i] + P_R[i] - RealPowerDemand[i];
		
		temp2 = 0;
		for(int j = 0; j< Nb; j++)
		{
			temp2 += X[2*NG + j]*( G_map[i][j]*cos(theta[i][j]) + B_map[i][j]*sin(theta[i][j]) );
		}
		result[i] = temp - temp2*X[2*NG + i];
	}
	
	// Reactive Power Balances
	for(int i = 0; i < Nb; i++ )
	{
		temp = Q_G[NG + i] + Q_R[NG + i] - ReactPowerDemand[NG + i];
		
		temp2 = 0;
		for(int j = 0; j< Nb; j++)
		{
			temp2 += X[2*NG + j]*( G_map[i][j]*sin(theta[i][j]) - B_map[i][j]*cos(theta[i][j]) );
		}
		result[i] = temp - temp2*X[2*NG + i];
	}
}
*/


void H(const Real* X, vector<double> RealPowerMax, vector<double> RealPowerMin, vector<double> ReactPowerMax, vector<double> ReactPowerMin, 
                        vector<double> VoltMax, vector<double> VoltMin, vector<double> FromBus, vector<double> ToBus,
                        vector<double> VoltAngMax, vector<double> VoltAngMin, Real* result) {
	// Variable Definitions
	Real temp = 0, theta_ij = 0;
	int Xplace = 0, counter = 0;
		
	// Function Operations
	
	for (int i=0; i < RealPowerMax.size(); i++) {
		// Formula 4
		//Max
		temp = X[Xplace] - RealPowerMax[i];
		result[counter++] = temp;
		//Min
		temp = -1*X[Xplace] + RealPowerMin[i];
		result[counter++] = temp;
		Xplace++;
	}
	
	for (int i=0; i < ReactPowerMax.size(); i++) {
	
		// Formula 5
		//Max
		temp = X[Xplace] - ReactPowerMax[i];
		result[counter++] = temp;
		//Min
		temp = -1*X[Xplace] + ReactPowerMin[i];
		result[counter++] = temp;
		Xplace++;
	}
	
	for (int i=0; i < VoltMax.size(); i++) {
	
		// Formula 7
		//Max
		temp = X[Xplace] - VoltMax[i];
		result[counter++] = temp;
		//Min
		temp = -1*X[Xplace] + VoltMin[i];
		result[counter++] = temp;
		Xplace++;
	}

	for (int i=0; i < FromBus.size(); i++) {
		
		// Formula 9
		theta_ij = X[Xplace+int(FromBus[i])-1] - X[Xplace+int(ToBus[i])-1];
		//Max
		temp =  theta_ij - VoltAngMax[i];
		//temp =  (X[Xplace+int(FromBus[i])-1] - X[Xplace+int(ToBus[i])-1]) - VoltAngMax[i];
		result[counter++] = temp;
		//Min
		temp = -1*theta_ij+ VoltAngMin[i];
		//temp = -1*(X[Xplace+int(FromBus[i])-1] - X[Xplace+int(ToBus[i])-1]) + VoltAngMin[i];
		result[counter++] = temp;
	}

}

// Derivative Declarations

vector<vector<Real>> forwardModeFirstDerivativeH(Real* X, int XSize, Real* Y, int YSize, vector<double> RealPowerMax, vector<double> RealPowerMin, 
								vector<double> ReactPowerMax, vector<double> ReactPowerMin, vector<double> VoltMax, 
								vector<double> VoltMin, vector<double> FromBus, vector<double> ToBus, vector<double> VoltAngMax, 
								vector<double> VoltAngMin, vector<vector<Real>> result) {
	//printf("Number of Equations: %d, Number of Variables: %d\n", YSize, XSize);
	for(int i = 0; i < YSize; ++i) {
		for(int j = 0; j < XSize; ++j){
			// Step 1: Set tangent seeding
			X[j].gradient() = 1.0;
			
			// Step 2: Evaluate function
			H(X, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin, Y);
			
			// Step 3: Access gradients
			result[i][j] = Y[i].gradient();
			//std::cout << result[i][j] << ' ';
			
			// Step 4: Reset tangent seeding
			X[j].gradient() = 0.0;
			//printf("%d ", j);
		}
		//printf("%d\n", i);
		//std::cout << std::endl;
	}
	//fflush(stdout);
	
	return result;
}

// Basic Lin Alg



#endif