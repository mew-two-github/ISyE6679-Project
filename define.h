#ifndef DEFINE_H_
#define DEFINE_H_

#include <iostream>
#include <vector>
#include <cmath>
#include<map>

using namespace std;

//function declarations

double f(vector<double> X, vector<double> a, vector<double> b, vector<double> c) {
		// Variable Definitions
		double result = 0;
		//double temp = 0;
		
		// Function Operations (reminder that Real Power are the first values of the X vector)
		// Formula 1
		for (int i=0; i < a.size(); i++) {
			//temp += a[i]*pow(X[i], 2) + b[i]*X[i] + c[i];
			//printf("a[%d]*X[%d]**2 + b[%d]*X[%d] + c[%d] = %f\n", i, i, i, i, i, temp);
			//result.push_back(temp);
			result += a[i]*pow(X[i], 2) + b[i]*X[i] + c[i];
		}
		
		// Return
		return result;
	}


void G(vector<double> X, vector<double> BusID_Gen, vector<double> BusID_REG, vector<int> FromBus, vector<int> ToBus,
				 double* PowerReg, vector<double> RealPowerDemand, vector<double> ReactPowerDemand, vector<double> G, vector<double> B,
				 vector<double> Real_P_G, vector<double> React_P_G, vector<double> Real_P_REG, vector<double> React_P_REG, double* y) 
{
		vector<double> result;
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
			y[i] = temp - temp2*X[2*NG + i];
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
			y[i] = temp - temp2*X[2*NG + i];
		}
}


vector<double> H(vector<double> X, vector<double> RealPowerMax, vector<double> RealPowerMin, 
				 vector<double> ReactPowerMax, vector<double> ReactPowerMin, vector<double> VoltMax, 
				 vector<double> VoltMin, vector<double> FromBus, vector<double> ToBus, vector<double> VoltAngMax, vector<double> VoltAngMin) 
{
		// Variable Definitions
		vector<double> result;
		double temp = 0;
		int Xplace = 0;
		
		// Function Operations
		// Formula 4
		for (int i=0; i < RealPowerMax.size(); i++) {
			//Max
			temp = X[Xplace] - RealPowerMax[i];
			result.push_back(temp);
			//Min
			temp = -1*X[Xplace] + RealPowerMin[i];
			result.push_back(temp);
			Xplace++;
		}
		// Formula 5
		for (int i=0; i < ReactPowerMax.size(); i++) {
			//Max
			temp = X[Xplace] - ReactPowerMax[i];
			result.push_back(temp);
			//Min
			temp = -1*X[Xplace] + ReactPowerMin[i];
			result.push_back(temp);
			Xplace++;
		}
		// Formula 7
		for (int i=0; i < VoltMax.size(); i++) {
			//Max
			temp = X[Xplace] - VoltMax[i];
			result.push_back(temp);
			//Min
			temp = -1*X[Xplace] + VoltMin[i];
			result.push_back(temp);
			Xplace++;
		}
		// Formula 9
		double theta_ij;
		for (int i=0; i < FromBus.size(); i++) {
			theta_ij = X[Xplace+FromBus[i]-1] - X[Xplace+ToBus[i]-1];
			//Max
			temp =  theta_ij - VoltAngMax[i];
			result.push_back(temp);
			//Min
			temp = -1*theta_ij + VoltAngMin[i];
			result.push_back(temp);

		}

		
		// Return
		return result;
	}

#endif