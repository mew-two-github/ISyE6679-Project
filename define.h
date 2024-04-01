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


vector<double> G(vector<double> X, vector<double> vector<double> G, vector<double> B) {
		// Variable Definitions
		vector<double> result;
		map<int,int> m;

		// Real Power Balances
		for(int i = 0; i < RealPowerMax.size(); i++ )
		{
			temp 
			if (m.find("f") == m.end()) {
			// not found
			} 
			else {
			// found
			}
		}
		
		// Reactive Power Balances
		for(int i = 0; i < RealPowerMax.size(); i++ )
		{
			temp 
			if (m.find("f") == m.end()) {
			// not found
			} 
			else {
			// found
			}
		}

		return result;
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
		for (int i=0; i < FromBus.size(); i++) {
			theta_ij = X[Xplace+FromBus[i]-1] - X[Xplace+ToBus[i]-1]
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


