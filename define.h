#ifndef DEFINE_H_

#define DEFINE_H_



#include <iostream>

#include <vector>

#include <cmath>

//function declarations



double f(std::vector<double> X, std::vector<double> a, std::vector<double> b, std::vector<double> c) {

		// Variable Definitions

		double result = 0;

		//double temp = 0;

		

		// Function Operations (reminder that Real Power are the first values of the X vector)

		// Formula 1

		for (int i=0; i < a.size(); i++) {

			//temp += a[i]*pow(X[i], 2) + b[i]*X[i] + c[i];

			//printf("a[%d]*X[%d]**2 + b[%d]*X[%d] + c[%d] = %f\n", i, i, i, i, i, temp);

			//result.insert(result.end(), temp);

			result += a[i]*pow(X[i], 2) + b[i]*X[i] + c[i];

		}

		

		// Return

		return result;

	}



/*

std::vector<double> G(std::vector<double> X) {

		// Variable Definitions

		std::vector<double> result;

		

		// Function Operations

		

		

		// Return

		return result;

	}

*/



std::vector<double> H(std::vector<double> X, std::vector<double> RealPowerMax, std::vector<double> RealPowerMin, std::vector<double> ReactPowerMax, std::vector<double> ReactPowerMin, std::vector<double> VoltMax, std::vector<double> VoltMin, std::vector<double> FromBus, std::vector<double> ToBus, std::vector<double> VoltAngMax, std::vector<double> VoltAngMin) {

		// Variable Definitions

		std::vector<double> result;

		double temp = 0;

		int Xplace = 0;

		

		// Function Operations

		// Formula 4

		for (int i=0; i < RealPowerMax.size(); i++) {

			//Max

			temp = X[Xplace] - RealPowerMax[i];

			result.insert(result.end(), temp);

			//Min

			temp = -1*X[Xplace] + RealPowerMin[i];

			result.insert(result.end(), temp);

			Xplace++;

		}

		// Formula 5

		for (int i=0; i < ReactPowerMax.size(); i++) {

			//Max

			temp = X[Xplace] - ReactPowerMax[i];

			result.insert(result.end(), temp);

			//Min

			temp = -1*X[Xplace] + ReactPowerMin[i];

			result.insert(result.end(), temp);

			Xplace++;

		}

		// Formula 7

		for (int i=0; i < VoltMax.size(); i++) {

			//Max

			temp = X[Xplace] - VoltMax[i];

			result.insert(result.end(), temp);

			//Min

			temp = -1*X[Xplace] + VoltMin[i];

			result.insert(result.end(), temp);

			Xplace++;

		}

		// Formula 9

		/*

		for (int i=0; i < FromBus.size(); i++) {

			//Max

			temp = X[Xplace] - VoltAngMax[i];

			result.insert(result.end(), temp);

			//Min

			temp = -1*X[Xplace] + VoltAngMin[i];

			result.insert(result.end(), temp);

			Xplace++;

		}

		*/

		

		// Return

		return result;

	}



#endif