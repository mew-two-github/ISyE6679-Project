#include "define.h"
#include <vector> 
#include <iostream>
#include <cmath>

std::vector<double> f(std::vector<double> X, std::vector<double> a, std::vector<double> b, std::vector<double> c) {
		// Variable Definitions
		std::vector<double> result;
		double temp;
		
		// Function Operations
		for (int i=0; i < X.size(); i++) {
			temp = a[i]*pow(X[i], 2) + b[i]*pow(X[i], 1) + c[i];
			result.insert(result.end(), temp);
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

std::vector<double> H(std::vector<double> X) {
		// Variable Definitions
		std::vector<double> result;
		
		// Function Operations
		
		
		// Return
		return result;
	}
*/
