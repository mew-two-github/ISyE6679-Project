#ifndef DEFINE_H2_
#define DEFINE_H2_

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <codi.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "read_data.h"
#define ucmd ublas::compressed_matrix<double>

using Real = codi::RealForward;
using namespace std;
namespace ublas = boost::numeric::ublas;


// Bus data
vector<double> Volt, VoltMax, VoltMin, VoltAng, RealPowerDemand, ReactPowerDemand;
// Generator Data
vector<double> RealPower, ReactPower, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, BusID_Gen;
// Generator Cost Data
vector<double>  a, b, c;
// Branch Data
vector<double> VoltAngMax, VoltAngMin, Gvec, Bvec, FromBus, ToBus;
// declare temp vectors for REG Data
vector<double> Real_P_REG, React_P_REG,	BusID_REG;
int Nb, sizeY;


/*****Function headers*****/
ublas::compressed_matrix<double> RealVectorToDoubleUBlas(vector<vector<Real>> input, int rowSize, int colSize);
ublas::compressed_matrix<double> RealPointerToDoubleUBlasVec(Real* input, int rowSize);
ublas::compressed_matrix<double> uBLASNaturalLog(ublas::compressed_matrix<double> input);
Real* DoubleUBlasToRealVector(ucmd input, int rowSize, int colSize);

/******************************************** Fuction Declarations ********************************************/

void f(const Real* X, vector<double> a, vector<double> b, vector<double> c, Real* result ) {

	// Function Operations (reminder that Real Power are the first values of the X vector)
	result[0] = 0;
	// Formula 1
	for (int i=0; i < a.size(); i++) {
		//temp += a[i]*pow(X[i], 2) + b[i]*X[i] + c[i];
		//printf("a[%d]*X[%d]**2 + b[%d]*X[%d] + c[%d] = %f\n", i, i, i, i, i, temp);
		//result.push_back(temp);
		result[0] += a[i]*pow(X[i], 2) + b[i]*X[i] + c[i];
		//std::cout << result << ' ';
	}
}


void G_func(const Real* X, vector<double> BusID_Gen, vector<double> BusID_REG, vector<double> FromBus, vector<double> ToBus,
	vector<double> RealPowerDemand, vector<double> ReactPowerDemand, vector<double> G, vector<double> B, vector<double> Real_P_REG,
	vector<double> React_P_REG, Real* result) {
	
	// Variable Definitions
	Real temp = 0;
	int Xplace = 0;
	//vector<double> result;
	vector<vector<Real> > G_map, B_map, theta;
	int Nb = RealPowerDemand.size();
	
	// Fill theta, G and B values
	for( int i = 0; i < Nb; i++ )
	{
		vector<Real> v(Nb, 0.0);
		G_map.push_back(v);
		B_map.push_back(v);
		theta.push_back(v);
	}
	int Nbranches = FromBus.size();
	for (int i=0; i < Nbranches; i++) 
	{
		theta[FromBus[i]-1][ToBus[i]-1] = X[3*Nb+int(FromBus[i])-1] - X[3*Nb+int(ToBus[i])-1];
		G_map[FromBus[i]-1][ToBus[i]-1] = G[i];
		B_map[FromBus[i]-1][ToBus[i]-1] = B[i];
	}
	
	// Fill Power terms
	vector<Real> P_G(Nb, 0.0), P_R(Nb, 0.0), Q_G(Nb, 0.0), Q_R(Nb, 0.0);
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

	Real temp2 = 0;
	int counter = 0;
	
	// Real Power Balances
	for(int i = 0; i < Nb; i++ )
	{
		temp = P_G[i] + P_R[i] - RealPowerDemand[i];
		
		temp2 = 0;
		for(int j = 0; j< Nb; j++)
		{
			temp2 += X[2*NG + j]*( G_map[i][j]*cos(theta[i][j]) + B_map[i][j]*sin(theta[i][j]) );
		}
		result[counter++] = temp - temp2*X[2*NG + i];
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
		result[counter++] = temp - temp2*X[2*NG + i];
	}
}


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


/*void Lag(const Real* X , ucmd lambda, ucmd mu, ucmd gamma, ucmd Z, ucmd &res )
{
	 Real Cost[1], GX[2*Nb], HX[sizeY];

	 f(X, a, b, c, Cost);
	 ublas::compressed_matrix<double> CostuBLAS(1, 1);
	 CostuBLAS = RealPointerToDoubleUBlasVec(Cost, 1);
	
	 G_func(X, BusID_Gen, BusID_REG, FromBus, ToBus, RealPowerDemand, ReactPowerDemand, Gvec, Bvec, Real_P_REG, React_P_REG, GX);
	 ublas::compressed_matrix<double> GXuBLAS(2*Nb, 1);
	 GXuBLAS = RealPointerToDoubleUBlasVec(GX, 2*Nb);
	
	 H(X, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin, HX);
	 ublas::compressed_matrix<double> HXuBLAS(sizeY, 1);
	 HXuBLAS = RealPointerToDoubleUBlasVec(HX, sizeY);
	 
	 res = CostuBLAS + ublas::prod(trans(lambda), GXuBLAS) + ublas::prod(trans(mu), HXuBLAS + Z)
	       - ublas::prod(trans(gamma), uBLASNaturalLog(Z));	

}

void Real_Lag(const Real* X , ucmd lambda, ucmd mu, ucmd gamma, ucmd Z, ucmd res, Real* Y)
{
	Lag(X, lambda, mu, gamma, Z, res);
	cout<<res(0,0)<<endl;
	Y[0].value()= res(0,0);

}
*/

void Lag(const Real* X , ucmd lambda, ucmd mu, ucmd gamma, ucmd Z, Real* res )
{
	 int Nb = lambda.size1(), sizeY = mu.size1();

	 Real Cost[1], GX[2*Nb], HX[sizeY];

	 f(X, a, b, c, Cost);
	
	 G_func(X, BusID_Gen, BusID_REG, FromBus, ToBus, RealPowerDemand, ReactPowerDemand, Gvec, Bvec, Real_P_REG, React_P_REG, GX);

	 H(X, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin, HX);

	 size_t j = 0;
	 int counter = 0;
	 res[0] = Cost[0];
	 for(size_t i=0 ; i < lambda.size1();i++)
	 {
	 	res[0] += lambda(i,j)*GX[counter];
	 	counter++;
	 }
	 counter = 0;
	 for(size_t i=0 ; i < mu.size1();i++)
	 {
	 	res[0] += mu(i,j)*HX[i];
	 	res[0] -= gamma(i,j)*Z(i,j);
	 }
	 ++counter;

	 //res = CostuBLAS + ublas::prod(trans(lambda), GXuBLAS) + ublas::prod(trans(mu), HXuBLAS + Z)
	 //      - ublas::prod(trans(gamma), uBLASNaturalLog(Z));	

}

/*
vector<vector<Real>> LagX(Real* X, int XSize, ucmd lambda, ucmd mu, ucmd gamma,
 						 ucmd Z, ucmd res, vector<vector<Real>> result, int YSize = 1)
{
 	Real Y[1];
	for(int i = 0; i < XSize; ++i) 
 	{
 		// Step 1: Set tangent seeding
 		X[i].gradient() = 1.0;
 		// Step 2: Evaluate function
 		Real_Lag( X , lambda, mu, gamma, Z, res, Y);
 		for(int j = 0; j < YSize; ++j){

 			// Step 3: Access gradients
 			result[j][i] = Y[j].getGradient();
 			cout<<Y[j].value()<<" ";

 		}
 		// Step 4: Reset tangent seeding
 		X[i].gradient() = 0.0;
 	}	
 	return(result);
}


vector<vector<Real>> LagX(Real* X, int XSize, ucmd lambda, ucmd mu, ucmd gamma,
						 ucmd Z, vector<vector<Real>> result, int YSize = 1)
{
	Real Y[1];
	for(int i = 0; i < XSize; ++i) 
	{
		// Step 1: Set tangent seeding
		X[i].gradient() = 1.0;
		// Step 2: Evaluate function
		Lag( X , lambda, mu, gamma, Z, Y);
		for(int j = 0; j < YSize; ++j){

			// Step 3: Access gradients
			result[j][i] = Y[j].getGradient();
			cout<<Y[j].value()<<" ";

		}
		// Step 4: Reset tangent seeding
		X[i].gradient() = 0.0;
	}	
	return(result);
}
*/

/******************************************** Derivative Declarations ********************************************/

vector<vector<Real>> fX(Real* X, int XSize, Real* Y, int YSize, vector<double> a, vector<double> b, vector<double> c, 
						vector<vector<Real>> result)
{
	for(int i = 0; i < XSize; ++i) 
	{
		// Step 1: Set tangent seeding
		X[i].gradient() = 1.0;
		// Step 2: Evaluate function
		f(X, a, b, c, Y);
		for(int j = 0; j < YSize; ++j){

			// Step 3: Access gradients
			result[j][i] = Y[j].getGradient();

		}
		// Step 4: Reset tangent seeding
		X[i].gradient() = 0.0;
	}	
	return(result);
}


vector<vector<Real>> forwardModeFirstDerivativeH(Real* X, int XSize, Real* Y, int YSize, vector<double> RealPowerMax, vector<double> RealPowerMin, 
								vector<double> ReactPowerMax, vector<double> ReactPowerMin, vector<double> VoltMax, 
								vector<double> VoltMin, vector<double> FromBus, vector<double> ToBus, vector<double> VoltAngMax, 
								vector<double> VoltAngMin, vector<vector<Real>> result) {
								
	//printf("Number of Equations: %d, Number of Variables: %d\n", YSize, XSize);
	for(int i = 0; i < XSize; ++i) {
		// Step 1: Set tangent seeding
		X[i].gradient() = 1.0;
		// Step 2: Evaluate function
		H(X, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin, Y);
		for(int j = 0; j < YSize; ++j){

			// Step 3: Access gradients
			result[j][i] = Y[j].getGradient();

		}
		// Step 4: Reset tangent seeding
		X[i].gradient() = 0.0;
		//printf("%d ", j);
	}
	//fflush(stdout);
	
	return result;
}


vector<vector<Real>> GX_func(Real* X, int XSize, Real* Y, int YSize, vector<double> BusID_Gen, vector<double> BusID_REG, vector<double> FromBus, 
			vector<double> ToBus, vector<double> RealPowerDemand, vector<double> ReactPowerDemand, vector<double> Gvec, vector<double> Bvec,
			vector<double> Real_P_REG, vector<double> React_P_REG, vector<vector<Real>> result) {
			
	//printf("Number of Equations: %d, Number of Variables: %d\n", YSize, XSize);
	for(int i = 0; i < XSize; ++i) {
		// Step 1: Set tangent seeding
		X[i].gradient() = 1.0;
		// Step 2: Evaluate function
		G_func(X, BusID_Gen, BusID_REG, FromBus, ToBus, RealPowerDemand, ReactPowerDemand, Gvec, Bvec, Real_P_REG, React_P_REG, Y);
		for(int j = 0; j < YSize; ++j){

			// Step 3: Access gradients
			result[j][i] = Y[j].getGradient();

		}
		// Step 4: Reset tangent seeding
		X[i].gradient() = 0.0;
		//printf("%d ", j);
	}
	//fflush(stdout);
	
	return result;
}

vector<vector<Real>> forwardModeFirstDerivativeL(Real* X, int XSize, ucmd lambda, ucmd mu, ucmd gamma, ucmd Z, vector<vector<Real>> result, int YSize = 1)
{
	Real Y[1];
	for(int i = 0; i < XSize; ++i) 
	{
		// Step 1: Set tangent seeding
		X[i].gradient() = 1.0;
		// Step 2: Evaluate function
		Lag( X , lambda, mu, gamma, Z, Y);

		for(int j = 0; j < YSize; ++j){

			// Step 3: Access gradients
			result[j][i] = Y[j].getGradient();
			// cout << Y[j].value() << " ";

		}
		// Step 4: Reset tangent seeding
		X[i].gradient() = 0.0;

	}	

	return(result);
}

/******************************************* Useful Funct Declarations *******************************************/
ublas::compressed_matrix<double> RealVectorToDoubleUBlas(vector<vector<Real>> input, int rowSize, int colSize) {

	ublas::compressed_matrix<double> output(rowSize, colSize);

	for(int i = 0; i < rowSize; ++i) {
		for(int j = 0; j < colSize; ++j) {
			output(i,j) = input[i][j].value();
			//std::cout << input[i][j].value() << ' ';
		}
		//std::cout << std::endl;
	}
	
	return output;
}

ublas::compressed_matrix<double> RealPointerToDoubleUBlasVec(Real* input, int rowSize) {

	ublas::compressed_matrix<double> output(rowSize, 1);

	for(int i = 0; i < rowSize; ++i) {
		for(int j = 0; j < 1; ++j) {
			output(i,j) = input[i].value();
			//std::cout << input[i][j].value() << ' ';
		}
		//std::cout << std::endl;
	}
	
	return output;
}

ublas::compressed_matrix<double> uBLASNaturalLog(ublas::compressed_matrix<double> input) {

	//ublas::compressed_matrix<double> output(rowSize, colSize);

	for(int i = 0; i < input.size1(); ++i) {
		for(int j = 0; j < input.size2(); ++j) {
			input(i,j) = log(input(i,j));
			//std::cout << input[i][j].value() << ' ';
		}
		//std::cout << std::endl;
	}
	
	return input;
}


// void DoubleUBlasToRealVector(ucmd input, int rowSize, int colSize, Real* Y) {

// 	Real output[rowSize][colSize];

// 	for(int i = 0; i < rowSize; ++i) {
// 		for(int j = 0; j < colSize; ++j) {
// 			input[i][j].value() = output(i,j) ;
// 		}
// 	}
	
// 	return output;
// }









#endif
