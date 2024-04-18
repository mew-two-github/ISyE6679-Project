#ifndef SECOND_DERIVS_H
#define SECOND_DERIVS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <codi.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "read_data.h"
#define ucmd ublas::compressed_matrix<double>

using namespace std;
namespace ublas = boost::numeric::ublas;




/******************************************** Fuction Declarations ********************************************/
template<typename T>
void Sf(const T* X, T* result ) {

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

template<typename T>
void SG_func(const T* X, T* result) {
	
	// Variable Definitions
	double temp = 0;
	int Xplace = 0;
	//vector<double> result;
	vector<vector<T> > G_map, B_map, theta;
	int Nb = RealPowerDemand.size();
	
	// Fill theta, G and B values
	for( int i = 0; i < Nb; i++ )
	{
		vector<T> v(Nb, 0.0);
		G_map.push_back(v);
		B_map.push_back(v);
		theta.push_back(v);
	}
	int Nbranches = FromBus.size();
	for (int i=0; i < Nbranches; i++) 
	{
		theta[FromBus[i]-1][ToBus[i]-1] = X[3*Nb+int(FromBus[i])-1] - X[3*Nb+int(ToBus[i])-1];
		G_map[FromBus[i]-1][ToBus[i]-1] = Gvec[i];
		B_map[FromBus[i]-1][ToBus[i]-1] = Bvec[i];
	}
	
	// Fill Power terms
	vector<T> P_G(Nb, 0.0), P_R(Nb, 0.0), Q_G(Nb, 0.0), Q_R(Nb, 0.0);
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

	double temp2 = 0;
	int counter = 0;
	
	// Real Power Balances
	for(int i = 0; i < Nb; i++ )
	{
		result[counter] = P_G[i] + P_R[i] - RealPowerDemand[i];
		
		// temp2 = 0;
		for(int j = 0; j< Nb; j++)
		{
			result[counter] -= X[2*NG + j]*( G_map[i][j]*cos(theta[i][j]) + B_map[i][j]*sin(theta[i][j]) );
		}
		counter++;
	}
	
	// Reactive Power Balances
	for(int i = 0; i < Nb; i++ )
	{
		result[counter] = Q_G[NG + i] + Q_R[NG + i] - ReactPowerDemand[NG + i];
		
		// temp2 = 0;
		for(int j = 0; j< Nb; j++)
		{
			result[counter] += X[2*NG + j]*( G_map[i][j]*sin(theta[i][j]) - B_map[i][j]*cos(theta[i][j]) );
		}
		counter++;
	}
}

template<typename T>
void SH(const T* X, T* result) {
                        
	// Variable Definitions
	T temp = 0, theta_ij = 0;
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


template<typename T>
void SLag(T const* X , T const* lambda, T const* mu, T const* gamma, T const* Z, T* res )
{

	 T Cost[1], GX[2*Nb], HX[sizeY];

	 Sf(X, Cost);
	
	 SG_func(X, GX);

	 SH(X, HX);

	 res[0] = Cost[0];
	 for(int i = 0; i < sizeX; ++i)
	 	cout<<X[i]<<" ";
	 for(int i=0 ; i < 2*Nb;i++)
	 {
	 	res[0] += lambda[i]*GX[i];
	 }
	 for(int i=0 ; i < sizeY;i++)
	 {
	 	res[0] += mu[i]*HX[i];
	 	res[0] -= gamma[i]*Z[i];
	 }

	 //res = CostuBLAS + ublas::prod(trans(lambda), GXuBLAS) + ublas::prod(trans(mu), HXuBLAS + Z)
	 //      - ublas::prod(trans(gamma), uBLASNaturalLog(Z));	

}

template<typename T>
void SLagWrap(vector<T> const &Xvec , vector<T> &res )
{
	T lambda[2*Nb], mu[sizeY], gamma[sizeY], Xarr[sizeX], Z[sizeY];
	int counter = 0;

	for(int i = 0;i < sizeX; i++ )
		Xarr[i] = Xvec[counter++];

	for(int i = 0; i < 2*Nb; i++ )
		lambda[i] = Xvec[counter++];

	for(int i = 0; i < sizeY; i++)
	{
		mu[i] = Xvec[counter];
		gamma[i] = Xvec[counter+sizeY];
		Z[i] = Xvec[counter+sizeY*2];
		counter++;
	}
	T y[1];
	SLag( Xarr, lambda, mu, gamma, Z, y );
	res[0] = y[0];
}

using EH = codi::EvaluationHelper;

auto hes = EH::createHessian(1, sizeX);

void DoubleUBlasVecToDoublePointer(ublas::compressed_matrix<double> input, double* output) {

	for(int i = 0; i < input.size1(); ++i) {
		for(int j = 0; j < input.size2(); ++j) {
			output[i] = input(i,j);
			//std::cout << input[i][j].value() << ' ';
		}
		//std::cout << std::endl;
	}
}

void second_deriv(Real* X , ucmd lambda, ucmd mu, ucmd gamma, ucmd Z)
{
	vector<double> y(1);
	vector<double> Xvec(sizeX+2*Nb+sizeY*3);
	double lambda_arr[2*Nb], mu_arr[sizeY], gamma_arr[sizeY], Z_arr[sizeY];
	DoubleUBlasVecToDoublePointer( lambda,lambda_arr );
	DoubleUBlasVecToDoublePointer( mu,mu_arr );
	DoubleUBlasVecToDoublePointer( gamma,gamma_arr );
	DoubleUBlasVecToDoublePointer( Z,Z_arr );

	int counter = 0;
	for(int i = 0;i < sizeX; i++ )
	{
		Xvec[counter++] = X[i].value();
	}
		

	for(int i = 0; i < 2*Nb; i++ )
		Xvec[counter++] = lambda_arr[i];

	for(int i = 0; i < sizeY; i++)
	{
		Xvec[counter] = mu_arr[i];
		Xvec[counter+sizeY] = gamma_arr[i];
		Xvec[counter+sizeY*2] = Z_arr[i];
		counter++;
	}
	SLagWrap(Xvec,y);
	cout<<y[0];
	cout<<endl;
	// EH::evalHessian(SLagWrap<EH::HessianComputationType>, Xvec, 1, hes);
	// for(size_t j = 0; j < hes.getN(); j += 1) 
	// {
	// std::cout << "  ";
	// for(size_t k = 0; k < hes.getN(); k += 1) 
	// {
	//   if(k != 0) 
	//   {
	//     std::cout << ", ";
	//   }
	//   std::cout << hes(0, j, k);
	// }
	// std::cout << "\n";
	// }
}





#endif
