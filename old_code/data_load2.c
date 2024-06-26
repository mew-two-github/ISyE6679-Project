#include "csv.h"
#include "define2.h"
#include <iostream> 
#include <vector>
#include "read_data.h"
#include <codi.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

using Real = codi::RealForward;
namespace ublas = boost::numeric::ublas;

void readREG(char filename[]){
  io::CSVReader<7> in(filename);
  in.read_header(io::ignore_extra_column, "Bus ID", "Real Power Output", "React Power Output", "Max Real Power Output", "Min Real Power Output", "Max React Power Output", "Min React Power Output");
  double Bus; double Po; double Qo; double Pmax; double Pmin; double Qmax; double Qmin;
  while(in.read_row(Bus, Po, Qo, Pmax, Pmin, Qmax, Qmin)){
    printf("Bus ID: %d, Real Power Output: %f, React Power Output: %f,\n\tMax Real Power Output: %f, Min Real Power Output: %f, Max React Power Output: %f, Max React Power Output: %f\n", int(Bus), Po, Qo, Pmax, Pmin, Qmax, Qmin);
  }
}

void readBranch(char filename[]){
  io::CSVReader<6> in(filename);
  in.read_header(io::ignore_extra_column, "From Bus", "To Bus", "Max Angle Difference", "Min Angle Difference", "G" , "B");
  double FromBus; double ToBus; double MaxAng; double MinAng; double G; double B;
  while(in.read_row(FromBus, ToBus, MaxAng, MinAng, G, B)){
    printf("From Bus: %d, To Bus: %d, Max Angle Difference: %f, Min Angle Difference: %f, Conductance: %f, Susecptance: %f\n", int(FromBus), int(ToBus), MaxAng, MinAng, G, B);
  }
}



int main(){

	load_data();
	
	/*********************************************** Declare X Vector ***********************************************/
	//contents {Real Power of Gen n, React Power of Gen n, Volt of Bus m, Volt Angle of Bus m}
	int sizeX = RealPower.size() + ReactPower.size() + Volt.size() + VoltAng.size();
	Real X[sizeX];
	int counter = 0;
	for (double i: RealPower)
		X[counter++] = i;
	for (double i: ReactPower)
		X[counter++] = i;
	for (double i: Volt)
		X[counter++] = i;
	for (double i: VoltAng)
		X[counter++] = i;
		
	
	/*********************************************** Declare/Construct Problem Variables ***********************************************/
	int Nb = FromBus.size();
	int sizeY = 2*RealPowerMax.size() + 2*ReactPowerMax.size() + 2*VoltMax.size() + 2*FromBus.size();
	Real Cost[1], GX[2*Nb], HX[sizeY], f_val[1], G_val[2*Nb];
	ublas::compressed_matrix<double> Z(sizeY, 1);
	// ublas::compressed_matrix<double> L(1, 1);
	
	ublas::compressed_matrix<double> lambda(2*Nb, 1), mu(sizeY, 1), gamma(sizeY, 1);
	
	//printf("Number of Equations: %d, Number of Variables: %d\n", sizeY, sizeX);
	
	for(size_t i = 0; i < mu.size1(); ++i) {
		for(size_t j = 0; j < mu.size2(); ++j) {
			mu(i,j) = 0;
		}
    	}
	
	
	/************************************
	 LOOP WILL GO BACK TO HERE
	*************************************/
	
	
	/*********************************************** Calculate Problem ***********************************************/
	f(X, a, b, c, Cost);
	ublas::compressed_matrix<double> CostuBLAS(1, 1);
	CostuBLAS = RealPointerToDoubleUBlasVec(Cost, 1);
	
	G_func(X, BusID_Gen, BusID_REG, FromBus, ToBus, RealPowerDemand, ReactPowerDemand, Gvec, Bvec, Real_P_REG, React_P_REG, GX);
	ublas::compressed_matrix<double> GXuBLAS(2*Nb, 1);
	GXuBLAS = RealPointerToDoubleUBlasVec(GX, 2*Nb);
	
	H(X, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin, HX);
	ublas::compressed_matrix<double> HXuBLAS(sizeY, 1);
	HXuBLAS = RealPointerToDoubleUBlasVec(HX, sizeY);
	
	/*********************************************** Declare Derivatives ***********************************************/
	vector<vector<Real>> dFdX(1, vector<Real>(sizeX));
	vector<vector<Real>> dHXdX(sizeY, vector<Real>(sizeX));
	vector<vector<Real>> dGXdX(2*Nb, vector<Real>(sizeX));


	/*********************************************** Calculate Derivatives ***********************************************/
	dFdX  = fX(X, sizeX, f_val, 1, a, b, c, dFdX);
	ublas::compressed_matrix<double> dFdXMatrix(dHXdX.size(), dHXdX[0].size());
	dFdXMatrix = RealVectorToDoubleUBlas(dFdX, dFdX.size(), dFdX[0].size());
	
	dHXdX = forwardModeFirstDerivativeH(X, sizeX, HX, sizeY, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, 
								VoltAngMax, VoltAngMin, dHXdX);
	ublas::compressed_matrix<double> dHdxMatrix(dHXdX.size(), dHXdX[0].size());
	dHdxMatrix = RealVectorToDoubleUBlas(dHXdX, dHXdX.size(), dHXdX[0].size());

	dGXdX  = GX_func(X, sizeX, G_val, 2*Nb, BusID_Gen, BusID_REG, FromBus, ToBus, RealPowerDemand, ReactPowerDemand, Gvec, Bvec, Real_P_REG,
					React_P_REG, dGXdX);
	ublas::compressed_matrix<double> dGXdXMatrix(dGXdX.size(), dGXdX[0].size());
	dGXdXMatrix = RealVectorToDoubleUBlas(dGXdX, dGXdX.size(), dGXdX[0].size());
	
	/*for(size_t i = 0; i < dGXdXMatrix.size1(); ++i) {
		for(size_t j = 0; j < dGXdXMatrix.size2(); ++j) {
			cout << fixed;
			cout << std::setprecision (1) <<dGXdXMatrix(i,j) << ' '; // prints rounded numbers
		}
		cout << '\n';
    	}
    	cout << scientific;*/
	
	
	/*********************************************** Declare Lagrangian Functions ***********************************************/
	Real L[0];
	vector<vector<Real>> dLdX(1, vector<Real>(sizeX));
	
	
	/*********************************************** Calculate Lagrangian Functions ***********************************************/
	// L = CostuBLAS + ublas::prod(trans(lambda), GXuBLAS) + ublas::prod(trans(mu), HXuBLAS + Z) - ublas::prod(trans(gamma), uBLASNaturalLog(Z));
	Lag(X, lambda, mu, gamma, Z, L);
	
	cout<<L[0].value() << endl;
	
	dLdX = forwardModeFirstDerivativeL(X, sizeX, lambda, mu, gamma, Z, dLdX);
	ublas::compressed_matrix<double> dLdXMatrix(dLdX.size(), dLdX[0].size());
	dLdXMatrix = RealVectorToDoubleUBlas(dLdX, dLdX.size(), dLdX[0].size());
	
	for(size_t i = 0; i < dLdXMatrix.size1(); ++i) {
		for(size_t j = 0; j < dLdXMatrix.size2(); ++j) {
			cout << fixed;
			cout << std::setprecision (4) <<dLdXMatrix(i,j) << ' '; // prints rounded numbers
		}
		cout << '\n';
    	}
    	cout << scientific;
    
    

	/************************************
	 Testing for f_X
	*************************************/
	// vector<vector<Real>> f_X(1, vector<Real>(sizeX));
	// cout<<"Cost ="<<Cost[0]<<"\n"<<"a.size()"<<a.size()<<"\n";
	// Real f_val[1];
	// f_X  = fX(X, sizeX, f_val, 1, a, b, c, f_X);
	// for (vector<Real> i: f_X)
	// {
	// 	for (Real j: i)
	// 		cout<<j << ' ';
	// 	cout<<"\n";
	// }

	/************************************
	 Testing for G, G_X
	*************************************/
	/*
	int Nb = FromBus.size();
	Real G_val[2*Nb];
	G_func( X, BusID_Gen, BusID_REG, FromBus, ToBus, RealPowerDemand, ReactPowerDemand, Gvec, Bvec, Real_P_REG, React_P_REG, G_val);
	vector<vector<Real>> G_X(2*Nb, vector<Real>(sizeX));
	for(Real val: G_val)
	{
		cout<<val<<" ";
	}
	cout<<endl;
	G_X = GX_func(X, sizeX, G_val, 2*Nb, BusID_Gen, BusID_REG, FromBus, ToBus, RealPowerDemand, ReactPowerDemand, Gvec, 
			 Bvec, Real_P_REG, React_P_REG,G_X );
	cout<<"First derivative"<<endl;
	for(vector<Real> v: G_X)
	{
		for(Real val: v)
			printf("%.1f ",val.value());
		cout<<endl;
	}*/
	
	
	
}








