#include "csv.h"
#include "define4.h"
#include <iostream> 
#include <vector>
#include "read_data.h"
#include "second_derivs.h"
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
	sizeX = RealPower.size() + ReactPower.size() + Volt.size() + VoltAng.size();
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
	
	H(X, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin, HX);
	ublas::compressed_matrix<double> HXuBLAS(sizeY, 1);
	HXuBLAS = RealPointerToDoubleUBlasVec(HX, sizeY);
	
	for(size_t i = 0; i < Z.size1(); ++i) {
		for(size_t j = 0; j < Z.size2(); ++j) {
			Z(i,j) = -1*HXuBLAS(i,j);
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
	//ublas::compressed_matrix<double> HXuBLAS(sizeY, 1);
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
	
	/*for(size_t i = 0; i < dFdXMatrix.size1(); ++i) {
		for(size_t j = 0; j < dFdXMatrix.size2(); ++j) {
			cout << fixed;
			cout << std::setprecision (4) <<dFdXMatrix(i,j) << ' '; // prints rounded numbers
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
	
	// cout << L[0].value() << endl;
	
	dLdX = forwardModeFirstDerivativeL(X, sizeX, lambda, mu, gamma, Z, dLdX);
	ublas::compressed_matrix<double> dLdXMatrix(dLdX.size(), dLdX[0].size());
	dLdXMatrix = RealVectorToDoubleUBlas(dLdX, dLdX.size(), dLdX[0].size());
	
	/*for(size_t i = 0; i < dLdXMatrix.size1(); ++i) {
		for(size_t j = 0; j < dLdXMatrix.size2(); ++j) {
			cout << fixed;
			cout << std::setprecision (4) <<dLdXMatrix(i,j) << ' '; // prints rounded numbers
		}
		cout << '\n';
    }
    cout << scientific;*/
	vector<double> SX(sizeX);
	//double SX[sizeX];
	for(int i = 0;i<sizeX;++i)
	{
		SX[i] = X[i].value();
	}
	second_deriv(SX , lambda, mu, gamma, Z);
    	
    	
	/*********************************************** Compute Matrix for cuSolver ***********************************************/
	ublas::compressed_matrix<double> muMatrix(sizeY, sizeY), ZMatrixInverse(sizeY, sizeY), M(sizeX, sizeX), N(sizeX,1);
	muMatrix = uBLASVectorToMatrix(mu);
	ZMatrixInverse = uBLASVectorToMatrix(Z);
	ZMatrixInverse = InvertDiagonalMatrix(ZMatrixInverse);
	
	// Computing N with temp variables
	ublas::compressed_matrix<double> temp1(sizeX, sizeY), temp2(sizeY,1);
	temp1 = ublas::prod(trans(dHdxMatrix), ZMatrixInverse);
	temp2 = gamma + ublas::prod(muMatrix, HXuBLAS);
	N = trans(dLdXMatrix) + ublas::prod(temp1, temp2);
	
	// Computing M with temp variables
	ublas::compressed_matrix<double> temp3(sizeY, sizeY);
	temp3 = ublas::prod(temp1, muMatrix);
	M = ublas::prod(temp3, dHdxMatrix);
	
	
	/*********************************************** Assign Values for cuSolver ***********************************************/
	//Create cuSolver matrix
	
	
	//Create cuSolver vector
	
	
	//printf("Number of Variables: %ld\n", ublas::prod(temp1, temp2).size2());
	/*for(size_t i = 0; i < M.size1(); ++i) {
		for(size_t j = 0; j < M.size2(); ++j) {
			cout << M(i,j) << ' ';
		}
		cout << '\n';
    }*/
    
    //NOTE*** size2 = rows size1 = columns
    
	/************************************
	CUSOLVER WILL GO HERE
	*************************************/
	
	/*********************************************** Compute Update Values ***********************************************/
	// We assume deltaX and deltalambda are given from the cuSolver output
	double chi = 0.99995, sigma = 0.1, alphaP, alphaD;
	ublas::compressed_matrix<double> deltalambda(2*Nb, 1), deltaX(sizeX, 1); // will be solved for on GPU
    ublas::compressed_matrix<double> deltaZ(sizeY, 1), deltamu(sizeY, 1), newGamma(1,1);
    
    /*for(size_t i = 0; i < deltaX.size1(); ++i) {
		for(size_t j = 0; j < deltaX.size2(); ++j) {
			deltaX(i,j) = 0;
		}
    }
    for(size_t i = 0; i < deltalambda.size1(); ++i) {
		for(size_t j = 0; j < deltalambda.size2(); ++j) {
			deltalambda(i,j) = 0;
		}
    }*/
    
    deltaZ = -1*HXuBLAS - Z - ublas::prod(dHdxMatrix, deltaX);
	ublas::compressed_matrix<double> temp4(sizeY,1);
	temp4 = gamma - ublas::prod(muMatrix, deltaZ);
    deltamu = -1*mu + ublas::prod(ZMatrixInverse, temp4);
    
    /*for(size_t i = 0; i < deltaZ.size1(); ++i) {
		for(size_t j = 0; j < deltaZ.size2(); ++j) {
			cout << deltaZ(i,j) << ' ';
		}
		cout << '\n';
    }
    cout << '\n';*/

	if(chi*-1*MaxFractionOverXiLTOne(Z, deltaZ) < 1) {
		alphaP = chi*-1*MaxFractionOverXiLTOne(Z, deltaZ);
	}
	else {
		alphaP = 1;
	}
	cout << alphaP << '\n';
	cout << '\n';
	if(chi*-1*MaxFractionOverXiLTOne(mu, deltamu) < 1) {
		alphaD = chi*-1*MaxFractionOverXiLTOne(mu, deltamu);
	}
	else {
		alphaD = 1;
	}
	cout << alphaD << '\n';
	cout << '\n';
	
	
	/*********************************************** Update Problem Variables ***********************************************/
	Real deltaXReal[sizeX];
	DoubleUBlasVecToRealPointer(alphaP*deltaX, deltaXReal); // changes vector to point
	RealPointerAdd(X, deltaXReal, sizeX); // adds deltaXReal to X
	Z = Z + alphaP*deltaZ;
	lambda = lambda + alphaD*deltalambda;
	mu = mu + alphaD*deltamu;
	
	newGamma = sigma*ublas::prod(trans(Z), mu);
	for(size_t i = 0; i < gamma.size1(); ++i) {
		for(size_t j = 0; j < gamma.size2(); ++j) {
			gamma(i,j) = newGamma(0,0);
		}
    }
    
    /*for(size_t i = 0; i < Z.size1(); ++i) {
		for(size_t j = 0; j < Z.size2(); ++j) {
			cout << Z(i,j) << ' ';
		}
		cout << '\n';
    }*/
}








