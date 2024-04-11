#include "csv.h"
#include "define2.h"
#include <iostream> 
#include <vector>
#include <codi.hpp>

using Real = codi::RealForward;

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
	char filename1[] = "./data/Case 14 Bus.csv";
	char filename2[] = "./data/Case 14 Generator.csv";
	char filename3[] = "./data/Case 14 Generator Cost.csv";
	char filename4[] = "./data/Case 14 Branch.csv";
	char filename5[] = "./data/Case 14 REG one scen.csv";
  
	// declare temp vectors for Bus Data
	std::vector<double> VoltAng;
	std::vector<double> Volt;
	std::vector<double> VoltMax;
	std::vector<double> VoltMin;
	vector<double> RealPowerDemand, ReactPowerDemand;
	// following code is to read Bus data
	{
		io::CSVReader<7> in(filename1);
		in.read_header(io::ignore_extra_column, "Bus ID", "Real Power Demand", "React Power Demand", "Volt Magnitude", "Max Volt Magnitude", "Min Volt Magnitude", "Volt Angle");
		double Bus; double Pd; double Qd; double V; double Vmax; double Vmin; double Vang;
		
		while(in.read_row(Bus, Pd, Qd, V, Vmax, Vmin, Vang)){
			// get Voltages and Voltage Angles
			VoltAng.insert(VoltAng.end(), Vang);
			Volt.insert(Volt.end(), V);
			VoltMax.insert(VoltMax.end(), Vmax);
			VoltMin.insert(VoltMin.end(), Vmin);
			RealPowerDemand.push_back(Pd);
			ReactPowerDemand.push_back(Qd);
		}
	}
  
	// declare temp vectors for Generator Data
	std::vector<double> RealPower;
	std::vector<double> ReactPower;
	std::vector<double> RealPowerMax;
	std::vector<double> RealPowerMin;
	std::vector<double> ReactPowerMax;
	std::vector<double> ReactPowerMin;
	vector<double> BusID_Gen;
  	// following code is to read Generator data
	{
		io::CSVReader<7> in(filename2);
		in.read_header(io::ignore_extra_column, "Bus ID", "Real Power Output", "React Power Output", "Max Real Power Output", "Min Real Power Output", "Max React Power Output", "Min React Power Output");
		double Bus; double Po; double Qo; double Pmax; double Pmin; double Qmax; double Qmin;
		while(in.read_row(Bus, Po, Qo, Pmax, Pmin, Qmax, Qmin)){
			// get Real and React power output
			RealPower.insert(RealPower.end(), Po);
			ReactPower.insert(ReactPower.end(), Qo);
			RealPowerMax.insert(RealPowerMax.end(), Pmax);
			RealPowerMin.insert(RealPowerMin.end(), Pmin);
			ReactPowerMax.insert(ReactPowerMax.end(), Qmax);
			ReactPowerMin.insert(ReactPowerMin.end(), Qmin);
			BusID_Gen.push_back(Bus);
		}
	}
	
	// declare temp vectors for Generator Cost Data
	std::vector<double> a;
	std::vector<double> b;
	std::vector<double> c;
  	// following code is to read Generator Cost data
	{
		io::CSVReader<3> in(filename3);
		in.read_header(io::ignore_extra_column, "Coef1", "Coef2", "Coef3");
		double coef1; double coef2; double coef3;
		while(in.read_row(coef1, coef2, coef3)){
			// get a, b, and c coefficients for cost function
			a.insert(a.end(), coef1);
			b.insert(b.end(), coef2);
			c.insert(c.end(), coef3);
		}
	}
	
	// declare temp vectors for Branch Data
	vector<double> VoltAngMax, VoltAngMin, Gvec, Bvec;
	vector<double> FromBus, ToBus;
  	// following code is to read Branch Data
	{
		io::CSVReader<6> in(filename4);
		in.read_header(io::ignore_extra_column, "From Bus", "To Bus", "Max Angle Difference", "Min Angle Difference", "G" , "B");
		int FBus, TBus; 
		double MaxAng, MinAng, Gij, Bij;
		while(in.read_row(FBus, TBus, MaxAng, MinAng, Gij, Bij)){
			// get From Bus, To Bus, and Voltage Angle Max and Min values
			FromBus.insert(FromBus.end(), FBus);
			ToBus.insert(ToBus.end(), TBus);
			VoltAngMax.insert(VoltAngMax.end(), MaxAng);
			VoltAngMin.insert(VoltAngMin.end(), MinAng);
			Gvec.push_back(Gij);
			Bvec.push_back(Bij);
		}
	}

	vector<double> Real_P_REG, React_P_REG;
	vector<double> BusID_REG;
	{
		io::CSVReader<7> in(filename5);
		in.read_header(io::ignore_extra_column, "Bus ID", "Real Power Output", "React Power Output", "Max Real Power Output", "Min Real Power Output", "Max React Power Output", "Min React Power Output");
		
		double REG_ID, n,o,m,l,	P, Q; 
		while(in.read_row(REG_ID, n,o,m,l,	P, Q)){
			BusID_REG.push_back(REG_ID); 
			Real_P_REG.push_back(P); 
			React_P_REG.push_back(Q); 

		}	
	}
    
    /*
	// declare X vector with contents {Real Power of Gen n, React Power of Gen n, Volt of Bus m, Volt Angle of Bus m}
	std::vector<double> X;
	X.insert( X.end(), RealPower.begin(), RealPower.end() );
	X.insert( X.end(), ReactPower.begin(), ReactPower.end() );
	X.insert( X.end(), Volt.begin(), Volt.end() );
	X.insert( X.end(), VoltAng.begin(), VoltAng.end() );
	*/
	
	// declare X vector with contents {Real Power of Gen n, React Power of Gen n, Volt of Bus m, Volt Angle of Bus m}
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

//	for (Real i: X)
//		std::cout << i << ' ';
//	printf("\n");
	
	// Declare Problem
	//double fX;
	std::vector<double> GX;
	//std::vector<double> HX;
	//std::vector<double> Z;
	
	Real Cost[1];
	printf("RealPower.Size() = %ld, Volt() = %ld, sizeX = %d\n", RealPower.size(), Volt.size(), sizeX);
	int sizeY = 2*RealPowerMax.size() + 2*ReactPowerMax.size() + 2*VoltMax.size() + 2*FromBus.size();
	Real HX[sizeY], Z[sizeY];
	vector<vector<Real>> dHXdX(sizeY, vector<Real>(sizeX));
	
	printf("Number of Equations: %d, Number of Variables: %d\n", sizeY, sizeX);
	
	// Construct Problem
	// f(X, a, b, c, Cost);
//	for (Real i: fX)
//		std::cout << i << ' ';
//	printf("\n");

	// H(X, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin, HX);
//	for (Real i: HX)
//		std::cout << i << std::endl;
//	printf("\n");

	// for (int i=0; i < sizeY; i++) {
	// 	Z[i] = HX[i]*-1.0;
	// }
//	for (Real i: Z)
//		std::cout << i << ' ';
//	printf("\n");

	// Calculate Derivatives
	/*dHXdX = forwardModeFirstDerivativeH(X, sizeX, HX, sizeY, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, 
								VoltAngMax, VoltAngMin, dHXdX);
	for (vector<Real> i: dHXdX) {
		for (Real j: i) {
			std::cout << j << ' ';
		}
		std::cout << std::endl;
	}*/

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
	}
	
	
}








