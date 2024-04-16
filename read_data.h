#ifndef READ_DATA_H_
#define READ_DATA_H_

#include <iostream> 
#include <vector>
#include "csv.h"
#define ucmd ublas::compressed_matrix<double>

using namespace std;

/***************************************************


Global Variables


*****************************************************/
// Bus data
extern vector<double> Volt, VoltMax, VoltMin, VoltAng, RealPowerDemand, ReactPowerDemand;
// Generator Data
extern vector<double> RealPower, ReactPower, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, BusID_Gen;
// Generator Cost Data
extern vector<double>  a, b, c;
// Branch Data
extern vector<double> VoltAngMax, VoltAngMin, Gvec, Bvec, FromBus, ToBus;
// declare temp vectors for REG Data
extern vector<double> Real_P_REG, React_P_REG,	BusID_REG;
extern int Nb, sizeY, sizeX;

namespace ublas = boost::numeric::ublas;



void load_data()
{
/*********************************************** Declare Filenames ***********************************************/
	char filename1[] = "./data/Case 14 Bus.csv";
	char filename2[] = "./data/Case 14 Generator.csv";
	char filename3[] = "./data/Case 14 Generator Cost.csv";
	char filename4[] = "./data/Case 14 Branch.csv";
	char filename5[] = "./data/Case 14 REG.csv";
  
  
  	/*********************************************** Read Bus Data ***********************************************/
	// declare temp vectors for Bus Data
	
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
	
  
  	/*********************************************** Read Generator Data ***********************************************/
	// declare temp vectors for Generator Data

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
	
	
	/*********************************************** Read Generator Cost Data ***********************************************/
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
	
	
	/*********************************************** Read Branch Data ***********************************************/
	// declare temp vectors for Branch Data

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
	
	
	/*********************************************** Read REG Data ***********************************************/

	// following code is to read REG Data
	{
		io::CSVReader<7> in(filename5);
		in.read_header(io::ignore_extra_column, "Bus ID", "Real Power Output", "React Power Output", "Max Real Power Output", "Min Real Power Output", "Max React Power Output", "Min React Power Output");
		
		double REG_ID, n,o,m,l,	P, Q; 
		while(in.read_row(REG_ID, n,o,m,l,	P, Q)){
			// get Bus ID, Real Power, and React Power values
			BusID_REG.push_back(REG_ID); 
			Real_P_REG.push_back(P); 
			React_P_REG.push_back(Q);
		}	
	}
    Nb = FromBus.size();
    sizeY = 2*RealPowerMax.size() + 2*ReactPowerMax.size() + 2*VoltMax.size() + 2*FromBus.size();
}



#endif