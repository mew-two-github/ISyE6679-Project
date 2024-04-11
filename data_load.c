#include<csv.h>

#include<define.h>

#include<iostream> 

#include<vector>



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

	char filename1[] = "Case 14 Bus.csv";

	char filename2[] = "Case 14 Generator.csv";

	char filename3[] = "Case 14 Generator Cost.csv";

	char filename4[] = "Case 14 Branch.csv";

	//char filename5[] = "Case 14 REG.csv";

  

	// declare temp vectors for Bus Data

	std::vector<double> VoltAng;

	std::vector<double> Volt;

	std::vector<double> VoltMax;

	std::vector<double> VoltMin;

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

		}

	}

  

	// declare temp vectors for Generator Data

	std::vector<double> RealPower;

	std::vector<double> ReactPower;

	std::vector<double> RealPowerMax;

	std::vector<double> RealPowerMin;

	std::vector<double> ReactPowerMax;

	std::vector<double> ReactPowerMin;

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

	std::vector<double> FromBus;

	std::vector<double> ToBus;

	std::vector<double> VoltAngMax;

	std::vector<double> VoltAngMin;

  	// following code is to read Generator Cost data

	{

		io::CSVReader<6> in(filename4);

		in.read_header(io::ignore_extra_column, "From Bus", "To Bus", "Max Angle Difference", "Min Angle Difference", "G" , "B");

		double FBus; double TBus; double MaxAng; double MinAng; double G; double B;

		while(in.read_row(FBus, TBus, MaxAng, MinAng, G, B)){

			// get From Bus, To Bus, and Voltage Angle Max and Min values

			FromBus.insert(FromBus.end(), FBus);

			ToBus.insert(ToBus.end(), TBus);

			VoltAngMax.insert(VoltAngMax.end(), MaxAng);

			VoltAngMin.insert(VoltAngMin.end(), MinAng);

		}

	}

    

	// declare X vector with contents {Real Power of Gen n, React Power of Gen n, Volt of Bus m, Volt Angle of Bus m}

	std::vector<double> X;

	X.insert( X.end(), RealPower.begin(), RealPower.end() );

	X.insert( X.end(), ReactPower.begin(), ReactPower.end() );

	X.insert( X.end(), Volt.begin(), Volt.end() );

	X.insert( X.end(), VoltAng.begin(), VoltAng.end() );



	for (double i: X)

		std::cout << i << ' ';

	printf("\n");

	

	// Declare Problem

	double fX;

	std::vector<double> GX;

	std::vector<double> HX;

	std::vector<double> Z;

	// Construct Problem

	fX = f(X, a, b, c);

	printf("Cost Sum: %f\n", fX);

	HX = H(X, RealPowerMax, RealPowerMin, ReactPowerMax, ReactPowerMin, VoltMax, VoltMin, FromBus, ToBus, VoltAngMax, VoltAngMin);

	for (double i: HX)

		std::cout << i << ' ';

	printf("\n");

	for (int i=0; i < HX.size(); i++) {

		Z.insert(Z.end(), HX[i]*-1.0);

	}

	for (double i: Z)

		std::cout << i << ' ';

	printf("\n");

}















