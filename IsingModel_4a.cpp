/* 
   Program to solve the two-dimensional Ising model 
   with zero external field and no parallelization
   The coupling constant J is set to J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis aolgorithm  is used as well as periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
   Run as
   ./executable Outputfile numberof spins number of MC cycles initial temp final temp tempstep
   ./test.x Lattice 100 10000000 2.1 2.4 0.01
   Compile and link as 
   c++ -O3 -std=c++11 -Rpass=loop-vectorize -o Ising.x IsingModel.cpp -larmadillo
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <random>
#include <armadillo>
#include <string>
using namespace  std;
using namespace arma;
// output file
ofstream ofile;

// inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void InitializeLattice(int, mat &, double&, double&);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, vec &, int&);

// Main program begins here

int main(int argc, char* argv[])
{
  string filename;
  int NSpins, MCcycles;
  double InitialTemp, FinalTemp, TempStep;
  if (argc <= 5) {
    cout << "Bad Usage: " << argv[0] << 
      " read output file, Number of spins, MC cycles, initial and final temperature and tempurate step" << endl;
    exit(1);
  }
  if (argc > 1) {
    filename=argv[1];
    NSpins = atof(argv[2]);
    MCcycles = atof(argv[3]);    
    InitialTemp = atof(argv[4]);
    FinalTemp = atof(argv[5]);
    TempStep = atof(argv[6]);
  }
  // Declare new file name and add lattice size to file name
  string fileout = filename;
  string argument = to_string(NSpins);
  fileout.append(argument);

  // Start Monte Carlo sampling by looping over the selcted Temperatures 
  // InitialTemp = 1.0, Finaltemp = 1.0, TempStep = 0.1
  for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
    vec ExpectationValues = zeros<mat>(5);
    int accept = 0;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(NSpins, MCcycles, Temperature, ExpectationValues, accept);
  }

  double T = 1.0;

  //Expectation value for L=2, found analytically
  double Z_exact = 2*exp(8.) + 2*exp(-8.) + 12; //partition function
  double E_mean_exact = (2*(-8)*exp(8.) + 2*8*exp(-8))/Z_exact; 
  double E_mean_sq_exact = (2*64*exp(8) + 2*64*exp(-8))/Z_exact;
  double M_mean_exact = 2.*(4*exp(8) + 8.)/Z_exact;
  double M_mean_sq_exact = (32*exp(8) + 32)/Z_exact;
  double sigmaE_sq = (E_mean_sq_exact - E_mean_exact*E_mean_exact);  //Standard deviations
  double sigmaM_sq = (M_mean_sq_exact - M_mean_exact*M_mean_exact);
  double Cv_exact = sigmaE_sq/T; //T = 1.0 kT/J
  double chi_exact = sigmaM_sq/T;

  cout << MCcycles << endl;
  cout<<"Exact <E> = " <<setprecision(15) <<E_mean_exact << endl;
  cout<<"Exact <|M|> = " <<setprecision(15) <<M_mean_exact << endl;
  cout << "Exact Cv = " << setprecision(15) << Cv_exact << endl;
  cout << "Exact chi = " << setprecision(15) << chi_exact << endl; 
  
  return 0;
}






// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void MetropolisSampling(int NSpins, int MCcycles, double Temperature, vec &ExpectationValues,int& accept)
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  // Initialize the lattice spin values
  mat SpinMatrix = zeros<mat>(NSpins,NSpins);
  //    initialize energy and magnetization
  double Energy = 0.;     
  double MagneticMoment = 0.;
  // initialize array for expectation values
  InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment);
  // setup array for possible energy changes (Metropolis algo)
  vec EnergyDifference = zeros<mat>(17); 
  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  // Start Monte Carlo cycles
  for (int cycles = 1; cycles <= MCcycles; cycles++){
    // The sweep over the lattice, looping over all spin sites
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
	int ix = (int) (RandomNumberGenerator(gen)*(double)NSpins);
	int iy = (int) (RandomNumberGenerator(gen)*(double)NSpins);
	int deltaE =  2*SpinMatrix(ix,iy)*
	  (SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
	   SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
	   SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
	   SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
	if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
    accept += 1;
	  SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
	  MagneticMoment += (double) 2*SpinMatrix(ix,iy);
	  Energy += (double) deltaE;
  }
}}
    // update expectation values  for local node
    ExpectationValues(0) += Energy;    
    ExpectationValues(1) += Energy*Energy;
    ExpectationValues(2) += MagneticMoment;    
    ExpectationValues(3) += MagneticMoment*MagneticMoment; 
    ExpectationValues(4) += fabs(MagneticMoment);
  }
//Numeric values for L=2
double norm = 1.0/((double) (MCcycles));  // divided by  number of cycles 
double E_ExpectationValues = ExpectationValues(0)*norm;
double E2_ExpectationValues = ExpectationValues(1)*norm;
double M_ExpectationValues = ExpectationValues(2)*norm;
double M2_ExpectationValues = ExpectationValues(3)*norm;
double Mabs_ExpectationValues = ExpectationValues(4)*norm;
// all expectation values are per spin, divide by 1/NSpins/NSpins
double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues);///NSpins/NSpins;
double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues);///NSpins/NSpins;

double T = 1.0;

cout <<"Numeric <E> = " << setprecision(15) << E_ExpectationValues << endl;///NSpins/NSpins << endl;
cout<<"Numeric <|M|> = " <<setprecision(15) << Mabs_ExpectationValues << endl; ///NSpins/NSpins << endl;
cout <<"Numeric Cv = " << setprecision(15) << Evariance/T/T << endl;
cout << "Numeric chi = " << setprecision(15) << Mvariance/T << endl;

} // end of Metropolis sampling over spins




// function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment)
{
  /*
  // setup spin matrix and initial magnetization: ALL SPINS UP
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      SpinMatrix(x,y) = 1.0; // spin orientation for the ground state, all spins up
      MagneticMoment +=  (double) SpinMatrix(x,y);
    }
  }
*/

  //Comment/uncomment this section 
  //setup spin matrix and initial magnetization: RANDOM ORIENTATION
  for (int x =0; x < NSpins; x++) {
     for (int y= 0; y < NSpins; y++){
       double a =  rand()/double(RAND_MAX);
       if (a < 0.5) {
       SpinMatrix(x,y) = -1; // spin orientation for the random state
       MagneticMoment += SpinMatrix(x,y);}
       else{
       SpinMatrix(x,y) = 1; 
       }
     }
  }

  // setup initial energy
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      Energy -=  (double) SpinMatrix(x,y)*
	(SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
	 SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
    }
  }
}// end function initialise



    
