//   Project 4 d)
//   Program to compute probability distribution (of energy) of steady state in
//   a 20x20 grid periodic Ising model. Histogram is plotted and computed 
//   with periodic boundary conditions. Values are plotted with python script. 

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <random>
#include "lib.h"
#include <armadillo>
using namespace arma;
using namespace  std;
ofstream ofile;

inline int PeriodicBoundary(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}

// Function to read in data from screen  
void read_input(int&, double&, int&);
// Function to initialise energy and magnetization
void initialize(int, mat , double&, double&);
// The metropolis algorithm 
void Metropolis(int, int, double, vec, int&);
// prints to file the results of the calculations  
void output(int, double);

int main(int argc, char* argv[])
{
  char *outfilename;
  long idum;
  int **spin_matrix, n_spins, mcs;
  double w[17], temperature, E;
  int steady_state_tolerance_cycles = 5E3;
  double meanE2, meanE, Evariance, norm;

  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] << 
      " read also output file on same line" << endl;
    exit(1);
  }
  else { outfilename=argv[1]; }
  read_input(n_spins, temperature, mcs);
  // Normalize with the number of spins in the grid:
  double norm2 = 1.0/((double) (n_spins*n_spins));
  // Write header in output file:
  ofile.open(outfilename);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "Final number of Monte Carlo trials: " << mcs << ", and temperature: " << temperature << endl;
  ofile << "Energy for each MC cycle after reached steady state. Computation started after ";
  ofile << steady_state_tolerance_cycles << " cycles." << endl;
  //    Read in initial values such as size of lattice, temp and cycles
  idum = -1; // random starting point
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  //    initialise energy and magnetization 
  E = 0.;
  meanE2 = meanE = 0.0;
  // setup array for possible energy changes
  for( int de =-8; de <= 8; de++) w[de+8] = 0;
  for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
  // initialise array for expectation values
  initialize(n_spins, SpinMatrix, Energy, MagneticMoment);
  // start Monte Carlo computation:
  for (int cycles = 1; cycles <= mcs; cycles++){
    Metropolis(n_spins, idum, spin_matrix, E, w);
    // update expectation values if tolerance cycle is passed:
    if (cycles >= steady_state_tolerance_cycles) {
      output(n_spins, E);
      // We calculate the averages per spin this time:
      meanE += E*norm2; meanE2 += E*E*norm2*norm2;
      }

  }
  // Calculate variance in energy for the cycles that occured after
  // reaching the steady state:
  norm = 1/((double) (mcs - steady_state_tolerance_cycles+1));
  meanE *= norm;
  meanE2 *= norm;
  
  // Print some results for verification:
  cout << "mean E per spin:   " << meanE << endl;
  cout << "mean E^2 per spin: " << meanE2 << endl;

  Evariance = meanE2 - meanE*meanE;
  ofile << " Final variance in energy per spin: ";
  ofile << setw(15) << setprecision(8) << Evariance << endl;
  // print results
  free_matrix((void **) spin_matrix); // free memory
  
  ofile.close();  // close output file
  return 0;
}

// read in input data
void read_input(int& n_spins, double& temperature, int& mcs)
{
  cout << "Final number of MC trials: ";
  cin >> mcs;
  cout << "Lattice size or number of spins (x and y equal): ";
  cin >> n_spins;
  cout << "Temperature with dimension energy: ";
  cin >> temperature;
} // end of function read_input


// function to initialise energy, spin matrix and magnetization
void initialize(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment)
{
  // setup spin matrix and initial magnetization: ALL SPINS UP
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      SpinMatrix(x,y) = 1.0; // spin orientation for the ground state, all spins up
      MagneticMoment +=  (double) SpinMatrix(x,y);
    }
  }
/*
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
*/

  // setup initial energy
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      Energy -=  (double) SpinMatrix(x,y)*
  (SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
   SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
    }
  }
}// end function initialise

void Metropolis(int NSpins, int MCcycles, double Temperature, vec &ExpectationValues,int& accept)
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
  initialize(NSpins, SpinMatrix, Energy, MagneticMoment);
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
      }
    }
    // update expectation values  for local node
    ExpectationValues(0) += Energy;    
    ExpectationValues(1) += Energy*Energy;
    ExpectationValues(2) += MagneticMoment;    
    ExpectationValues(3) += MagneticMoment*MagneticMoment; 
    ExpectationValues(4) += fabs(MagneticMoment);
  }
} // end of Metropolis sampling over spins

void output(int n_spins, double E)
{
  double norm2 = 1.0/(n_spins*n_spins);   // divide by the total number of spins
  double energy = E*norm2;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  ofile << setw(15) << setprecision(8) << energy << endl;
} // end output function