
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

using namespace  std;
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to read in data from screen  
void read_input(int&, double&, int&);
// Function to initialise energy and magnetization
void initialize(int, double, int **, double&);
// The metropolis algorithm 
void Metropolis(int, long&, int **, double&, double *);
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
  initialize(n_spins, temperature, spin_matrix, E);
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
void initialize(int n_spins, double temperature, int **spin_matrix, 
                double& E)
{
  // Setup spin matrix and intial magnetization; all spins up.
  // Commment/uncomment this section depending on use
  // Ground state configuration gives good convergence for low temperatures
  
  // for(int y =0; y < n_spins; y++) {
  //   for (int x= 0; x < n_spins; x++){
  //     spin_matrix[y][x] = 1; // spin orientation for the ground state
  //   }
  // }
  
  // Setup spin matrix and intial magnetization; random start configuration.
  // Commment/uncomment this section depending on use.
  // Random spin configuration gives good convergence for high temperatures

  long idum_dum = 2;
  for(int y =0; y < n_spins; y++) {
     for (int x= 0; x < n_spins; x++){
       double a = (double) ran1(&idum_dum);
       int r = 1;
       if (a < 0.5) {r = -1;}
       spin_matrix[y][x] = r; // spin orientation for the random state
     }
   }

  // setup initial energy:
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
  (spin_matrix[periodic(y,n_spins,-1)][x] +
   spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}// end function initialise

void Metropolis(int n_spins, long& idum, int **spin_matrix, double &E, double *w)
{
  // loop over all spins
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      int ix = (int) (ran1(&idum)*(double)n_spins);
      int iy = (int) (ran1(&idum)*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
  (spin_matrix[iy][periodic(ix,n_spins,-1)]+
   spin_matrix[periodic(iy,n_spins,-1)][ix] +
   spin_matrix[iy][periodic(ix,n_spins,1)] +
   spin_matrix[periodic(iy,n_spins,1)][ix]);
      if ( ran1(&idum) <= w[deltaE+8] ) {
  spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
        E += (double) deltaE;
      }
    }
  }
} // end of Metropolis sampling over spins

void output(int n_spins, double E)
{
  double norm2 = 1.0/(n_spins*n_spins);   // divide by the total number of spins
  double energy = E*norm2;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  ofile << setw(15) << setprecision(8) << energy << endl;
} // end output function