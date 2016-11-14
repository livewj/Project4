/* 
   TO DO: look up parallel RNG, check ranlux48
   TO DO: performance analysis and vectorization, find bottlenecks
   TO DO: make smarter array for energy differences
   Program to solve the two-dimensional Ising model 
   with zero external field and no parallelization
   Parallel version using MPI
   The coupling constant J is set to J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis aolgorithm  is used as well as periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
   
   Compile and link as 
   mpicxx -O3 -std=c++11 -Rpass=loop-vectorize -o Ising.x IsingModel.cpp -larmadillo
   
   Run as
  (mpirun -n 4 ./executable Outputfile numberof spins number of MC cycles initial temp final temp tempstep)
   mpirun -n 4 ./Ising.x Lattice 100 10000000 2.1 2.4 0.01

*/

#include "mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
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
void MetropolisSampling(int, long&, int, double, vec &);
// prints to file the results of the calculations  
void WriteResultstoFile(int, int, double, vec);
// ran2 for uniform deviates, initialize with negative seed.
double ran2(long *);

// Main program begins here

int main(int argc, char* argv[])
{
  string filename;
  int NSpins, MonteCarloCycles;
  double InitialTemp, FinalTemp, TempStep;
  int NProcesses, RankProcess;
  long idum;
  //  MPI initializations
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &NProcesses);
  MPI_Comm_rank (MPI_COMM_WORLD, &RankProcess);
  if (RankProcess == 0 && argc <= 5) {
    cout << "Bad Usage: " << argv[0] << 
      " read output file, Number of spins, MC cycles, initial and final temperature and tempurate step" << endl;
    exit(1);
  }
  if ((RankProcess == 0) && (argc > 1)) {
    filename=argv[1];
    NSpins = atoi(argv[2]);
    MonteCarloCycles = atoi(argv[3]);    
    InitialTemp = atof(argv[4]);
    FinalTemp = atof(argv[5]);
    TempStep = atof(argv[6]);
  }
  // Declare new file name and add lattice size to file name, only master node opens file
  if (RankProcess == 0) {
    string fileout = filename;
    string argument = to_string(NSpins);
    fileout.append(argument);
    ofile.open(fileout);
    //Write header to file
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "      T:            <E>/L^2:           Cv:         <M>/L^2:        X:           <|M|>/L^2: " << endl;
  }



  /*
  Determine number of intervall which are used by all processes
  myloop_begin gives the starting point on process my_rank
  myloop_end gives the end point for summation on process my_rank
  */
  int no_intervalls = MonteCarloCycles/NProcesses;
  int myloop_begin = RankProcess*no_intervalls + 1;
  int myloop_end = (RankProcess+1)*no_intervalls;
  if ( (RankProcess == NProcesses-1) &&( myloop_end < MonteCarloCycles) ) myloop_end = MonteCarloCycles;


  // broadcast to all nodes common variables since only master node reads from command line
  MPI_Bcast (&MonteCarloCycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&NSpins, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&InitialTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&FinalTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&TempStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // every node has its own seed for the random numbers, this is important else
  // if one starts with the same seed, one ends with the same random numbers
  idum = -1 - RankProcess; // random starting point

  // Start Monte Carlo sampling by looping over the selected Temperatures
  double  TimeStart, TimeEnd, TotalTime;
  TimeStart = MPI_Wtime();
  for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
    vec LocalExpectationValues = zeros<mat>(5);
    // Start Monte Carlo computation and get local expectation values
    MetropolisSampling(NSpins, idum, MonteCarloCycles, Temperature, LocalExpectationValues);
    // Find total average
    vec TotalExpectationValues = zeros<mat>(5);
    for( int i =0; i < 5; i++){
      MPI_Reduce(&LocalExpectationValues[i], &TotalExpectationValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if ( RankProcess == 0) WriteResultstoFile(NSpins, MonteCarloCycles*NProcesses, Temperature, TotalExpectationValues);
  }
  if(RankProcess == 0)  ofile.close();  // close output file
  TimeEnd = MPI_Wtime();
  TotalTime = TimeEnd-TimeStart;
  if ( RankProcess == 0) {
    cout << "Time = " <<  TotalTime  << " sec, on number of processors: "  << NProcesses  << endl;
  }
  // End MPI
  MPI_Finalize (); 
  return 0;
}


// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void MetropolisSampling(int NSpins, long& idum, int MonteCarloCycles, double Temperature, vec &ExpectationValues)
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
  // setup array for possible energy changes
  vec EnergyDifference = zeros<mat>(17); 
  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  // Start Monte Carlo experiments
  int AllSpins = NSpins*NSpins;
  for (int cycles = 1; cycles <= MonteCarloCycles; cycles++){
    // The sweep over the lattice, looping over all spin sites
    for(int Spins =0; Spins < AllSpins; Spins++) {
      int ix = (int) (ran2(&idum)*(double)NSpins);
      int iy = (int) (ran2(&idum)*(double)NSpins);
      int deltaE =  2*SpinMatrix(ix,iy)*
	(SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
	 SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
	 SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
	 SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
      if ( ran2(&idum) <= EnergyDifference(deltaE+8) ) {
	SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
	MagneticMoment += 2.0*SpinMatrix(ix,iy);
	Energy += (double) deltaE;
      }
    }
    // update expectation values  for local node after a sweep through the lattice
    ExpectationValues(0) += Energy;    
    ExpectationValues(1) += Energy*Energy;
    ExpectationValues(2) += MagneticMoment;    
    ExpectationValues(3) += MagneticMoment*MagneticMoment; 
    ExpectationValues(4) += fabs(MagneticMoment);
  }
} // end of Metropolis sampling over spins

// function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment)
{
  // setup spin matrix and initial magnetization using cold start, all spins pointing up or down
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      SpinMatrix(x,y) = 1.0; // spin orientation for the ground state
      MagneticMoment +=  (double) SpinMatrix(x,y);
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
}// end function initialize



void WriteResultstoFile(int NSpins, int MonteCarloCycles, double temperature, vec ExpectationValues)
{
  double norm = 1.0/((double) (MonteCarloCycles));  // divided by  number of cycles 
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm;
  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double AllSpins = 1.0/((double) NSpins*NSpins);
  double HeatCapacity = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)*AllSpins/temperature/temperature;
  double MagneticSusceptibility = (M2_ExpectationValues - M_ExpectationValues*M_ExpectationValues)*AllSpins/temperature;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << E_ExpectationValues*AllSpins; //<E>
  ofile << setw(15) << setprecision(8) << HeatCapacity;                 //Cv
  ofile << setw(15) << setprecision(8) << M_ExpectationValues*AllSpins; //<M>
  ofile << setw(15) << setprecision(8) << MagneticSusceptibility;       //X
  ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues*AllSpins << endl; //<|M|>
} // end output function



/*
** The function 
**         ran2()
** is a long periode (> 2 x 10^18) random number generator of 
** L'Ecuyer and Bays-Durham shuffle and added safeguards.
** Call with idum a negative integer to initialize; thereafter,
** do not alter idum between sucessive deviates in a
** sequence. RNMX should approximate the largest floating point value
** that is less than 1.
** The function returns a uniform deviate between 0.0 and 1.0
** (exclusive of end-point values).
*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int            j;
  long           k;
  static long    idum2 = 123456789;
  static long    iy=0;
  static long    iv[NTAB];
  double         temp;

  if(*idum <= 0) {
    if(-(*idum) < 1) *idum = 1;
    else             *idum = -(*idum);
    idum2 = (*idum);
    for(j = NTAB + 7; j >= 0; j--) {
      k     = (*idum)/IQ1;
      *idum = IA1*(*idum - k*IQ1) - k*IR1;
      if(*idum < 0) *idum +=  IM1;
      if(j < NTAB)  iv[j]  = *idum;
    }
    iy=iv[0];
  }
  k     = (*idum)/IQ1;
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if(*idum < 0) *idum += IM1;
  k     = idum2/IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j     = iy/NDIV;
  iy    = iv[j] - idum2;
  iv[j] = *idum;
  if(iy < 1) iy += IMM1;
  if((temp = AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran2()


