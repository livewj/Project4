#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <math.h>
#include <iomanip>
#include <fstream>

using namespace arma;
using namespace std;
ofstream ofile;


int magnet(mat& S, int L);
double energy(mat& S, double E, int L);
void ising_create(mat& S, int L);


int main(){
	int seed=1;
	srand(seed);
	int L = 20;    //Number of spins in each dimension (2 dimensions)
	mat S;
	mat S_old;
	double Energy_i;		//Momentary energy
	double deltaE;			//Energy difference due to state transition
	int M_i;			    //Momentary magnetiztion
	double M_mean;
	double maxMCC = 1E5;        
	int accept = 0; //count # of accepted configurations

	S = mat(L+1,L+1,fill::ones);	

	//write numeric energy and magnetization to file for different values of maxMCC:
	string outfile; //write to this file
	outfile = "T_random.txt";
	ofile.open(outfile);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "T:      Accept:" << endl;

	for (double T = 0.; T <= 4.; T+=0.1) {
		accept = 0;

		double M = 0; 				//Mean abs magnetization
		double M_sq = 0; 			//Mean squared magnetization
		double E = 0; 		    	//Total energy
		double E_mean = 0;          //Mean energy
		double EnergySum = 0;		//Needed for calculating mean energy (expectation value)	

		ising_create(S,L); // initial, random matrix. comment if u want up matrix
		S_old=S;
		for (int MCC = 0; MCC < maxMCC; MCC++) {		
			ising_create(S,L);
			deltaE =energy(S,E,L) - energy(S_old,E,L);
			if (deltaE > 0){
				if (rand()*1./RAND_MAX < exp(-deltaE/T)){
					accept += 1;                         //accept a higher energystate
					Energy_i =energy(S,E,L);			 //Momentary energy
					EnergySum += Energy_i; 				 // N*<E>
					M_i = magnet(S,L);				 	 //Momentary magnetization
					M+= M_i;							 //<|M|> mean magnetization

					S_old = S;
					}
				else{
					Energy_i =energy(S_old,E,L);		 //Momentary energy
					EnergySum += Energy_i; 				 // N*<E>

					M_i = magnet(S_old,L);				 
					M+= M_i;							 //<|M|> mean magnetization
					}
				}
			else{
				accept += 1;
				Energy_i =energy(S,E,L);
				EnergySum += Energy_i; 				 // N*<E> energy
				M_i = magnet(S,L);				 
				M+= M_i;							 //<|M|> mean magnetization

				S_old = S;
				}
			}

		E_mean = EnergySum / (1.+maxMCC);
		M_mean = M/(1.+maxMCC);

		ofile << setw(6) << setprecision(8) << T << "     " << accept << endl;

	}

	ofile.close();

	return 0;

}  //end of main



double energy(mat& S, double E, int L){
	E=0;
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			E += S(i,j)*S(i,j+1)+S(i,j)*S(i+1,j);
		}
	}
	return -E;
}

int magnet(mat& S,int L){
	int M_i = 0; //returns M_i to zero, or else it accumulates ALL spin values
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			M_i += S(i,j);
		}
	}
	if(M_i < 0){
		M_i = -M_i;
		return M_i;
	}
	else{
		return M_i;
		}
}


void ising_create(mat& S, int L){
	for(int i = 0; i < L; i++){	
		for(int j = 0; j < L; j++){	
			int I = rand()/RAND_MAX*(L-1);
			int J = rand()/RAND_MAX*(L-1);
			S(I,J) *= -1;
		}
	S(i,L) = S(i,0);
	S(L,i) = S(0,i);
	}
}
