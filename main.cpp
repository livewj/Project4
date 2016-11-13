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

	int M_i;					//Momentary magnetiztion
	double M = 0; 				//Mean abs magnetization
	double M_sq = 0; 			//Mean squared magnetization
	double M_mean;
	double M_mean_sq;
	
	double E = 0; 			//Total energy
	double Energy_i;		//Momentary energu
	double E_mean = 0;      //Mean energy
	double E_mean_sq = 0;   //Mean squared energy
	double Cv;
	double chi;


	double EnergySum = 0;		//Needed for calculating mean energy (expectation value)
	double EnergySum_sq = 0;	//For mean sqrd energy, needed for sigmaE and heat capacity 

	//Expectation value for L=2, found analytically
	double T = 1.;			//Units: k_b*T/J
	double Z_exact = 2*exp(8.) + 2*exp(-8.) + 12; //partition function
	double E_mean_exact = (2*(-8)*exp(8.) + 2*8*exp(-8))/Z_exact;	
	double E_mean_sq_exact = (2*64*exp(8) + 2*64*exp(-8))/Z_exact;
	double M_mean_exact = 2.*(4*exp(8) + 8.)/Z_exact;
	double M_mean_sq_exact = (32*exp(8) + 32)/Z_exact;
	double sigmaE_sq = E_mean_sq_exact - E_mean_exact*E_mean_exact;  //Standard deviations
	double sigmaM_sq = M_mean_sq_exact - M_mean_exact*M_mean_exact;
	double Cv_exact = sigmaE_sq/T; //T = 1.0 kT/J
	double chi_exact = sigmaM_sq/T;


	double deltaE;			//Energy difference due to state transition

	int maxMCC = 10;

	mat S;
	mat S_old;
	S = mat(L+1,L+1,fill::ones);		

	ising_create(S,L); // initial
	S_old=S;
	for (int MCC = 0; MCC < maxMCC; MCC++) {		
		ising_create(S,L);
		deltaE =energy(S,E,L) - energy(S_old,E,L);
		if (deltaE > 0){
			if (rand()*1./RAND_MAX < exp(-deltaE/T)){

				Energy_i =energy(S,E,L);			 //Momentary energy
				EnergySum += Energy_i; 				 // N*<E>
				EnergySum_sq += Energy_i*Energy_i; 	 //N*<E^2> 

				M_i = magnet(S,L);				 	 //Momentary magnetization
				M+= M_i;							 //<|M|> mean magnetization
				M_sq += M_i*M_i;					 //<M^2> magnetization squared

				S_old = S;
				}
			else{
				Energy_i =energy(S_old,E,L);		 //Momentary energy
				EnergySum += Energy_i; 				 // N*<E>
				EnergySum_sq += Energy_i*Energy_i; 	 //N*<E^2> 

				M_i = magnet(S_old,L);				 
				M+= M_i;							 //<|M|> mean magnetization
				M_sq += M_i*M_i;					 //<M^2> magnetization squared
				}
			}
		else{
			Energy_i =energy(S,E,L);
			EnergySum += Energy_i; 				 // N*<E> energy
			EnergySum_sq += Energy_i*Energy_i; 	 //N*<E^2> 

			M_i = magnet(S,L);				 
			M+= M_i;							 //<|M|> mean magnetization
			M_sq += M_i*M_i;					 //<M^2> magnetization squared

			S_old = S;
			}
		}

	E_mean = EnergySum / (1.+maxMCC);
	E_mean_sq = EnergySum_sq / (1.+maxMCC);
	M_mean = M/(1.+maxMCC);
	M_mean_sq = M_sq/(1.+maxMCC);

	cout << endl;
	cout<<"Numeric Energy = "<<setprecision(15) <<E_mean << endl;


/*
	cout<<"Exact Energy = " <<setprecision(15) <<E_mean_exact << endl;
	cout<<"Numeric Expected Absolute Magnetization = " <<setprecision(15) <<M*1./maxMCC << endl;
	cout<<"Exact Expected Absolute Magnetization = " <<setprecision(15) <<M_mean_exact << endl;
	cout <<"Numeric Cv = " << setprecision(15) << (E_mean_sq - E_mean*E_mean)/T << endl;
	cout << "Exact Cv = " << setprecision(15) << Cv_exact << endl;
	cout << "Numeric chi = " << setprecision(15) << (M_mean_sq - M_mean*M_mean) << endl;
	cout << "Exact chi = " << setprecision(15) << chi_exact << endl; 
*/

	//write numeric energy and magnetization to file for different values of maxMCC:
	string outfile; //write to this file
	outfile = "energy_vs_maxMCC.txt";
	ofile.open(outfile);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "Mean Energy:       Mean Magnetization:        maxMCC:  " << omega << endl;
	for (int i=0; i<n-1; i++) {
		ofile << setw(10) << setprecision(8) << E_mean;
		ofile << setprecision(8) << M_mean << setprecision(8) <<  maxMCC << endl;
	}

	ofile.close();

	return 0;

}//end of main



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
			if (rand() < RAND_MAX/2.){
				S(i,j) *= -1;
			}
		}
		S(i,L) = S(i,0);
		S(L,i) = S(0,i);
	}

}
