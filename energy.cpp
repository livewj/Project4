#include <iostream>
#include <armadillo>
#include <cstdlib>
using namespace arma;
using namespace std;

void energy(mat& S, double E, double M, int L, int iter);
void ising_create(mat& S, double L);

int main(){
	int L = 2;
	int maxMCC = 10000; //Number of ising models created/MCsycles
	mat S = mat(L+1,L+1,fill::ones);


	double E = 0; //Total energy
	double M = 0; //Total magnetization
	double E_mean = 0; //Mean energy
	double EnergySum = 0;
	int iter = 0; //count number of MC cycles
	double k = 2;
	ising_create(S,L); //create ising matrix
    energy(S,E,M,L,iter); //compute total spin energy

/*
mat S_prev = ising_create(S,L); //initial matrix

	for (int MCC = 0; MCC < maxMCC; MCC++) {  //LES 435 i kompendiet
		mat S = mat(L+1,L+1,fill::ones);  //generate spin-matrix
		ising_create(S,L);


		// Foreslå spin flip     <-- METROPOLIS
		// Accept kanskje
		// Hvis accept:
		EnergySum += energy(newState);
		// Hvis ikke:
		EnergySum += energy(oldState);
	}
	double EnergyExpectationValue = EnergySum / maxMCC;
	cout << EnergyExvalue << endl;
*/
}


void energy(mat& S, double E, double M, int L, int iter) {

	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			E += S(i,j)*S(i,j+1); //neighbour energy (beneath)
			E += S(i,j)*S(i+1,j); //neighbour energy (besides)
			M += S(i,j);		  //total magnetization
			iter += 1;
		}
	}
cout << "Energy: " << E << endl;
cout << "Mean magnetization: " << abs(M)/iter << endl;
cout << "# of MC cycles: " << iter << endl;
}

void ising_create(mat& S, double L) {
    //initialize random seed
    srand (time(NULL));

	for (int i = 0; i < L+1; i++) {
		for (int j = 0; j < L+1; j++) {
			if (rand() < RAND_MAX/2.) {  //RAND_MAX = largest number on the computer
			    S(i,j) = -1; 
			}
			//periodic boundary conditions
		    S(L,j) = S(0,j); 
		}
		S(i,L) = S(i,0); 
	} 
	cout << S << endl;
}

