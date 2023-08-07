#include <cstdio>
#include <iostream>
#include <random>
#include <iomanip>

using namespace std;

#include "common.h"
#include "mc_Ising.h"

// Classical Monte Carlo for 2D Ising model

int main()
{
	size_t time0= time(0);

	int Lx, Ly;
	GetParaFromInput_int("input.in", "Lx", Lx);
	GetParaFromInput_int("input.in", "Ly", Ly);

	cout << "==================================================================" << endl;
	cout << "Monte Carlo simulations for the 2D Ising model on square lattices." << endl;
	cout << "Lx = " << Lx << endl;
	cout << "Ly = " << Ly << endl;

	// generate basis
	Square_Lattice lattice(Lx,Ly);
	MC_Ising mc(&lattice);

	size_t time1 = time(0);
	cout << endl << "==============================================================" << endl;
	cout << "Total time cost(s):" << time1 - time0 << endl;

	return 0;
}