#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

#include "common.h"
#include "mc_Ising.h"

// lattice geometry
Square_Lattice::Square_Lattice(int lx, int ly) :
	Lx(lx),
	Ly(ly)
{
}

void Square_Lattice::print_spin(int* spins)
{
	for (int y = 0; y < Ly; y++)
	{
		for (int x = 0; x < Lx; x++)
		{
			std::cout << std::setw(4) << spins[x + y * Lx];
		}
		std::cout << std::endl;
	}
}

void Square_Lattice::NN_sites(int site, int* nn_sites)
{
	int x = site % Lx;
	int y = site / Ly;
	// up 
	nn_sites[0] = (x + 0 + Lx) % Lx + ((y + 1 + Ly) % Ly) * Lx;
	// right 
	nn_sites[1] = (x + 1 + Lx) % Lx + ((y + 0 + Ly) % Ly) * Lx;
	// down 
	nn_sites[2] = (x + 0 + Lx) % Lx + ((y - 1 + Ly) % Ly) * Lx;
	// left  
	nn_sites[3] = (x - 1 + Lx) % Lx + ((y + 0 + Ly) % Ly) * Lx;
}

// MC functions 
MC_Ising::MC_Ising(Square_Lattice* _lattice) :
	lattice(_lattice),
	Lx(lattice->get_Lx()),
	Ly(lattice->get_Ly())
{
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> ran_0or1(0, 1);

	// lattice parameters 
	N = Lx * Ly;

	// Read Monte Carlo parameters from file "input.in"
	GetParaFromInput_int("input.in", "n_warmup", n_warmup);
	GetParaFromInput_int("input.in", "n_measure", n_measure);
	GetParaFromInput_int("input.in", "n_bins", n_bins);
	GetParaFromInput_real("input.in", "T_min", T_min);
	GetParaFromInput_real("input.in", "T_inc", T_inc);
	GetParaFromInput_real("input.in", "n_T", n_T);

	// run MC simulations on different temperatures
	T_vec = new double[n_T];
	for (int iT = 0; iT < n_T; iT++) T_vec[iT] = T_min + T_inc * iT;
	ofstream of1("T_vec.txt");
	for (int iT = 0; iT < n_T; iT++)
	{
		of1 << setw(8) << iT << setw(24) << setprecision(14) << T_vec[iT] << endl;
	}
	of1.close();

	// loop of temperature
	for (int iT = 0; iT < n_T; iT++)
	{
		cout << endl << "iT = " << iT << ", T = " << T_vec[iT] << endl;
		double T = T_vec[iT];
		
		// random initial spin configuration
		int* Spin = new int[N];
		for (int i = 0; i < N; i++) { Spin[i] = ran_0or1(gen) * 2 - 1; }
		// 
		cout << "MC_warmup..." << endl;
		size_t time0 = time(0);
		MC_warmup(iT, T, Spin);
		size_t time1 = time(0);
		cout << "  time cost(s):" << time1 - time0 << endl;
		// 
		cout << "MC_measure..." << endl;
		MC_measure(iT, T, Spin);
		size_t time2 = time(0);
		cout << "  time cost(s):" << time2 - time1 << endl;

		delete[]Spin;
	}

	// 
	delete[]T_vec;
}

MC_Ising::~MC_Ising()
{}

void MC_Ising::MC_warmup(int iT, double T, int* spin)
{
	// in MC_warmup, we save energy of each MC step to check convergence
	double* E_vec = new double[n_warmup];
	for (int im = 0; im < n_warmup; im++)
	{
		E_vec[im] = MC_step_Metropolis(T, spin);
	}
	char fname[80];
	sprintf(fname, "MC_warmup_E_iT%d.bin", iT);
	Vec_fwrite_double(fname, E_vec, n_warmup);
	delete[]E_vec;
}

void MC_Ising::MC_measure(int iT, double T, int* spins)
{
	MC_Measurements* bin_obs = new MC_Measurements[n_bins];
	for (int ibin = 0; ibin < n_bins; ibin++)
	{
		// print/store all configurations for haming distance related calculations 

		MC_Measurements* step_obs = new MC_Measurements[n_measure];
		for (int im = 0; im < n_measure; im++)
		{
			MC_step_Metropolis(T, spins);
			MC_step_measure(T, spins, step_obs[im]);
		}
		cal_mean_measurements(n_measure, step_obs, bin_obs[ibin]);
		delete[]step_obs;
	}
	// print bin_obs, and post-calculate other observables
	char ofname[80];
	sprintf(ofname, "bins_measures_E_E2_M_M2_M4_iT%d.txt", iT);
	ofstream of0(ofname);
	for (int ibin = 0; ibin < n_bins; ibin++)
	{
		of0 << setw(24) << setprecision(14) << bin_obs[ibin].E;
		of0 << setw(24) << setprecision(14) << bin_obs[ibin].E2;
		of0 << setw(24) << setprecision(14) << bin_obs[ibin].M;
		of0 << setw(24) << setprecision(14) << bin_obs[ibin].M2;
		of0 << setw(24) << setprecision(14) << bin_obs[ibin].M4;
		of0 << endl;
	}
	of0.close();

	delete[]bin_obs;
}

void MC_Ising::MC_step_measure(double T, int* spin, MC_Measurements& obs)
{
	obs.E = calc_H(spin);
	obs.E2 = obs.E * obs.E;
	obs.M = abs(calc_M(spin));
	obs.M2 = obs.M * obs.M;
	obs.M4 = obs.M2 * obs.M2;
}


double MC_Ising::MC_step_Metropolis(double T, int* spin)
{
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> ran_double_01(0, 1);

	for (int i = 0; i < N; i++)
	{
		int* NN_sites = new int[4];
		lattice->NN_sites(i, NN_sites);
		// flip energy 
		double delta_E = 2 * spin[i] * (spin[NN_sites[0]] + spin[NN_sites[1]]
			+ spin[NN_sites[2]] + spin[NN_sites[3]]);
		if (ran_double_01(gen) < exp(-delta_E / T))
		{
			spin[i] = -spin[i];
		}
		delete[]NN_sites;
	}
	return calc_H(spin);
}

double MC_Ising::calc_H(int* spin)
{
	int N = Lx * Ly;
	double aux = 0;
	for (int i = 0; i < N; i++)
	{
		int* NN_sites = new int[4];
		lattice->NN_sites(i, NN_sites);
		// count 2 bonds: up and right
		aux += -spin[i] * (spin[NN_sites[0]] + spin[NN_sites[1]]);
		delete[]NN_sites;
	}
	return aux / N;
}

double MC_Ising::calc_M(int* spin)
{
	int N = Lx * Ly;
	double aux = 0;
	for (int i = 0; i < N; i++)
	{
		aux += spin[i];
	}
	return (double)aux / (double)N;
}

void MC_Ising::cal_mean_measurements(int dlen, MC_Measurements* step_obs, MC_Measurements& obs)
{
	obs.E = 0;
	obs.E2 = 0;
	obs.M = 0;
	obs.M2 = 0;
	obs.M4 = 0;
	for (int i = 0; i < dlen; i++)
	{
		obs.E += step_obs[i].E;
		obs.E2 += step_obs[i].E2;
		obs.M += step_obs[i].M;
		obs.M2 += step_obs[i].M2;
		obs.M4 += step_obs[i].M4;
	}
	obs.E = obs.E / dlen;
	obs.E2 = obs.E2 / dlen;
	obs.M = obs.M / dlen;
	obs.M2 = obs.M2 / dlen;
	obs.M4 = obs.M4 / dlen;
}