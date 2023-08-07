#pragma once

#include <random>

class Square_Lattice
{
private:
	int Lx, Ly;
public:
	Square_Lattice(int lx, int ly);

	int get_Lx(void) { return Lx; }
	int get_Ly(void) { return Ly; }
	void NN_sites(int site, int* nn_sites);
	void print_spin(int* spins);
};

// we print basic measurements of each 'bin', and post-compute observables   
class MC_Measurements
{
public:
	double E, E2, M, M2, M4;
};

class MC_Ising
{
private:
	Square_Lattice* lattice;

public:
	MC_Ising(Square_Lattice* _lattice);
	~MC_Ising();

	// lattice information 
	int Lx, Ly, N;

	// MC parameters
	int n_warmup, n_measure, n_bins;
	double T_min, T_inc, n_T;
	double* T_vec;
	
	// MC functions 
	void MC_warmup(int iT, double T, int* spin);
	void MC_measure(int iT, double T, int* spin);
	
	void MC_step_measure(double T, int *spin, MC_Measurements &obs);
	double MC_step_Metropolis(double T, int* spin);
	double calc_H(int* spin);
	double calc_M(int* spin);
	void cal_mean_measurements(int dlen, MC_Measurements *step_obs, MC_Measurements &obs);

	// Intrinsic Dimension
	int n_id_samples;	// no. id samples in each bin



	std::random_device rd;  // Will be used to obtain a seed for the random number engine
};
