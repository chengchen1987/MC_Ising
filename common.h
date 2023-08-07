#pragma once

#pragma warning(disable : 4996)

#include <boost/dynamic_bitset.hpp>

// parse input 
void GetParaFromInput_int(const char* fname, const char* string_match, int& para);
void GetParaFromInput_real(const char* fname, const char* string_match, double& para);
void GetParaFromInput_char(const char* fname, const char* string_match, char& para);

void Vec_fwrite_double(const char* fname, double* data, const int& dsize);
void Vec_fread_double(const char* fname, double* data, const int& dsize);

double vec_mean(int len, double* dat);
double vec_mean(int len, int* dat);
double vec_var(int len, double* dat, double avg);
double vec_std(int len, double* dat, double avg);

// dynamic_bit set functions
