#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <fstream>

using namespace std;

#include "common.h"

void GetParaFromInput_int(const char* fname, const char* string_match, int& para) {
	FILE* f_in = fopen(fname, "r");
	char testchar[40], line[80];
	sprintf(testchar, string_match);
	int len = strlen(testchar);
	while (fgets(line, 200, f_in) != NULL) {
		if (strncmp(line, testchar, len) == 0)
		{
			char* p = strtok(line, "=");
			stringstream ss;
			p = strtok(NULL, "=");
			ss << p;
			ss >> para;
			std::cout << "GetParaFromInput: " << string_match << " " << para << std::endl;
			break;
		}
	}
	fclose(f_in);
	f_in = NULL;
}

void GetParaFromInput_real(const char* fname, const char* string_match, double& para) {
	FILE* f_in = fopen(fname, "r");
	char testchar[40], line[80];
	sprintf(testchar, string_match);
	int len = strlen(testchar);
	while (fgets(line, 200, f_in) != NULL) {
		if (strncmp(line, testchar, len) == 0)
		{
			char* p = strtok(line, "=");
			stringstream ss;
			p = strtok(NULL, "=");
			ss << p;
			ss >> para;
			std::cout << "GetParaFromInput: " << string_match << " " << para << std::endl;
			break;
		}
	}
	fclose(f_in);
	f_in = NULL;
}

void GetParaFromInput_char(const char* fname, const char* string_match, char& para) {
	FILE* f_in = fopen(fname, "r");
	char testchar[40], line[80];
	sprintf(testchar, string_match);
	int len = strlen(testchar);
	while (fgets(line, 200, f_in) != NULL) {
		if (strncmp(line, testchar, len) == 0)
		{
			char* p = strtok(line, "=");
			stringstream ss;
			p = strtok(NULL, "=");
			ss << p;
			ss >> para;
			std::cout << "GetParaFromInput: " << string_match << " " << para << std::endl;
			break;
		}
	}
	fclose(f_in);
	f_in = NULL;
}

void Vec_fwrite_double(const char* fname, double* data, const int& dsize)
{
	FILE* f_out;
	f_out = fopen(fname, "wb");
	fwrite(data, sizeof(double), dsize, f_out);
	fclose(f_out);
}

void Vec_fread_double(const char* fname, double* data, const int& dsize)
{
	FILE* f_in;
	f_in = fopen(fname, "rb");
	fread(data, sizeof(double), dsize, f_in);
	fclose(f_in);
}

double vec_mean(int len, double* dat)
{
	double aux = 0;
	for (int i = 0; i < len; i++)
	{
		aux += dat[i];
	}
	return aux / len;
}

double vec_mean(int len, int* dat)
{
	double aux = 0;
	for (int i = 0; i < len; i++)
	{
		aux += dat[i];
	}
	return aux / (double)len;
}

double vec_var(int len, double* dat, double avg)
{
	double aux = 0;
	for (int i = 0; i < len; i++)
	{
		aux += (dat[i] - avg) * (dat[i] - avg);
	}
	return aux;
}

double vec_std(int len, double* dat, double avg)
{
	double aux = 0;
	for (int i = 0; i < len; i++)
	{
		aux += (dat[i] - avg) * (dat[i] - avg);
	}
	return sqrt(aux / len);
}