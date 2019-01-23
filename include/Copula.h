#pragma once

#include <string>
#include <vector>
#include "Vector.h"
#include "Matrix.h"

enum copulaTypes { gaussian };

class Copula
{
public:
	string copula;
	vector<double> parameters;
	double c(vector <double> u, double k, CMatrix &distanceMatrix);
	double static logc(vector <double> u, double k, CMatrix &distanceMatrix);
//	CMatrix CCopula::copula_matrix(int nx, int ny);
//	CMatrix CCopula::copula_matrix_int(int nx, int ny);
	Copula();
	~Copula();
	copulaTypes copulaType = gaussian;
};

