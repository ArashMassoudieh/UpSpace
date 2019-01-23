
#pragma once
#include<vector>
#include<string>
#include "Vector.h"
#include "Matrix.h"

using namespace std;

class CDistribution
{
public:
	CDistribution(void);
	CDistribution(string name);
	~CDistribution(void);
	vector<double> params;
	string name;
	double evaluate(double x);
	double evaluate_CDF(double x);
	double pi;
	int n;
	vector<double> s;
	vector<double> e;
	CDistribution(int nn);
	CDistribution(const CDistribution &C);
	CDistribution operator = (const CDistribution &C);
	int GetRand();
	double inverseCDF(double u);
	double map_normal_to(double z);
	
};

//double erf(double x);
//double erfc(double x);
double Gammapdf(double x, double k, double theta);
double gamma(double x);
double unitrandom();
double getstdnormalrand();
double getnormalrand(double mu, double std);
CVector getnormal(CVector &mu, CMatrix &sigma);
CMatrix getnormal(int m, int n, double mu, double std);
CVector getnormal(int m, double mu, double std);
double getlognormalrand(double mu, double std);
CVector getlognormal(CVector &mu, CMatrix &sigma);
CMatrix getlognormal(int m, int n, double mu, double std);
CVector getlognormal(int m, double mu, double std);
double stdnormal_cdf(double u);
double GetRndUniF(double xmin, double xmax);
double std_normal_phi_inv(double u);