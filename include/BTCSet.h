#pragma once
#include "BTC.h"
#include <vector>
#include "Vector.h"
#include "2DMap.h"

#define CBTCSet CTimeSeriesSet

class CTimeSeriesSet
{
public:
	CTimeSeriesSet(void);
	CTimeSeriesSet(int n);
	CTimeSeriesSet(int numberofBTCs, int sizeofBTCvector);
	CTimeSeriesSet(const CTimeSeriesSet &BTC);
	CTimeSeriesSet(const CBTC &BTC);
	CTimeSeriesSet(string filename, bool varytime);
	int nvars;
	vector<CTimeSeries> BTC;
	void writetofile(char outputfile[]);
	int maxnumpoints();
	CTimeSeriesSet& operator = (const CTimeSeriesSet &C);
	vector<string> names;
	bool unif;
	void writetofile(string outputfile, bool writeColumnHeaders = false);
	void writetofile(string outputfile, int writeinterval);
	vector<double> interpolate(double t);
	vector<double> interpolate(double t, int fnvars);
	void getfromfile(string filename, bool varytime);
	double maxtime();
	double mintime();
	vector<double> getrandom();
	vector<double> percentile(double x);
	vector<double> mean(int limit);
	vector<double> mean(int limit, vector<int> index);
	vector<double> std(int limit);
	vector<double> std(int limit, vector<int> index);
	CMatrix correlation(int limit, int n);
	vector<double> integrate();
	vector<double> average();
	vector<double> percentile(double x, int limit);
	vector<double> percentile(double x, int limit, vector<int> index);
	vector<double> getrandom(int burnin);
	void append(double t, vector<double> c, double weight=1);
	CBTC add(vector<int> ii);
	CBTC add_mult(vector<int> ii, vector<double> mult);
	CBTC add_mult(vector<int> ii, CTimeSeriesSet &mult);
	CBTC divide(int ii, int jj);
	CTimeSeriesSet make_uniform(double increment);
	CTimeSeriesSet getpercentiles(vector<double> percents);
	CVector out_of_limit(double limit);
	CTimeSeriesSet distribution(int n_bins, int n_columns, int limit=0, double smoothing_factor=0);
	CTimeSeriesSet add_noise(vector<double> std, bool logd);
	void clear();
	vector<double> max_wiggle();
	vector<double> max_wiggle_corr(int _n = 10);
	vector<int> max_wiggle_sl(int ii, double tol);
	CTimeSeriesSet getflow (double A);
	void knockout(double t);
	int lookup(string S);
	vector<double> getrow(int a);
	void setname(int i, string name);
    void append(string name, double time, double c, double weight=1);
	CVector get_kappa_gamma(double delta_x);


	bool file_not_found=false;
	CBTC &operator[](int index);
	CBTC &operator[](string BTCName);

	CTimeSeriesSet(vector < vector<double> > &, int writeInterval = 1);
	int indexOf(const string& name) const;
	void pushBackName(string name);
	void append(CBTC &BTC, string name = "");
	CTimeSeriesSet sort(int burnOut = 0);
	CTimeSeriesSet detivative();
	~CTimeSeriesSet(void);
	TDMap get2DMap(int number_of_bins_x, int number_of_bins_y, double low_lim_x, double up_lim_x, double low_lim_y, double up_lim_y);
	TDMap get2DMap(int number_of_bins, double low_lim, double up_lim);
	TDMap get2DMap(int number_of_bins);
};

double diff(CTimeSeriesSet B1, CTimeSeriesSet B2);
CTimeSeriesSet operator * (const CTimeSeriesSet &BTC, const double &C);
CVector norm2dif(CTimeSeriesSet &A, CTimeSeriesSet &B);
CTimeSeriesSet merge(CTimeSeriesSet A, const CTimeSeriesSet &B);
CTimeSeriesSet merge(vector<CTimeSeriesSet> &A);
CVector sum_interpolate(vector<CTimeSeriesSet> &BTC, double t);
double sum_interpolate(vector<CTimeSeriesSet> &BTC, double t, string name);
int max_n_vars(vector<CTimeSeriesSet> &BTC);
