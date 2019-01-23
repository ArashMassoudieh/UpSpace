
#pragma once


#include <string>
#include <vector>
#include "QuickSort.h"
#include "NormalDist.h"

#define GNUplot

#define CBTC CTimeSeries



using namespace std;


class CTimeSeries
{
public:
    bool structured;
    CTimeSeries();
    CTimeSeries(int n);
    virtual ~CTimeSeries();
    int n;
    vector<double> t;
    vector<double> C;

    string name;
    string unit;
    string defaultUnit;
    vector<string> unitsList;
    vector<double> D;

    double interpol(double x);
    double interpol_D(double x);
    CTimeSeries interpol(vector<double> x);
    CTimeSeries interpol(CTimeSeries &x);
    CTimeSeries(const CTimeSeries &C);
    CTimeSeries(string Filename);
    CTimeSeries& operator = (const CTimeSeries &C);
    void readfile(string);
    void writefile(string Filename);
    double maxC();
    double minC();
    void setnumpoints(int);
    CTimeSeries Log();
    CTimeSeries XLog();
    CTimeSeries Log(double min);
    double std();
    double mean();
    double percentile(double x);
    double percentile(double x, int limit);
    double mean(int limit);
    double std(int nlimit);
    double mean_log(int limit);
    double integrate();
    double integrate(double t);
    double integrate(double t1, double t2);
    int lookupt(double t);
    double average();
    double average(double t);
    double slope(double tt);
    CTimeSeries distribution(int n_bins = 40, double smoothing_span=0, int limit=0);
    void append(double x);
    void append(double tt, double xx, double weight=1);
    void append(CTimeSeries &CC);
    CTimeSeries& operator+=(CTimeSeries &v);
    CTimeSeries& operator%=(CTimeSeries &v);
    CTimeSeries make_uniform(double increment);
    CTimeSeries extract(double t1, double t2);
    vector<double> trend();
    double mean_t();
    CTimeSeries add_noise(double std, bool);
    void assign_D();
    void clear();
    double wiggle();
    double wiggle_corr(int _n=10);
    bool wiggle_sl(double tol);
    double maxfabs();
    double max_fabs;
    void knock_out(double t);
    double AutoCor1(int i=0);
    bool file_not_found = false;
    CTimeSeries getcummulative();
    CTimeSeries getcummulative_direct(int number_of_bins);
    bool isweighted();
    CTimeSeries Exp();
    CTimeSeries fabs();
    CTimeSeries derivative();
    //GUI
    //QList <QMap <QVariant, QVariant>> CTimeSeries::compact() const;

    CTimeSeries(double a, double b, const vector<double>&x);
    CTimeSeries(double a, double b, const CTimeSeries &btc);
    CTimeSeries(const vector<double> &t, const vector<double> &C);
    CTimeSeries(vector<double>&, int writeInterval = 1);
    bool error = false;
    double rank(int i);
    CTimeSeries rank();
    CTimeSeries rank_bd(int nintervals = 100);
    CTimeSeries map_to_standard_normal(int nintervals);
    CTimeSeries uniform_cummulative(int nintervals=100);
    CTimeSeries unlog();
    CTimeSeries distribution_fw(int n_bins, double smoothing_span = 0, int limit=0, string s="");
    CTimeSeries make_flux_weighted(string log);
    CTimeSeries distribution_log(int n_bins, double smoothing_span = 0, int limit=0);
    CTimeSeries standardize();
    vector<double> weight;
    bool weighted = false;
    double autocorrelation();
    CBTC normalize_by_max();


};

double diff(CTimeSeries &BTC_p, CTimeSeries &BTC_d);
double diff_abs(CTimeSeries &BTC_p, CTimeSeries &BTC_d);
double diff_log(CTimeSeries &BTC_p, CTimeSeries &BTC_d, double lowlim);
double diff_norm(CTimeSeries &BTC_p, CTimeSeries &BTC_d);
double diff(CTimeSeries BTC_p, CTimeSeries BTC_d, int scale);
double diff(CTimeSeries BTC_p, CTimeSeries BTC_d, CTimeSeries Q);
double diff2(CTimeSeries BTC_p, CTimeSeries BTC_d);
double diff_mixed(CTimeSeries &BTC_p, CTimeSeries &BTC_d, double lowlim, double std_n, double std_ln);
double ADD(CTimeSeries &BTC_p, CTimeSeries &BTC_d);
double diff_relative(CTimeSeries &BTC_p, CTimeSeries &BTC_d, double m);
double R2(CTimeSeries BTC_p, CTimeSeries BTC_d);
double R(CTimeSeries BTC_p, CTimeSeries BTC_d, int nlimit);
CTimeSeries operator*(double, CTimeSeries);
CTimeSeries operator*(CTimeSeries, double);
CTimeSeries operator+(CTimeSeries, double);
CTimeSeries operator-(CTimeSeries, double);
CTimeSeries operator*(CTimeSeries&, CTimeSeries&);
CTimeSeries operator/(CTimeSeries, CTimeSeries);
CTimeSeries operator/(CTimeSeries BTC1, double x);
CTimeSeries operator+(CTimeSeries, CTimeSeries);
CTimeSeries operator-(CTimeSeries, CTimeSeries);
CTimeSeries operator%(CTimeSeries, CTimeSeries);
CTimeSeries operator&(CTimeSeries, CTimeSeries);
CTimeSeries operator>(CTimeSeries BTC1, CTimeSeries BTC2);
double XYbar(CTimeSeries BTC_p, CTimeSeries BTC_d);
double X2bar(CTimeSeries BTC_p, CTimeSeries BTC_d);
double Y2bar(CTimeSeries BTC_p, CTimeSeries BTC_d);
double Ybar(CTimeSeries BTC_p, CTimeSeries BTC_d);
double Xbar(CTimeSeries BTC_p, CTimeSeries BTC_d);
CTimeSeries operator+(CTimeSeries v1, CTimeSeries v2);
double prcntl(vector<double> C, double x);
vector<double> prcntl(vector<double> C, vector<double> x);
double sgn(double val);
int sgn(int val);
double sum_interpolate(vector<CTimeSeries>, double t);
double R2_c(CTimeSeries BTC_p, CTimeSeries BTC_d);
double norm2(CTimeSeries BTC1);
CTimeSeries max(CTimeSeries A, double b);
double pos(double x);
double neg(double x);
namespace TS {
    double correlation(const CTimeSeries &TS1, const CTimeSeries &TS2);
};


