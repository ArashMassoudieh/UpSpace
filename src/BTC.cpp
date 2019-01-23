// BTC.cpp: implementation of the CTimeSeries class.
//
//////////////////////////////////////////////////////////////////////

#include "BTC.h"
#include "math.h"
#include "string.h"
#include <iostream>
#include <fstream>
#include "StringOP.h"
#include "NormalDist.h"
#ifdef Qt_version
#include <qdebug.h>
#endif // Qt_version
#define GNUplot

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTimeSeries::CTimeSeries()
{
	n = 0;
	structured = true;
	max_fabs = 0;

}

CTimeSeries::CTimeSeries(int n1)
{
	n=n1;
	t.resize(n);
	C.resize(n);
        weight.resize(n);
	structured = true;
	max_fabs = 0;
}

CTimeSeries::CTimeSeries(vector<double> &data, int writeInterval)
{
    n = 0;
    structured = 0;
    for (int i = 0; i < int(data.size()); i++)
	if (i%writeInterval == 0)
	{
            n++;
            t.push_back(i);
            C.push_back(data[i]);
            weight.push_back(1.0);
	}
}
double CTimeSeries::rank(int ii)
{
	int xx = 0;
	for (int i = 0; i < n; i++)
		if (C[i] <= ii) xx++;
	return double(xx / double(n));

}
CTimeSeries CTimeSeries::rank()
{
	CTimeSeries out;
	for (int i = 0; i < n; i++)
		out.append(t[i], rank(i));

	return out;
}
CTimeSeries CTimeSeries::rank_bd(int nintervals)
{
	CTimeSeries X = getcummulative_direct(nintervals).XLog();
    X.structured = true;

	CTimeSeries out;
        out.weighted = true;
	if (weight.size())
            for (int i = 0; i < n; i++)
            {
                out.append(double(i), X.interpol(log(C[i])),weight[i]);
            }
        else
            for (int i = 0; i < n; i++)
            {
                out.append(double(i), X.interpol(log(C[i])),1);
            }
	return out;
}
CTimeSeries CTimeSeries::map_to_standard_normal(int nintervals)
{
    CTimeSeries dist = distribution(nintervals, 0);
    CTimeSeries X = dist.uniform_cummulative(nintervals);
    CTimeSeries out;
    out.weighted = true;
    if (weight.size())
        for (int i = 0; i < n; i++)
            out.append(double(i), stdnormal_inv(X.interpol(C[i])),weight[i]);
    else
        for (int i = 0; i < n; i++)
            out.append(double(i), stdnormal_inv(X.interpol(C[i])),1);

    return out.standardize();
}

CTimeSeries CTimeSeries::uniform_cummulative(int nintervals)
{
    CTimeSeries X = getcummulative();
    return X.make_uniform(t[n - 1] / double(nintervals));

}
void CTimeSeries::setnumpoints(int n1)
{
    n = n1;
    t.resize(n);
    C.resize(n);
    weight.resize(n);
}

CTimeSeries::~CTimeSeries()
{

}

CTimeSeries::CTimeSeries(const CTimeSeries &CC)
{
	n=CC.n;
	t = CC.t;
	D = CC.D;
	C = CC.C;
        weight = CC.weight;
        weighted = CC.weighted;
	structured = CC.structured;
	name = CC.name;
	unit = CC.unit;
	defaultUnit = CC.defaultUnit;
	unitsList = CC.unitsList;
	error = CC.error;
}

CTimeSeries::CTimeSeries(string Filename)
{
	n = 0;
	t.clear();
	C.clear();
	ifstream file(Filename.c_str());
	if (file.good() == false)
	{
            file_not_found = true;
            error = true;
            return;
	}

	vector<string> s;
	structured = true;
	if (file.good())
	while (file.eof()== false)
	{
		s = getline(file);
		if (s.size() == 1)
		{
			error = true;
//			file.close();
//			return;
		}
		if (s.size()>=2)
		if ((s[0].substr(0,2)!="//") && (tolower(s[0])!="names"))
		{
			t.push_back(atof(s[0].c_str()));
			C.push_back(atof(s[1].c_str()));
                        weight.push_back(1.0);
			n++;
			if (t.size()>2)
				if ((t[t.size()-1]-t[t.size()-2])!=(t[t.size()-2]-t[t.size()-3]))
					structured = false;

		}
	}
	error = (n == 0) ? true : false;
	file.close();
}

/*CTimeSeries CTimeSeries::operator = (const CTimeSeries &CC)
{
	n=CC.n;
	t = new double[n];
	C = new double[n];
	for (int i=0; i<n; i++)
	{
		t[i] = CC.t[i];
		C[i] = CC.C[i];
	}

	return *this;
}*/

CTimeSeries& CTimeSeries::operator = (const CTimeSeries &CC)
{
	n=CC.n;
	t = CC.t;
	D = CC.D;
	C = CC.C;
        weight = CC.weight;
        weighted = CC.weighted;
	structured = CC.structured;
	name = CC.name;
	unit = CC.unit;
	defaultUnit = CC.defaultUnit;
	unitsList = CC.unitsList;
	error = CC.error;

	return *this;
}

CTimeSeries CTimeSeries::Log()
{
    CTimeSeries _BTC;
    _BTC.weighted = weighted;
    for (int i=0; i<n; i++)
    {
        if (C[i]>0)
           _BTC.append(t[i],log(C[i]));

    }
    return _BTC;
}

CTimeSeries CTimeSeries::XLog()
{
    CTimeSeries _BTC;
    _BTC.weighted = weighted;
    for (int i=0; i<n; i++)
    {
        if (C[i]>0)
           _BTC.append(log(t[i]),C[i]);

    }
    return _BTC;
}


CTimeSeries CTimeSeries::Log(double m)
{
	CTimeSeries _BTC(n);
        _BTC.weighted = weighted;
	for (int i=0; i<n; i++)
	{
		_BTC.t[i] = t[i];
		_BTC.C[i] = log(max(C[i],m));
	}
	return _BTC;
}


double CTimeSeries::interpol(double x)
{
    double r=0;
    if (n>1)
    {

        if (structured == false)
        {	for (int i=0; i<n-1; i++)
            {
                    if (t[i] <= x && t[i+1] >= x)
                            r=(C[i+1]-C[i])/(t[i+1]-t[i])*(x-t[i]) + C[i];
            }
            if (x>t[n-1]) r=C[n-1];
            if (x<t[0]) r=C[0];
        }
        else
        {
            if (x < t[0]) return C[0];
            if (x > t[n - 1]) return C[n - 1];
            double dt = t[1]-t[0];
            int i = int((x-t[0])/dt);
            if (i>=n-1) r=C[n-1];
            else if (i<0) r=C[0];
            else r=(C[i+1]-C[i])/(t[i+1]-t[i])*(x-t[i]) + C[i];
        }
    }
    else
            r = C[0];
    return r;

}

double CTimeSeries::interpol_D(double x)
{
    double r=0;
    if (n>1)
    {

        if (structured == false)
        {	for (int i=0; i<n-1; i++)
            {
                if (t[i] <= x && t[i+1] >= x)
                        r=(D[i+1]-D[i])/(t[i+1]-t[i])*(x-t[i]) + D[i];
            }
            if (x>t[n-1]) r=D[n-1];
            if (x<t[0]) r=D[0];
        }
        else
        {
            double dt = t[1]-t[0];
            int i = int((x-t[0])/dt);
            if (x>=t[n-1]) r=D[n-1];
            else if (x<=t[0]) r=D[0];
            else r=(D[i+1]-D[i])/(t[i+1]-t[i])*(x-t[i]) + D[i];
        }
    }
    else
            r = D[0];
    return r;

}

CTimeSeries CTimeSeries::interpol(vector<double> x)
{
    CTimeSeries BTCout;
    for (int i=0; i<int(x.size()); i++)
        BTCout.append(x[i],interpol(x[i]));
    return BTCout;

}

CTimeSeries CTimeSeries::interpol(CTimeSeries &x)
{
    CTimeSeries BTCout;
    for (int i=0; i<x.n; i++)
            BTCout.append(x.t[i],interpol(x.t[i]));
    return BTCout;

}

double ADD(CTimeSeries &BTC_p, CTimeSeries &BTC_d)
{
    double sum = 0;
    for (int i=0; i<BTC_d.n; i++)
            if (abs(BTC_d.C[i]) < 1e-3)
                    sum += abs(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]));
            else
                    sum += abs(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i])) /BTC_d.C[i];

    return sum/BTC_d.n;
}

double diff_relative(CTimeSeries &BTC_A, CTimeSeries &BTC_B, double m)
{
    double sum = 0;
    for (int i=0; i<min(BTC_A.n,BTC_B.n); i++)
        if (abs(BTC_A.C[i]) < m)
                sum += abs(BTC_B.C[i] - BTC_A.interpol(BTC_B.t[i]));
        else
                sum += abs(BTC_B.C[i] - BTC_A.interpol(BTC_B.t[i])) /BTC_A.C[i];

    return sum;
}


double diff(CTimeSeries BTC_p, CTimeSeries BTC_d, int scale)
{
    double sum = 0;
    for (int i=0; i<BTC_d.n; i++)
    {
        if (BTC_d.C[i] > BTC_p.interpol(BTC_d.t[i]))
                sum += scale*pow(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]),2)/sqrt(1.0+scale*scale);
        else
                sum += pow(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]),2)/sqrt(1.0+scale*scale);
    }
    return sum;
}

double diff(CTimeSeries &BTC_p, CTimeSeries &BTC_d)
{
    double sum = 0;
    double a;
    for (int i=0; i<BTC_d.n; i++)
    {
        a = BTC_p.interpol(BTC_d.t[i]);
        sum += pow(BTC_d.C[i] - a,2);
    }

    return sum;
}

double diff_abs(CTimeSeries &BTC_p, CTimeSeries &BTC_d)
{
    double sum = 0;

    for (int i=0; i<BTC_d.n; i++)
    {
            sum += abs(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]));
    }

    return sum;
}

double diff_log(CTimeSeries &BTC_p, CTimeSeries &BTC_d, double lowlim)
{
    double sum = 0;
    double a;
    for (int i=0; i<BTC_d.n; i++)
    {
        a = BTC_p.interpol(BTC_d.t[i]);
        sum += pow(log(max(BTC_d.C[i],lowlim)) - log(max(a,lowlim)),2);

    }

    return sum;
}


double diff2(CTimeSeries BTC_p, CTimeSeries BTC_d)
{
    double sum = 0;
    double sumvar1 = 0;
    for (int i=0; i<BTC_d.n; i++)
    {
        sum += pow(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]),2);
        sumvar1 += BTC_d.C[i]*BTC_d.C[i];
    }

    return sum;
}


double R2(CTimeSeries BTC_p, CTimeSeries BTC_d)
{
    double sumcov = 0;
    double sumvar1 = 0;
    double sumvar2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    for (int i=0; i<BTC_d.n; i++)
    {
        double x2 = BTC_p.interpol(BTC_d.t[i]);
        sumcov += BTC_d.C[i]*x2/BTC_d.n;
        sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
        sumvar2 += x2*x2/BTC_d.n;
        sum1 += BTC_d.C[i]/BTC_d.n;
        sum2 += x2/BTC_d.n;
    }

    return pow(sumcov-sum1*sum2,2)/(sumvar1-sum1*sum1)/(sumvar2-sum2*sum2);
}

double R(CTimeSeries BTC_p, CTimeSeries BTC_d, int nlimit)
{
    double sumcov = 0;
    double sumvar1 = 0;
    double sumvar2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    int N = BTC_d.n - nlimit;

    for (int i=nlimit; i<BTC_d.n; i++)
    {
        double x1 = BTC_d.C[i];
        double x2 = BTC_p.C[i];
        sumcov += x1*x2/N;
        sumvar1 += x1*x1/N;
        sumvar2 += x2*x2/N;
        sum1 += x1/N;
        sum2 += x2/N;
    }

    double R_x1x2 = (sumcov-sum1*sum2)/pow(sumvar1-sum1*sum1,0.5)/pow(sumvar2-sum2*sum2,0.5);

    return R_x1x2;
}

double XYbar(CTimeSeries BTC_p, CTimeSeries BTC_d)
{
    double sumcov = 0;
    double sumvar1 = 0;
    double sumvar2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    for (int i=0; i<BTC_d.n; i++)
    {
        double x2 = BTC_p.interpol(BTC_d.t[i]);
        sumcov += BTC_d.C[i]*x2/BTC_d.n;
        sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
        sumvar2 += x2*x2/BTC_d.n;
        sum1 += BTC_d.C[i]/BTC_d.n;
        sum2 += x2/BTC_d.n;
    }

    return sumcov;
}

double X2bar(CTimeSeries BTC_p, CTimeSeries BTC_d)
{
    double sumcov = 0;
    double sumvar1 = 0;
    double sumvar2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    for (int i=0; i<BTC_d.n; i++)
    {
        double x2 = BTC_p.interpol(BTC_d.t[i]);
        sumcov += BTC_d.C[i]*x2/BTC_d.n;
        sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
        sumvar2 += x2*x2/BTC_d.n;
        sum1 += BTC_d.C[i]/BTC_d.n;
        sum2 += x2/BTC_d.n;
    }

    return sumvar1;
}

double Y2bar(CTimeSeries BTC_p, CTimeSeries BTC_d)
{
    double sumcov = 0;
    double sumvar1 = 0;
    double sumvar2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    for (int i=0; i<BTC_d.n; i++)
    {
        double x2 = BTC_p.interpol(BTC_d.t[i]);
        sumcov += BTC_d.C[i]*x2/BTC_d.n;
        sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
        sumvar2 += x2*x2/BTC_d.n;
        sum1 += BTC_d.C[i]/BTC_d.n;
        sum2 += x2/BTC_d.n;
    }

    return sumvar2;
}

double Ybar(CTimeSeries BTC_p, CTimeSeries BTC_d)
{
    double sumcov = 0;
    double sumvar1 = 0;
    double sumvar2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    for (int i=0; i<BTC_d.n; i++)
    {
        double x2 = BTC_p.interpol(BTC_d.t[i]);
        sumcov += BTC_d.C[i]*x2/BTC_d.n;
        sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
        sumvar2 += x2*x2/BTC_d.n;
        sum1 += BTC_d.C[i]/BTC_d.n;
        sum2 += x2/BTC_d.n;
    }

    return sum2;
}

double Xbar(CTimeSeries BTC_p, CTimeSeries BTC_d)
{
    double sumcov = 0;
    double sumvar1 = 0;
    double sumvar2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    for (int i=0; i<BTC_d.n; i++)
    {
        double x2 = BTC_p.interpol(BTC_d.t[i]);
        sumcov += BTC_d.C[i]*x2/BTC_d.n;
        sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
        sumvar2 += x2*x2/BTC_d.n;
        sum1 += BTC_d.C[i]/BTC_d.n;
        sum2 += x2/BTC_d.n;
    }

    return sum1;
}

double diff_norm(CTimeSeries &BTC_p, CTimeSeries &BTC_d)
{
    double sum = 0;
    double sumvar1 = 0;
    double sumvar2 = 0;
    double a;
    for (int i=0; i<BTC_d.n; i++)
    {
        a = BTC_p.interpol(BTC_d.t[i]);
        sum += pow(BTC_d.C[i] - a,2)/BTC_d.n;
        sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
        sumvar2 += pow(a,2)/BTC_d.n;
    }
    //cout<<sum<<endl;
    return sum/sqrt(sumvar1*sumvar2);

}


double diff(CTimeSeries BTC_p, CTimeSeries BTC_d, CTimeSeries Q)
{
    double sum = 0;
    for (int i=0; i<BTC_d.n; i++)
    {
        sum += pow(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]),2)*pow(Q.interpol(BTC_d.t[i]),2);
    }
    return sum;
}

void CTimeSeries::readfile(string Filename)
{
    ifstream file(Filename.c_str());
    vector<string> s;
    if (file.good() == false)
    {
        file_not_found = true;
        return;
    }

    if (file.good())
    while (file.eof()== false)
    {
        s = getline(file);
        if (s.size()>0)
        if (s[0].substr(0,2)!="//")
        {
            t.push_back(atof(s[0].c_str()));
            C.push_back(atof(s[1].c_str()));
            n++;
            if (t.size()>2)
                if (t[t.size()-1]-t[t.size()-2]!=t[t.size()-2]-t[t.size()-3])
                    structured = false;

            }
    }
    file.close();

}

/*void CTimeSeries::readfile(CString Filename)
{
	FILE *FILEBTC;
	FILEBTC = fopen(Filename, "r");
	if (FILEBTC == NULL)
		double e=1;
	int numpoints = 0;
	double tt, CC;
	while (feof(FILEBTC)==false)
	{
		fscanf(FILEBTC, "%lf, %lf\n", &tt, &CC);
		numpoints++;
	}
	//numpoints--;
	fclose(FILEBTC);

	n=numpoints;
	t = new double[numpoints];
	C = new double[numpoints];

	FILEBTC = fopen(Filename, "r");
	for (int i=0; i<numpoints; i++)
	{
		fscanf(FILEBTC, "%lf, %lf", &t[i], &C[i]);
	}
	fclose(FILEBTC);



}*/

void CTimeSeries::writefile(string Filename)
{
    FILE *FILEBTC;
    FILEBTC = fopen(Filename.c_str(), "w");
#ifndef GNUplot
        fprintf(FILEBTC, "n %i, BTC size %i\n", n, int(C.size()));
	for (int i=0; i<n; i++)
		fprintf(FILEBTC, "%lf, %le\n", t[i], C[i]);
#else
    if (!weighted)
    {
        fprintf(FILEBTC, "# n %i\t  BTC size %i\n", n, int(C.size()));
        for (int i=0; i<n; i++)
                fprintf(FILEBTC, "%lf\t %le\n", t[i], C[i]);
    }
    else
    {
        fprintf(FILEBTC, "# n %i\t  BTC size %i\n", n, int(C.size()));
        for (int i=0; i<n; i++)
                fprintf(FILEBTC, "%lf\t %le \t %le\n", t[i], C[i], weight[i]);
    }
#endif
    fclose(FILEBTC);

}

/*void CTimeSeries::writefile(CString Filename)
{
	FILE *FILEBTC;
	FILEBTC = fopen(Filename, "w");
	for (int i=0; i<n; i++)
		fprintf(FILEBTC, "%lf, %le\n", t[i], C[i]);

	fclose(FILEBTC);

}*/

/*double CTimeSeries::GetS0(CTimeSeries &M)
{
	double sumprod = 0;
	double sumsqr = 0;
	for (int i = 0; i<M.n; i++)
	{
		sumprod += M.C[i]*interpol(M.t[i]);
		sumsqr += interpol(M.t[i])*interpol(M.t[i]);
	}
	double S0 = sumprod/sumsqr;
	return S0;
}

double CTimeSeries::GetS0(CTimeSeries &M, CTimeSeries &Q)
{
	double sumprod = 0;
	double sumsqr = 0;
	for (int i = 0; i<M.n; i++)
	{
		sumprod += M.C[i]*interpol(M.t[i])*pow(Q.interpol(M.t[i]),2);
		sumsqr += interpol(M.t[i])*interpol(M.t[i])*pow(Q.interpol(M.t[i]),2);
	}
	double S0 = sumprod/sumsqr;
	return S0;
}*/

CTimeSeries operator*(double alpha, CTimeSeries CTimeSeries_T)
{
    CTimeSeries S(CTimeSeries_T.n);
    for (int i=0; i<CTimeSeries_T.n; i++)
    {
        S.t[i] = CTimeSeries_T.t[i];
        S.C[i] = alpha*CTimeSeries_T.C[i];
        S.weight[i] = CTimeSeries_T.weight[i];
    }

    return S;
}

CTimeSeries operator*(CTimeSeries CTimeSeries_T, double alpha)
{
    CTimeSeries S = CTimeSeries_T;
    for (int i=0; i<CTimeSeries_T.n; i++)
    {
	//S.t[i] = CTimeSeries_T.t[i];
	S.C[i] = alpha*CTimeSeries_T.C[i];
    }


    return S;
}

CTimeSeries operator+(CTimeSeries CTimeSeries_T, double alpha)
{
    CTimeSeries S = CTimeSeries_T;
    for (int i=0; i<CTimeSeries_T.n; i++)
    {
	//S.t[i] = CTimeSeries_T.t[i];
	S.C[i] = alpha+CTimeSeries_T.C[i];
    }


    return S;
}

CTimeSeries operator-(CTimeSeries CTimeSeries_T, double alpha)
{
    CTimeSeries S = CTimeSeries_T;
    for (int i=0; i<CTimeSeries_T.n; i++)
    {
	//S.t[i] = CTimeSeries_T.t[i];
	S.C[i] = CTimeSeries_T.C[i] - alpha;
    }


    return S;
}

CTimeSeries operator/(CTimeSeries BTC1, CTimeSeries BTC2)
{
    CTimeSeries S = BTC1;
    for (int i=0; i<BTC1.n; i++)
        S.C[i] = BTC1.C[i]/BTC2.interpol(BTC1.t[i]);

    return S;

}

CTimeSeries operator/(CTimeSeries BTC1, double x)
{
    CTimeSeries S = BTC1;
    for (int i=0; i<BTC1.n; i++)
        S.C[i] = BTC1.C[i]/x;

    return S;

}

CTimeSeries operator-(CTimeSeries BTC1, CTimeSeries BTC2)
{
    CTimeSeries S = BTC1;
    for (int i=0; i<BTC1.n; i++)
        S.C[i] = BTC1.C[i]-BTC2.interpol(BTC1.t[i]);

    return S;
}


CTimeSeries operator*(CTimeSeries &BTC1, CTimeSeries &BTC2)
{
    CTimeSeries S = BTC1;
    for (int i=0; i<BTC1.n; i++)
        S.C[i] = BTC1.C[i]*BTC2.interpol(BTC1.t[i]);

    return S;
}

CTimeSeries operator%(CTimeSeries BTC1, CTimeSeries BTC2)
{
    CTimeSeries S = BTC1;
    for (int i=0; i<BTC1.n; i++)
        S.C[i] = BTC1.C[i]/BTC2.C[i];

    return S;
}
CTimeSeries operator&(CTimeSeries BTC1, CTimeSeries BTC2)
{
    CTimeSeries S = BTC1;
    for (int i=0; i<BTC1.n; i++)
        S.C[i] = BTC1.C[i]+BTC2.C[i];

    return S;


}

/*double CTimeSeries::EMC(CTimeSeries &M)
{
	double sum = 0;
	double sumflow = 0;
	for (int i=0; i<n; i++)
	{
		sum += C[i]*M.interpol(t[i]);
		sumflow += M.interpol(t[i]);
	}
	if (sumflow == 0.0)
		return 0;
	else
		return sum/sumflow;
}

double CTimeSeries::Calculate_load(CTimeSeries &M)
{
	double sum = 0;
	double sumflow = 0;
	for (int i=0; i<n; i++)
	{
		sum += C[i]*M.interpol(t[i])*(t[2]-t[1]);

	}

	return sum;
}*/

double CTimeSeries::maxC()
{
    double max = -1e32;
    for (int i=0; i<n; i++)
    {	if (C[i]>max)
            max = C[i];
    }
    return max;
}

double CTimeSeries::maxfabs()
{
    if (max_fabs>0)
            return max_fabs;
    else
    {
        double max = -1e32;
        for (int i=0; i<n; i++)
        {	if (std::fabs(C[i])>max)
                        max = std::fabs(C[i]);
        }
        return max;
    }

}

double CTimeSeries::minC()
{
    double min = 1e32;
    for (int i=0; i<n; i++)
    {	if (C[i]<min)
            min = C[i];
    }
    return min;
}

double CTimeSeries::std()
{
    double sum = 0;
    double m = mean();
    for (int i=0; i<n; i++)
        sum+= pow(C[i]-m,2);

    return sqrt(sum/n);
}

double CTimeSeries::std(int nlimit)
{
    double sum = 0;
    double m = mean(nlimit);
    for (int i=nlimit; i<n; i++)
        sum+= pow(C[i]-m,2);

    return sqrt(sum/n);
}

double CTimeSeries::mean()
{
    double sum = 0;
    for (int i=0; i<n; i++)
    {
        sum+= C[i];
    }
    if (n>0)
        return sum/n;
    else
        return 0;
}

double CTimeSeries::integrate()
{
    double sum = 0;
    for (int i=1; i<n; i++)
        sum+= (C[i]+C[i-1])/2.0*(t[i]-t[i-1]);

    return sum;
}

double CTimeSeries::integrate(double tt)
{
    double sum = 0;
    for (int i = 1; i<n; i++)
        if (t[i]<=tt) sum += (C[i] + C[i - 1]) / 2.0*(t[i] - t[i - 1]);

    return sum;
}

double CTimeSeries::integrate(double t1, double t2)
{
    if (structured)
    {
	int i1 = int(t1 - t[0]) / (t[1] - t[0]);
	int i2 = int(t1 - t[0]) / (t[1] - t[0]);
	double sum=0;
	for (int i = i1; i <= i2; i++)
            sum += C[i] / (i2+1 - i1)*(t2-t1);
            return sum;
    }
    else
    {
        return 0; // fill in later
    }
}

int CTimeSeries::lookupt(double _t)
{
	for (int i = 0; i < n - 1; i++)
		if ((t[i]<_t) && (t[i + 1]>_t))
			return i;
        return -1;
}

double CTimeSeries::average()
{
	if (n>0)
		return integrate()/(t[n-1]-t[0]);
	else
		return 0;
}

double CTimeSeries::average(double tt)
{
	if (n>0)
		return integrate(tt) / (max(tt,t[n - 1]) - t[0]);
	else
		return 0;
}

double CTimeSeries::slope(double tt)
{
   return (C[n - 1] - C[n - 2]) / (t[n - 1] - t[n - 2]);
}



double CTimeSeries::percentile(double x)
{
	vector<double> X = QSort(C);
	int i = int(x*X.size());
	return X[i];

}

double CTimeSeries::percentile(double x, int limit)
{
	vector<double> C1(C.size()-limit);
	for (int i=0; i<int(C1.size()); i++)
		C1[i] = C[i+limit];
	vector<double> X = bubbleSort(C1);
	//vector<double> X = bubbleSort(C1);
//	vector<double> X = C1;
	int ii = int(x*double(X.size()));
	return X[ii];

}

double CTimeSeries::mean(int limit)
{
	double sum = 0;
	for (int i=limit; i<n; i++)
		sum += C[i];
	return sum/double(n-limit);
}

double CTimeSeries::mean_log(int limit)
{
	double sum = 0;
	for (int i=limit; i<n; i++)
		sum += log(C[i]);
	return sum/double(n-limit);
}

void CTimeSeries::append(double x)
{
	n++;
	t.push_back(0);
	C.push_back(x);
	max_fabs = max(max_fabs,std::fabs(x));

}

void CTimeSeries::append(double tt, double xx, double _weight)
{
    n++;
    t.push_back(tt);
    C.push_back(xx);
    if (weighted)
        weight.push_back(_weight);


    if (n>2)
        if (t[n-1]-t[n-2]!=t[n-2]-t[n-3])
            structured = false;
    max_fabs = max(max_fabs,std::fabs(xx));
}
void CTimeSeries::append(CTimeSeries &CC)
{
    if (!CC.weighted)
        for (int i = 0; i<CC.n; i++) append(CC.t[i], CC.C[i]);
    else
        for (int i = 0; i<CC.n; i++) append(CC.t[i], CC.C[i], CC.weight[i]);
}

CTimeSeries& CTimeSeries::operator+=(CTimeSeries &v)
{
	for (int i=0; i<n; ++i)
		C[i] += v.interpol(t[i]);
	return *this;
}

CTimeSeries& CTimeSeries::operator%=(CTimeSeries &v)
{
	for (int i=0; i<n; ++i)
		C[i] += v.C[i];
	return *this;

}

CTimeSeries operator+(CTimeSeries v1, CTimeSeries v2)
{
	return v1 += v2;
}

CTimeSeries CTimeSeries::make_uniform(double increment)
{
	CTimeSeries out;
	assign_D();

	if (t.size() >1 && C.size() > 1)
	{
		out.append(t[0], C[0]);
		out.D.push_back(0);
		for (int i = 0; i < n - 1; i++)
		{
			int i1 = int((t[i] - t[0]) / increment);
			int i2 = int((t[i + 1] - t[0]) / increment);
			for (int j = i1 + 1; j <= i2; j++)
			{
                            double x = j*increment + t[0];
                            double CC = (x - t[i]) / (t[i + 1] - t[i])*(C[i + 1] - C[i]) + C[i];
                            double DD = (x - t[i]) / (t[i + 1] - t[i])*(D[i + 1] - D[i]) + D[i];
                            out.append(x, CC);
                            out.D.push_back(DD);

			}
		}
	}
	out.structured = true;

	return out;

}

double prcntl(vector<double> C, double x)
{
	vector<double> X = QSort(C);
	int ii = int(x*double(X.size()));
	return X[ii];

}

vector<double> prcntl(vector<double> C, vector<double> x)
{
	vector<double> X = QSort(C);
	vector<double> Xout = x;
	for(int j =0; j< int(x.size()); j++)
	{
		int ii = int(x[j]*double(X.size()));
		Xout[j] = X[ii];
	}

	return Xout;
}

CTimeSeries CTimeSeries::extract(double t1, double t2)
{
	CTimeSeries out;
	for (int i=0; i<n; i++)
		if ((t[i]>=t1) && (t[i]<=t2))
			out.append(t[i], C[i]);

	return out;
}


CTimeSeries CTimeSeries::distribution(int n_bins, double smoothing_span, int limit)
{
    CTimeSeries out(n_bins+2);

    CVector C1(C.size()-limit);
    for (int i=0; i<C1.num; i++)
            C1[i] = C[i+limit];

    double p_start = min(C1);
    double p_end = max(C1)*1.001;
    double dp = abs(p_end - p_start)/n_bins;
    //cout << "low limit: " << p_start << " up limit: " << p_end << " increment: " << dp << std::endl;
    if (dp == 0) return out;
    out.t[0] = p_start - dp/2;
    out.C[0] = 0;
    for (int i=0; i<n_bins+1; i++)
    {
        out.t[i+1] = out.t[i] + dp;
        out.C[i+1] = out.C[i];
    }

    if (smoothing_span==0)
    {
        if (!weighted)
        {   for (int i=0; i<C1.num; i++)
                out.C[int((C1[i]-p_start)/dp)+1] += 1.0/C1.num/dp;
            return out;
        }
        else
        {   for (int i=0; i<C1.num; i++)
                out.C[int((C1[i]-p_start)/dp)+1] += 1.0/C1.num/dp*weight[i];
            return out/out.integrate();
        }
    }
    else
    {
        int span_count = smoothing_span/dp;
        if (!weighted)
        {   for (int i=0; i<C1.num; i++)
            {
                int center = int((C1[i]-p_start)/dp)+1;
                for (int j=max(0,center-3*span_count); j<=min(n_bins,center+3*span_count); j++)
                {
                    double l_bracket = p_start + (j-1)*dp;
                    double r_bracket = p_start + (j)*dp;
                    double eff_smoothing_span = max(min(min(C[i]-p_start,smoothing_span),p_end-C[i]),dp/10.0);
                    double portion = (exp((C1[i]-l_bracket)/eff_smoothing_span)/(1+exp((C1[i]-l_bracket)/eff_smoothing_span)) - exp((C1[i]-r_bracket)/eff_smoothing_span)/(1+exp((C1[i]-r_bracket)/eff_smoothing_span)));
                    out.C[j] += 1.0/C1.num/dp*portion;
                }
            }
            return out/out.integrate();
        }
        else
        {   for (int i=0; i<C1.num; i++)
            {
                int center = int((C1[i]-p_start)/dp)+1;
                for (int j=max(0,center-3*span_count); j<=min(n_bins,center+3*span_count); j++)
                    {
                        double l_bracket = p_start + (j-1)*dp;
                        double r_bracket = p_start + (j)*dp;
                        double eff_smoothing_span = max(min(min(C[i]-p_start,smoothing_span),p_end-C[i]),dp/10.0);
                        double portion = (exp((C1[i]-r_bracket)/eff_smoothing_span)/(1+exp((C1[i]-r_bracket)/eff_smoothing_span)) - exp((C1[i]-l_bracket)/eff_smoothing_span)/(1+exp((C1[i]-l_bracket)/eff_smoothing_span)));
                        out.C[j] += 1.0/C1.num/dp*portion*weight[i];
                    }
            }
            return out/out.integrate();
        }
    }

}

vector<double> CTimeSeries::trend()
{
	double x_bar = mean_t();
	double y_bar = mean();
	double sum_num = 0;
	double sum_denom = 0;
	for (int i=0; i<n; i++)
	{
		sum_num+=(t[i]-x_bar)*(C[i]-y_bar);
		sum_denom+=(t[i]-x_bar)*(t[i]-x_bar);
	}
	vector<double> out(2);
	out[1] = sum_num/sum_denom;
	out[0] = y_bar-out[1]*x_bar;
	return out;

}

double CTimeSeries::mean_t()
{
	double sum = 0;
	for (int i=0; i<n; i++)
		sum += t[i];
	return sum/double(n);

}

int sgn(int val) {
    return (int(0) < val) - (val < int(0));
}

double sgn(double val) {
    return double(double(0) < val) - (val < double(0));
}

CTimeSeries CTimeSeries::add_noise(double std, bool logd)
{
	CTimeSeries X(n);
	for (int i=0; i<n; i++)
	{
		X.t[i] = t[i];
		if (logd==false)
			X.C[i] = C[i]+getnormalrand(0,std);
		else
			X.C[i] = C[i]*exp(getnormalrand(0,std));
	}
	return X;

}

double sum_interpolate(vector<CTimeSeries> BTC, double t)
{
	double sum=0;
	for (int i=0; i<int(BTC.size()); i++)
	{
		sum+=BTC[i].interpol(t);
	}
	return sum;
}


void CTimeSeries::assign_D()
{
	for (int i = 0; i<n; i++)
	{
		double counter = 0;
		for (int j = i + 1; j<n; j++)
		{
			if (C[j] == C[i]) counter += (t[j] - t[j - 1]);
			if (C[j] != C[i]) break;
		}
		D.push_back(counter);
	}
}

void CTimeSeries::clear()
{
	C.clear();
	t.clear();
	D.clear();
	n = 0;
}

double CTimeSeries::wiggle()
{
	if (n>2)
		return 3*(std::fabs(C[n-1])*(t[n-2]-t[n-3])-std::fabs(C[n-2])*(t[n-1]-t[n-3])+std::fabs(C[n-3])*(t[n-1]-t[n-2]))/(t[n-1]-t[n-3])/max(maxfabs(),1e-7);
	else
		return 0;

}

double CTimeSeries::wiggle_corr(int _n)
{
	if (n < _n) return 0;
	double sum=0;
	double var = 0;
	double C_m=0;
	for (int i = 0; i < _n; i++)
	{
		C_m += C[n - i-1] / double(_n);
	}
	for (int i = 0; i < _n-1; i++)
	{
		sum += (C[n - i-1] - C_m)*(C[n - i - 2] - C_m);
	}
	for (int i = 0; i < _n ; i++)
	{
		var += pow(C[n - i-1] - C_m,2);
	}
	if (var == 0)
		return 0;
	else
		return sum / var;
}

bool CTimeSeries::wiggle_sl(double tol)
{
	if (n < 4) return false;
	double mean = std::fabs(C[n - 1] + C[n - 2] + C[n - 3] + C[n - 4]) / 4.0+tol/100;
	double slope1 = (C[n - 1] - C[n - 2]) / (t[n - 1] - t[n - 2])/mean;
	double slope2 = (C[n - 2] - C[n - 3]) / (t[n - 2] - t[n - 3])/mean;
	double slope3 = (C[n - 3] - C[n - 4]) / (t[n - 3] - t[n - 4])/mean;
	if (std::fabs(slope1) < tol && std::fabs(slope2) < tol && std::fabs(slope3) < tol) return false;
	if ((slope1*slope2 < 0) && (slope2*slope3 < 0))
		return true;
	else
		return false;
}

void CTimeSeries::knock_out(double tt)
{
	int i=n-1;
	if (n>0)
	{	while (t[i]>tt)
		{	if (i<n)
			{	t.pop_back();
				C.pop_back();
				n--;
				i--;
			}
		}
        }
}

double CTimeSeries::AutoCor1(int k)
{
	if (k == 0) k = n;
	double sum_product = 0;
	double sum_sq = 0;
	double mean1 = mean();
	for (int i = n - k; i < n - 1; i++)
	{
		sum_product += (C[i] - mean1)*(C[i + 1] - mean1);
		sum_sq += (C[i] - mean1);
	}
	return sum_product / sum_sq;

}
/*
vector<double> CTimeSeries::trend()
{
	double x_bar = mean_t();
	double y_bar = mean();
	double sum_num = 0;
	double sum_denom = 0;
	for (int i = 0; i<n; i++)
	{
		sum_num += (t[i] - x_bar)*(C[i] - y_bar);
		sum_denom += (t[i] - x_bar)*(t[i] - x_bar);
	}
	vector<double> out(2);
	out[1] = sum_num / sum_denom;
	out[0] = y_bar - out[1] * x_bar;
	return out;

}
*/

CTimeSeries CTimeSeries::getcummulative()
{
	CTimeSeries X(n);
	X.t = t;
	X.C[0] = 0;
	for (int i = 1; i<n; i++)
            X.C[i] = X.C[i - 1] + (X.t[i] - X.t[i - 1])*0.5*(C[i] + C[i - 1]);

	return X;
}

bool CTimeSeries::isweighted()
{
    if (weight.size())
        return true;
    else
        return false;
}

CTimeSeries CTimeSeries::getcummulative_direct(int number_of_bins)
{
	CTimeSeries X(number_of_bins+1);
	for (int j=0; j<number_of_bins+1; j++)
        X.t[j] = exp(log(minC()) + j*(log(maxC())-log(minC()))/(number_of_bins));

	for (int i = 0; i<n; i++)
        for (int j=number_of_bins; j>=0; j--)
        {
            if (C[i]<X.t[j])
            {
                if (weight.size())
                {
                    X.C[j] += weight[i];
                    //cout<<"weighted "<<n << "  " << weight[i]<<endl;
                }
                else
                    X.C[j] += 1;
            }
            else
                j=-1;
        }

	X.structured = true;
	cout<<"inside getcum direct: "<<X.C[number_of_bins]<<","<<X.C[0]<<endl;
	return X/X.C[number_of_bins];
}

CTimeSeries CTimeSeries::Exp()
{
	CTimeSeries BTC(n);
	for (int i = 0; i<n; i++)
	{
		BTC.t[i] = t[i];
		BTC.C[i] = exp(C[i]);
	}
	return BTC;
}

CTimeSeries CTimeSeries::fabs()
{
	CTimeSeries BTC = CTimeSeries(n);
	for (int i = 0; i<n; i++)
	{
		BTC.t[i] = t[i];
		BTC.C[i] = std::fabs(C[i]);
	}
	return BTC;
}

CTimeSeries CTimeSeries::derivative()
{
	CTimeSeries out;
	for (int i = 0; i < n - 1; i++)
	{
		out.C.push_back((C[i + 1] - C[i]) / (t[i + 1] - t[i]));
		out.t.push_back((t[i + 1] + t[i])/2);
		out.n++;
	}

	return out;
}

double R2_c(CTimeSeries BTC_p, CTimeSeries BTC_d)
{
	double sumcov = 0;
	double sumvar1 = 0;
	double sumvar2 = 0;
	double sum1 = 0;
	double sum2 = 0;
	double totcount = min(BTC_d.n, BTC_p.n);
	for (int i = 0; i<totcount; i++)
	{
		sumcov += fabs(BTC_d.C[i])*fabs(BTC_p.C[i]) / totcount;
		sumvar1 += BTC_d.C[i] * BTC_d.C[i] / totcount;
		sumvar2 += BTC_p.C[i] * BTC_p.C[i] / totcount;
		sum1 += fabs(BTC_d.C[i]) / totcount;
		sum2 += fabs(BTC_p.C[i]) / totcount;
	}

	return pow(sumcov - sum1*sum2, 2) / (sumvar1 - sum1*sum1) / (sumvar2 - sum2*sum2);
}
double norm2(CTimeSeries BTC1)
{
	double sum = 0;
	for (int i = 0; i<BTC1.n; i++)
		sum += pow(BTC1.C[i], 2);

	return sum;
}
CTimeSeries max(CTimeSeries A, double b)
{
	CTimeSeries S = A;
	for (int i = 0; i<A.n; i++)
		S.C[i] = max(A.C[i], b);
	return S;
}
double pos(double x)
{
	if (x > 0) return x; else return 0;
}
double neg(double x)
{
	return pos(-x);
}
CTimeSeries operator>(CTimeSeries BTC1, CTimeSeries BTC2)
{
	CTimeSeries S = BTC1;
	for (int i = 0; i<min(BTC1.n, BTC2.n); i++)
		S.C[i] = BTC1.C[i] - BTC2.C[i];

	return S;
}

CTimeSeries CTimeSeries::unlog()
{
	CTimeSeries out(n);
	for (int i = 0; i < n; i++)
	{
		out.C[i] = C[i] / exp(t[i]);
		out.t[i] = exp(t[i]);
	}
	return out;

}

CTimeSeries CTimeSeries::distribution_fw(int n_bins, double smoothing_span, int limit, string s)
{
	CTimeSeries out(n_bins + 2);

	CVector C1(C.size() - limit);
	for (int i = 0; i<C1.num; i++)
		C1[i] = C[i + limit];

	double p_start = min(C1);
	double p_end = max(C1)*1.001;
	double dp = abs(p_end - p_start) / n_bins;
	if (dp == 0) return out;
	out.t[0] = p_start - dp / 2;
	out.C[0] = 0;
	for (int i = 0; i<n_bins + 1; i++)
	{
		out.t[i + 1] = out.t[i] + dp;
		out.C[i + 1] = out.C[i];
	}

    if (smoothing_span==0)
    {

        if (weight.size()==0)
        {
            for (int i = 0; i < C1.num; i++)
            {
                if (s=="log")
                    out.C[int((C1[i] - p_start) / dp) + 1] += 1.0 / dp*exp(C1[i]);
                else
                    out.C[int((C1[i] - p_start) / dp) + 1] += 1.0 / dp*C1[i];
            }
        }
        else
        {
            for (int i = 0; i < C1.num; i++)
            {
                if (s=="log")
                    out.C[int((C1[i] - p_start) / dp) + 1] += 1.0 / dp*exp(C1[i])*weight[i];
                else
                    out.C[int((C1[i] - p_start) / dp) + 1] += 1.0 / dp*C1[i]*weight[i];

            }
        }
        if (s == "log")
            return (1/C1.Exp().sum())*out;
        else
            return (1 / C1.sum())*out;
    }
    else
    {
        int span_count = smoothing_span/dp;
        cout<<"span_count"<< span_count<<endl;
        if (weight.size()==0)
            {
                for (int i = 0; i < C1.num; i++)
                {
                    int center = int((C1[i]-p_start)/dp)+1;
                    for (int j=max(0,center-span_count); j<=min(n_bins,center+span_count); j++)
                    {
                        double l_bracket = p_start + (j-1)*dp;
                        double r_bracket = p_start + (j)*dp;
                        //cout << "l:" << l_bracket<< "r:" << r_bracket << endl;
                        double eff_smoothing_span = max(min(min(C[i]-p_start,smoothing_span),p_end-C[i]),dp/10.0);
                        double portion = exp((C1[i]-l_bracket)/eff_smoothing_span)/(1.0+exp((C1[i]-l_bracket)/eff_smoothing_span)) - exp((C1[i]-r_bracket)/eff_smoothing_span)/(1+exp((C1[i]-r_bracket)/eff_smoothing_span));
                        //cout << "portion:" << portion << endl;
                        if (s=="log")
                            out.C[j] += 1.0 / dp*exp(C1[i])*portion;
                        else
                            out.C[j] += 1.0 / dp*C1[i]*portion;
                    }
                }
            }
            else
                for (int i = 0; i < C1.num; i++)
                {
                    int center = int((C1[i]-p_start)/dp)+1;
                    for (int j=max(0,center-3*span_count); j<=min(n_bins,center+3*span_count); j++)
                    {
                        double l_bracket = p_start + (j-1)*dp;
                        double r_bracket = p_start + (j)*dp;
                        double eff_smoothing_span = max(min(min(C[i]-p_start,smoothing_span),p_end-C[i]),dp/10.0);
                        double portion = exp((C1[i]-l_bracket)/eff_smoothing_span)/(1.0+exp((C1[i]-l_bracket)/eff_smoothing_span)) - exp((C1[i]-r_bracket)/eff_smoothing_span)/(1+exp((C1[i]-r_bracket)/eff_smoothing_span));
                        if (s=="log")
                            out.C[j] += 1.0 / dp*weight[i]*exp(C1[i])*portion;
                        else
                            out.C[j] += 1.0 / dp*C1[i]*weight[i]*portion;
                    }
                }

        if (s == "log")
            return (1/C1.Exp().sum())*out;
        else
            return (1 / C1.sum())*out;
    }
}

CTimeSeries CTimeSeries::make_flux_weighted(string log)
{
	CTimeSeries out(n);
	for (int i = 0; i < n; i++)
	{
		if (log=="log")
			out.C[i] = C[i] * exp(t[i]);
		else
			out.C[i] = C[i] * t[i];
		out.t[i] = t[i];
	}

	return (1.0/out.integrate())*out;
}

CTimeSeries CTimeSeries::distribution_log(int n_bins, double smooting_span, int limit)
{
	CTimeSeries out(n_bins + 2);

	CVector C1;
	for (unsigned int i = 0; i<C.size()-limit; i++)
		if (C[i+limit]>0)
			C1.append(log(C[i + limit]));

	double p_start = min(C1);
	double p_end = max(C1)*1.001;
	double dp = abs(p_end - p_start) / n_bins;
	if (dp == 0) return out;
	out.t[0] = p_start - dp / 2;
	out.C[0] = 0;
	for (int i = 0; i<n_bins + 1; i++)
	{
		out.t[i + 1] = out.t[i]* exp(dp);
		out.C[i + 1] = out.C[i];
	}

	for (int i = 0; i < C1.num; i++)
	{
		int k = int((C1[i] - p_start) / dp)+1;
		out.C[k] += 1.0 / C1.num / (out.t[k + 1] - out.t[k]);

	}

	return out;
}

CTimeSeries CTimeSeries::standardize()
{
    double _std = std();
    double _mean = mean();
    CTimeSeries out = (*this - _mean)/_std;
    return out;
}

double CTimeSeries::autocorrelation()
{
    double sum2=0;
    double sumprod=0;
    for (int i=1; i<n-1; i++)
    {
        sum2+=0.5*(pow(C[i],2)+pow(C[i+1],2));
        sumprod += C[i]*C[i+1];
    }
    return sumprod/sum2;
}

double TS::correlation(const CTimeSeries &TS1, const CTimeSeries &TS2)
{
    double sum2=0;
    double sumprod=0;
    for (int i=0; i<min(TS1.n,TS2.n); i++)
    {
        sum2+=0.5*(pow(TS1.C[i],2)+pow(TS2.C[i],2));
        sumprod += TS1.C[i]*TS2.C[i];
    }
    return sumprod/sum2;
}

CBTC CBTC::normalize_by_max()
{
    double _max = maxC();
    CBTC out;
    for (int i=0; i<n; i++)
    {
        out.append(t[i], C[i]/_max);
    }
    return out;
}

