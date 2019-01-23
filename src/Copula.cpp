#include "Copula.h"
#include "NormalDist.h"


Copula::Copula()
{
}


Copula::~Copula()
{
}

double Copula::c(vector <double> u, double k, CMatrix &distanceMatrix)
{
	double c;
	//if (copulaType == gaussian)
	//{
		CMatrix R(u.size());
		for (int i = 0; i < u.size(); i++)
			for (int j = 0; j < u.size(); j++)
				R[i][j] = exp(-k*distanceMatrix[i][j]);
		CVector uu = CVector(u);
                CMatrix mm = inv(R) - Identity(u.size());
                c = 1 / sqrt(abs(R.det())) * exp(-0.5*dotproduct(uu, mm*uu));
	//}
	return c;
}
double Copula::logc(vector <double> u, double k, CMatrix &distanceMatrix)
{
	double c;
	//if (copulaType == gaussian)
	//{
		CMatrix R(int(u.size()));
		for (int i = 0; i < int(u.size()); i++)
			for (int j = 0; j < int(u.size()); j++)
				R[i][j] = exp(-k*(distanceMatrix[i][j]));
		double a = det(R);
		CMatrix b = inv(R);
		if (!b.getnumrows())
			return -30000;

		CMatrix mm = b - Identity(u.size());
                CVector uu = CVector(u);
                CVector e = mm*uu;
		double d = dotproduct(CVector(u), e);
		c = -0.5 * log(a) - 0.5 * d;
	//	c = - 0.5 * log(abs(R.det())) -0.5*dotproduct(CVector(u), ((inv(R) - Identity(u.size()))*CVector(u)));
	//}
	return c;
}


