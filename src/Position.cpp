#include "Position.h"



CPosition::CPosition()
{
	v = CVector(2);
	y = 0;
	x = 0;
	t = 0;
	u = 0;
	z = 0;
        weight = 1;
}

CPosition::CPosition(int n_coor)
{
	v = CVector(n_coor);
	y = 0;
	x = 0;
	t = 0;
	u = 0;
	z = 0;
    weight = 1;
}

CPosition::CPosition(const CPosition & P)
{
	t = P.t;
	x = P.x;
	y = P.y;
	u = P.u;
	z = P.z;
	v = P.v;
        weight = P.weight;
}


CPosition& CPosition::operator=(const CPosition & P)
{
	t = P.t;
	x = P.x;
	y = P.y;
	u = P.u;
	z = P.z;
	v = P.v;
        weight = P.weight;
	return *this;
}

CPosition::~CPosition()
{
}

double CPosition::getvar(string var)
{
	if (var == "t") return t;
	if (var == "u") return u;
	if (var == "z") return z;
	if (var == "x") return x;
	if (var == "y") return y;
	if (var == "vx") return v[0];
	if (var == "vy") return v[1];
	if (var == "v_eff") {if (t > 0) return x / t; else return 0;};
	if (var == "t_eff") {if (x > 0) return t / x; else return 0;};
        if (var == "weight") return weight;
	return 0;
}

CPosition operator+(const CPosition p1, const CPosition p2)
{
	CPosition p;
	p.v = p1.v + p2.v;
	p.y = p1.y + p2.y;
	p.u = p1.u + p2.u;
	p.z = p1.z + p2.z;
	p.x = p1.x + p2.x;
	p.t = p1.t + p2.t;
        p.weight = (p1.weight + p2.weight)/2.0;
	return p;
}

CPosition operator-(const CPosition p1, const CPosition p2)
{
    CPosition p;
    p.v = p1.v - p2.v;
    p.y = p1.y - p2.y;
    p.u = p1.u - p2.u;
    p.z = p1.z - p2.z;
    p.x = p1.x - p2.x;
    p.t = p1.t - p2.t;
    p.weight = (p1.weight + p2.weight)/2.0;
    return p;
}

CPosition operator*(double d, CPosition p1)
{
	CPosition p;
	p.v = d*p1.v;
	p.y = d*p1.y;
	p.u = d*p1.u;
	p.z = d*p1.z;
	p.x = d*p1.x;
	p.t = d*p1.t;
        p.weight = p1.weight;
	return p;
}

CPosition operator*(CPosition p1, double d)
{
	return d*p1;
}

CPosition operator/(CPosition p1, double d)
{
	return p1*(1 / d);
}


