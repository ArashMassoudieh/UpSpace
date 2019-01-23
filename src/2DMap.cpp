#include "2DMap.h"
#include <fstream>
#include <string>
#include "StringOP.h"

TDMap::TDMap()
{
    //ctor
}

TDMap::~TDMap()
{
    //dtor
}

TDMap::TDMap(unsigned int number_of_bins, double low_lim, double up_lim)
{
    reset(number_of_bins, number_of_bins, low_lim, up_lim, low_lim, up_lim);
}

TDMap::TDMap(unsigned int number_of_bins_x, unsigned int number_of_bins_y, double low_lim_x, double up_lim_x, double low_lim_y, double up_lim_y)
{
    reset(number_of_bins_x, number_of_bins_y, low_lim_x, up_lim_x, low_lim_y, up_lim_y);
}

TDMap::TDMap(const TDMap& other)
{
    val = other.val;
    x_bin = other.x_bin;
    y_bin = other.y_bin;
}

TDMap& TDMap::operator=(const TDMap& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    val = rhs.val;
    x_bin = rhs.x_bin;
    y_bin = rhs.y_bin;
    return *this;
}

void TDMap::reset(unsigned int number_of_bins_x, unsigned int number_of_bins_y, double _low_lim_x, double _up_lim_x, double _low_lim_y, double _up_lim_y)
{
    low_lim_x = _low_lim_x;
    up_lim_x = _up_lim_x;
    low_lim_y = _low_lim_y;
    up_lim_y = _up_lim_y;
    val.resize(number_of_bins_x);
    for (unsigned int i=0; i<number_of_bins_x; i++)
        val[i].resize(number_of_bins_y);

    x_bin.resize(number_of_bins_x+1);
    y_bin.resize(number_of_bins_y+1);

    for (unsigned int i=0; i<number_of_bins_x+1; i++)
        x_bin[i] = low_lim_x + (up_lim_x-low_lim_x)/number_of_bins_x*i;


    for (unsigned int i=0; i<number_of_bins_y+1; i++)
        y_bin[i] = low_lim_y + (up_lim_y-low_lim_y)/number_of_bins_y*i;
}

void TDMap::set_val(unsigned int i, unsigned int j, double value)
{
    val[i][j] = value;
}

void TDMap::add_val(unsigned int i, unsigned int j, double value)
{
    val[i][j] += value;
}

void TDMap::add_val(double x, double y, double value)
{
    int i = int((x-low_lim_x)/(up_lim_x-low_lim_x))*val.size();
    int j = int((y-low_lim_y)/(up_lim_y-low_lim_y))*val.size();
    if (i>=0 && j>=0 && i<val.size() && j<val[0].size())
        val[i][j] += value;
}

double TDMap::sum()
{
    double sum = 0;
    for (unsigned int i=0; i<val.size(); i++)
        for (unsigned int j=0; j<val[i].size(); j++)
            sum+=val[i][j];

    return sum;
}

double TDMap::marginal_x(unsigned int i)
{
    double sum = 0;
    for (unsigned int j=0; j<val[i].size(); j++)
        sum+=val[i][j];

    return sum;
}

double TDMap::marginal_y(unsigned int j)
{
    if (j>val[0].size())
        return 0;
    double sum = 0;
    for (unsigned int i=0; i<val.size(); i++)
        sum+=val[i][j];

    return sum;
}

vector<double> TDMap::marginal_x()
{
    vector<double> out;
    for (unsigned int i=0; i<val.size(); i++)
        out.push_back(marginal_x(i));

    return out;
}

vector<double> TDMap::marginal_y()
{
    vector<double> out;
    for (unsigned int j=0; j<val[0].size(); j++)
        out.push_back(marginal_y(j));

    return out;
}


void TDMap::normalize()
{
    double s = sum();
    for (unsigned int i=0; i<val.size(); i++)
        for (unsigned int j=0; j<val[i].size(); j++)
            val[i][j] = val[i][j]/s*val.size()*val[0].size()/(up_lim_x-low_lim_x)/(up_lim_y-low_lim_y);

}

void TDMap::normalise_x()
{
    vector<double> s = marginal_x();
    for (unsigned int i=0; i<val.size(); i++)
        for (unsigned int j=0; val[i].size(); j++)
            val[i][j] = val[i][j]/s[i];

}

void TDMap::normalize_y()
{
    vector<double> s = marginal_y();
    for (unsigned int i=0; i<val.size(); i++)
        for (unsigned int j=0; val[i].size(); j++)
            val[i][j] = val[i][j]/s[j];
}

double TDMap::get_val(unsigned int i, unsigned int j)
{
    if (i<val.size() && j<val[0].size())
        return val[i][j];
    return 0;
}

void TDMap::writetofile(string filename)
{
    ofstream file(filename);

    for (unsigned int i=0; i<val.size(); i++)
    {
        file << "," << (x_bin[i] + x_bin[i+1])/2;
    }
    file << endl;
    for (unsigned int j=0; j<val[0].size(); j++)
    {
        file << (y_bin[j] + y_bin[j+1])/2 << ",";
        for (unsigned int i=0; i<val.size(); i++)
            file << val[i][j] << ",";
        file << endl;
    }
    file.close();
}

void TDMap::writetofile_GNU(string filename, string pngfilename, string xlabel, string ylabel, string title, bool logscale)
{
    if (pngfilename=="")
        pngfilename = split(filename,'.')[0] + ".png";
    ofstream file(filename);
    if (logscale)
    {
        file << "set xrange ["<<low_lim_x+0.49*(up_lim_x-low_lim_x)/val.size()<<":"<<up_lim_x-0.49*(up_lim_x-low_lim_x)/val.size()<<"]"<<endl;
        file << "set yrange ["<<low_lim_y+0.49*(up_lim_y-low_lim_y)/val[0].size()<<":"<<up_lim_y-0.49*(up_lim_y-low_lim_y)/val[0].size()<<"]"<<endl;
    }
    else
    {
        file << "set xrange ["<<low_lim_x+0*0.5*(up_lim_x-low_lim_x)/val.size()<<":"<<up_lim_x-0*0.5*(up_lim_x-low_lim_x)/val.size()<<"]"<<endl;
        file << "set yrange ["<<low_lim_y+0*0.5*(up_lim_y-low_lim_y)/val[0].size()<<":"<<up_lim_y-0*0.5*(up_lim_y-low_lim_y)/val[0].size()<<"]"<<endl;
    }

    file << "show xrange" << endl;
    file << "show yrange" << endl;
    file << "set xlabel \"" << xlabel << "\"" << endl;
    file << "show xlabel" << endl;
    file << "set ylabel \"" << ylabel << "\"" << endl;
    file << "show ylabel" << endl;
    file << "set title \"" << title << "\"" << endl;
    file << "show title" << endl;
    file << "set palette rgb -21,-22,-23" << endl;
    file << "set logscale cb" << endl;
    if (logscale) file << "set logscale xy" << endl;
    file << "set format cb \"10^{%T}\"" << endl;
    file << "set pm3d map interpolate 5,5"<<endl;
    file << "set size square" << endl;
    file << "splot \"-\" matrix using ($1*"<<(up_lim_x-low_lim_x)/val.size()<<"+"<<low_lim_x+0.5*(up_lim_x-low_lim_x)/val.size()<<"):($2*"<<(up_lim_y-low_lim_y)/val.size()<<"+"<<low_lim_x+0.5*(up_lim_y-low_lim_y)/val.size()<<"):3 notitle"<<endl;
    for (unsigned int j=0; j<val[0].size(); j++)
    {
        for (unsigned int i=0; i<val.size(); i++)
            file << val[i][j] << " ";
        file << endl;
    }
    file<<"e"<<endl;
    file<<"e"<<endl;
    file<<"set term png"<<endl;
    file<<"set output \"" << pngfilename << "\""<<endl;
    file<<"replot"<<endl;
    file<<"set term x11"<<endl;
    file.close();
}
