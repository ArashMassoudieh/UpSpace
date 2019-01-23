#include "Grid.h"
#include "StringOP.h"
#include "NormalDist.h"
#ifdef QT_version
#include <qtextbrowser.h>
#include "heteval.h"
#include "qdebug.h"
#include <QCursor>
#endif // Qt_version


#include <sys/resource.h>

//#define Qt_version

CGrid::CGrid()
{
}

CGrid::~CGrid()
{
}

vector<ijval> CGrid::get_closest_K_dets(int i, int j, int n)
{
	vector<ijval> out;
	double max_dist = 0;
	for (int k = 1; k < max(GP.nx, GP.ny); k++)
	{
		double min_dist = 1e12;
		int jj = j - k;
		if (jj>=0)
			for (int ii = max(i - k,0); ii <= min(i + k,GP.nx-1); ii++)
			{
				double dist2 = ((i - ii)*GP.dx*(i - ii)*GP.dx + (j - jj)*GP.dy*(j - jj)*GP.dy);
				min_dist = min(min_dist, dist2);
				if (p[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
						ijval pp;
						pp.i = ii;
						pp.j = jj;
						pp.val = sqrt(dist2);
						out.push_back(pp);
						max_dist = max(dist2, max_dist);
					}
				}

			}
		jj = j + k;
		if (jj < GP.ny)
			for (int ii = max(i - k, 0); ii <= min(i + k, GP.nx - 1); ii++)
			{
				double dist2 = ((i - ii)*GP.dx*(i - ii)*GP.dx + (j - jj)*GP.dy*(j - jj)*GP.dy);
				min_dist = min(min_dist, dist2);
				if (p[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
						ijval pp;
						pp.i = ii;
						pp.j = jj;
						pp.val = sqrt(dist2);
						out.push_back(pp);
						max_dist = max(dist2, max_dist);
					}
				}

			}
		int ii = i - k;
		if (ii >= 0)
			for (int jj = max(j - k+1, 0); jj <= min(j + k-1, GP.ny - 1); jj++)
			{
				double dist2 = ((i - ii)*GP.dx*(i - ii)*GP.dx + (j - jj)*GP.dy*(j - jj)*GP.dy);
				min_dist = min(min_dist, dist2);
				if (p[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
						ijval pp;
						pp.i = ii;
						pp.j = jj;
						pp.val = sqrt(dist2);
						out.push_back(pp);
						max_dist = max(dist2, max_dist);
					}
				}

			}
		ii = i + k;
		if (ii < GP.nx)
			for (int jj = max(j - k + 1, 0); jj <= min(j + k - 1, GP.ny - 1); jj++)
			{
				double dist2 = ((i - ii)*GP.dx*(i - ii)*GP.dx + (j - jj)*GP.dy*(j - jj)*GP.dy);
				min_dist = min(min_dist, dist2);
				if (p[ii][jj].k_det)
				{
					if ((dist2 < max_dist) || (int(out.size()) < n))
					{
                                            ijval pp;
                                            pp.i = ii;
                                            pp.j = jj;
                                            pp.val = sqrt(dist2);
                                            out.push_back(pp);
                                            max_dist = max(dist2, max_dist);
					}
				}

			}
		if ((int(out.size()) >= n) && min_dist > max_dist)
		{
                    k = max(GP.nx, GP.ny);
                    ijval pp;
                    pp.i = i;
                    pp.j = j;
                    pp.val = 0;
                    out.push_back(pp);
                    return out;
		}
	}
	ijval pp;
	pp.i = i;
	pp.j = j;
	pp.val = 0;
	out.push_back(pp);
	return out;
}

correl_mat_vec CGrid::get_correll_matrix_vec(int i, int j)
{
	correl_mat_vec M;
	vector<ijval> ij = get_top_n(get_closest_K_dets(i, j, min(n_k_dets,field_gen.max_correl_n)+1), min(n_k_dets, field_gen.max_correl_n)+1);
#ifdef  arma

	M.M_22 = CMatrix_arma(ij.size() - 1);
	M.V_21 = CVector_arma(ij.size() - 1);
	M.V_RHS = CVector_arma(ij.size() - 1);
#else
	M.M_22 = CMatrix(ij.size() - 1);
	M.V_21 = CVector(ij.size() - 1);
	M.V_RHS = CVector(ij.size() - 1);
#endif //  arma
	for (int ii = 1; ii < int(ij.size()); ii++)
	{
		M.V_21[ii-1] = exp(-sqrt((i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x*(i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x + (j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y*(j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y) );
		M.V_RHS[ii - 1] = p[ij[ii].i][ij[ii].j].K_gauss[0];
		for (int jj = 1; jj < int(ij.size()); jj++)
		{
#ifdef arma
		M.M_22(ii-1,jj-1) = exp(-sqrt((ij[jj].i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x*(ij[jj].i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x + (ij[jj].j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y*(ij[jj].j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y));
#else
		M.M_22[ii - 1][jj - 1] = exp(-sqrt((ij[jj].i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x*(ij[jj].i - ij[ii].i)*GP.dx/ field_gen.k_correlation_lenght_scale_x + (ij[jj].j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y*(ij[jj].j - ij[ii].j)*GP.dy/ field_gen.k_correlation_lenght_scale_y));
#endif // arma
		}
	}
	return M;
}

void CGrid::assign_K_gauss(int i, int j)
{
	CNormalDist ND;
	correl_mat_vec M = get_correll_matrix_vec(i, j);
	double mu;
	double sigma;
	if (M.V_RHS.num == 0)
	{
		mu = 0;
		sigma = 1;
	}
	else
	{
            CMatrix_arma M_inv = inv(M.M_22);
            mu = dotproduct(M_inv*M.V_21, M.V_RHS);
            sigma = 1.0 - dotproduct(M_inv*M.V_21, M.V_21);
	}

	double K_gauss = getnormalrand(mu, sigma);
	p[i][j].k_det = true;
	n_k_dets++;
	p[i][j].K_gauss[0] = K_gauss;
	p[i][j].K[0] = map_to_KCDF(getnormalcdf(K_gauss));
}

void CGrid::assign_K_gauss()
{
	int n_filled = 0;
	srand(time(NULL));
	while (n_filled<GP.nx*GP.ny)
	{
        if (n_filled%100==0)
		set_progress_value(double(n_filled) / double(GP.nx * GP.ny));
            int i = unitrandom()*(GP.nx-1) + 0.5;
            int j = unitrandom()*(GP.ny-1) + 0.5;
            if (!p[i][j].k_det)
            {
                assign_K_gauss(i, j);
                n_filled++;
            }
    }
    cout<<endl;

}

bool CGrid::getfromfile(string filename, int _nx, int _ny)
{

        bool file_not_found;
	GP.nx = _nx;
	GP.ny = _ny;
	p.resize(_nx);
	for (int i = 0; i < _nx; i++)
	{
		p[i].resize(_ny);
		for (int j = 0; j < GP.ny; j++)
			p[i][j].V = CVector(2);
	}

	ifstream file(filename);
	vector<string> s;
	if (file.good() == false)
	{
		file_not_found = true;
		return false;
	}

	double x_0, y_0;
	if (file.good())
		for (int j = 0; j < GP.ny; j++)
			for (int i = 0; i < GP.nx; i++)
			{

				s = getline(file);
				if ((i == 0) && (j == 0))
				{
					x_0 = atof(s[0].c_str());
					y_0 = atof(s[1].c_str());
				}
				if (i == 1)
					GP.dx = atof(s[0].c_str()) - x_0;

				if (j == 1)
					GP.dy = atof(s[1].c_str()) - y_0;

				p[i][j].V[0] = atof(s[2].c_str());
				p[i][j].V[1] = atof(s[3].c_str());
			}
	file.close();
	return true;

}

bool CGrid::getKfromfile(string filename, int _nx, int _ny, int _nx_act, int _ny_act)
{
	bool file_not_found;
	GP.nx = _nx;
	GP.ny = _ny;
	p.resize(_nx);
	for (int i = 0; i < _nx; i++)
	{
		p[i].resize(_ny);
		for (int j = 0; j < GP.ny; j++)
		{
			p[i][j].V = CVector(2);
			p[i][j].K = CVector(2);
		}
	}

	ifstream file(filename);
	vector<string> s;
	if (file.good() == false)
	{
		file_not_found = true;
		return false;
	}

	double x_0, y_0;
	if (file.good())
		for (int j = 0; j < _ny_act; j++)
			for (int i = 0; i < _nx_act; i++)
			{

				s = getline(file);
				if ((i == 0) && (j == 0))
				{
					x_0 = atof(s[0].c_str());
					y_0 = atof(s[1].c_str());
				}
				if (i == 1)
					GP.dx = atof(s[0].c_str()) - x_0;

				if (j == 1)
					GP.dy = atof(s[1].c_str()) - y_0;

				//p[i][j].K[0] = atof(s[2].c_str());
				if ((i<GP.nx) && (j<GP.ny)) p[i][j].K[0] = atof(s[3].c_str());
			}
	file.close();
	return true;

}

bool CGrid::createuniform(double K, int _nx, int _ny)
{
	GP.nx = _nx;
	GP.ny = _ny;
	p.resize(_nx);
	for (int i = 0; i < _nx; i++)
	{
		p[i].resize(_ny);
		for (int j = 0; j < GP.ny; j++)
		{
			p[i][j].V = CVector(2);
			p[i][j].K = CVector(2);
			p[i][j].K[0] = K;
			p[i][j].K[1] = K;
		}

	}

	return false;
}




CGrid::CGrid(string filename)
{
	ifstream file(filename);
	if (!file.good())
        {
#if QT_version
            show_in_window("File " + filename + "was not found!");
#endif
            cout << "File " + filename + "was not found!" << endl;
            return;
        }
        else
        {
#if QT_version
            show_in_window("Reading " + filename);
#endif
            cout << "Reading " + filename + "..." << endl;

        }
	vector<string> s;

	int n_for = 1;
	int last_for_location = 0;
	while (!file.eof())
	{
		s = getline(file);
		if (s.size())
		{

			if (tolower(s[0]) == "pathout")
            {
                pathout = s[1].c_str();
                pathout.erase( std::remove(pathout.begin(), pathout.end(), '\r'), pathout.end() );
                cout<<"Output path is: " + pathout << endl;
            }
			if (tolower(s[0]) == "for")
			{
				cout << "For loop encountered"<<endl;
				n_for = atoi(s[1].c_str());
				last_for_location = commands.size();
			}
			if (tolower(s[0]) == "endfor")
			{
                _command command;
                command.command = "clear_all";
                commands.push_back(command);
                int end_for_location = commands.size();
                for (int i = 0; i < n_for-1; i++)
                {
                    _command write_commend;
                    write_commend.command = "write";
                    write_commend.parameters["content"] = "starting realization " + numbertostring(i + 1);
                    commands.push_back(write_commend);
                    for (int j = last_for_location; j < end_for_location; j++)
                    {
                        _command new_command = commands[j];
                        if (new_command.parameters["filename"].size() > 0)
                                new_command.parameters["filename"] = insert_counter_in_file_name(new_command.parameters["filename"], i+1);
                        if (new_command.parameters["filename_x"].size() > 0)
                                new_command.parameters["filename_x"] = insert_counter_in_file_name(new_command.parameters["filename_x"], i+1);
                        if (new_command.parameters["filename_y"].size() > 0)
                                new_command.parameters["filename_y"] = insert_counter_in_file_name(new_command.parameters["filename_y"], i+1);
                        if (new_command.parameters["filename_mag"].size() > 0)
                                new_command.parameters["filename_mag"] = insert_counter_in_file_name(new_command.parameters["filename_mag"], i+1);
                        if (new_command.parameters["dist_filename"].size() > 0)
                                new_command.parameters["dist_filename"] = insert_counter_in_file_name(new_command.parameters["dist_filename"], i+1);
                        if (new_command.parameters["ranks_filename"].size() > 0)
                                new_command.parameters["ranks_filename"] = insert_counter_in_file_name(new_command.parameters["ranks_filename"], i+1);
                        if (new_command.parameters["ranks_filename"].size() > 0)
                                new_command.parameters["ranks_filename"] = insert_counter_in_file_name(new_command.parameters["normal_filename"], i+1);
                        if (new_command.parameters["OU_parameters_filename"].size() > 0)
                                new_command.parameters["OU_parameters_filename"] = insert_counter_in_file_name(new_command.parameters["OU_parameters_filename"], i+1);

                        commands.push_back(new_command);

                    }
                    _command command;
                    command.command = "clear_all";
                    commands.push_back(command);

                }
                n_for = 1;

            }
            if (tolower(s[0]) == "k_dist")
            {

                marginal_K_dist_type = s[1];
                for (int j = 2; j < s.size(); j++)
                    marginal_K_dist_params.push_back(atof(s[j].c_str()));
                cout << "Hydraulic conductivity marginal distribution initialized ..." << endl;
            }
            if (tolower(s[0]) == "command")
            {
                _command command;
                command.command = tolower(s[1]);
                for (int j = 2; j < s.size(); j++)
                {   if (split(s[j], '=').size()>1)
                        command.parameters[split(s[j], '=')[0]] = split(s[j], '=')[1];

                }
                cout << "command: " << s[1] << endl;
                commands.push_back(command);
            }
            if (tolower(s[0]) == "x0_trajs") trajectory_params.x0_trajs = atof(s[1].c_str());
		}
	}


}

void CGrid::creategrid(int _nx, int _ny, double _dx, double _dy)
{
	GP.nx = _nx;
	GP.ny = _ny;
	GP.dx = _dx;
	GP.dy = _dy;
	p.resize(GP.nx);
	for (int i = 0; i < GP.nx; i++)
	{
		p[i].resize(GP.ny);
		for (int j = 0; j < GP.ny; j++)
		{
			p[i][j].V = CVector(2);
			p[i][j].K = CVector(2);
			p[i][j].K_gauss = CVector(2);
		}
	}

}

void CGrid::writeasmatrix(string filename, int component)
{
	ofstream file(filename);
	for (int i = 0; i < GP.nx; i++)
	{
		for (int j = 0; j < GP.ny; j++)
			file << p[i][j].V[component]<<",";
		file << endl;
	}


}

void CGrid::writeasmatrixK(string filename, int component)
{
	ofstream file(filename);
	for (int i = 0; i < GP.nx; i++)
	{
		for (int j = 0; j < GP.ny; j++)
			file << p[i][j].K[component] << ",";
		file << endl;
	}
}

CVector CGrid::getvelocity(point pp)
{
	int i_floar = int(pp.x / GP.dx);
	int j_floar = int(pp.y / GP.dy);
	if ((i_floar > GP.nx - 2) || (j_floar>GP.ny - 2) || (i_floar<0) || (j_floar<0))
        {   //qDebug() << "Empty CVector returned";
            return CVector();

        }
	CVector V1 = p[i_floar][j_floar].V + (1.0 / GP.dx*(pp.x - GP.dx*i_floar))*(p[i_floar + 1][j_floar].V - p[i_floar][j_floar].V);
	CVector V2 = p[i_floar][j_floar+1].V + (1.0 / GP.dx*(pp.x - GP.dx*i_floar))*(p[i_floar + 1][j_floar+1].V - p[i_floar][j_floar+1].V);
	CVector V = V1 + (1.0 / GP.dy*(pp.y - GP.dy*j_floar))*(V2 - V1);
	return V;
}

double CGrid::interpolate_K(double x, double y)
{
	point pp;
	pp.x = x;
	pp.y = y;
	int i_floar = int(pp.x / GP.dx);
	int j_floar = int(pp.y / GP.dy);
	if ((i_floar > GP.nx - 2) || (j_floar>GP.ny - 2) || (i_floar<0) || (j_floar<0)) return 0;
	CVector K1 = p[i_floar][j_floar].K + (1.0 / GP.dx*(pp.x - GP.dx*i_floar))*(p[i_floar + 1][j_floar].K - p[i_floar][j_floar].K);
	CVector K2 = p[i_floar][j_floar + 1].K + (1.0 / GP.dx*(pp.x - GP.dx*i_floar))*(p[i_floar + 1][j_floar + 1].K - p[i_floar][j_floar + 1].K);
	CVector K = K1 + (1.0 / GP.dx*(pp.y - GP.dy*j_floar))*(K2 - K1);
	return K[0];
}

vector<double> CGrid::interpolate_V(double x, double y)
{
	point pp;
	pp.x = x;
	pp.y = y;
	return getvelocity(pp).vec;
}

CPathway CGrid::gettrajectory(point pp, double dt, double t_end)
{
	point pt = pp;
	CPathway Trajectory;
        Trajectory.weight = getvelocity(pt)[0];
	double t = 0;

	while (t < t_end)
	{
            CVector V = getvelocity(pt);
            if (V.getsize() == 0) return Trajectory;
            CPosition ps;
            ps.x = pt.x;
            ps.y = pt.y;
            ps.t = t;
            ps.v = V;
            Trajectory.append(ps);
            pt.x += V[0] * dt;
            pt.y += V[1] * dt;
            t += dt;
	}

	return Trajectory;
}

CPathway CGrid::gettrajectory_fix_dx(point pp, double dx, double x_end)
{
    point pt = pp;
    CPathway Trajectory;
    if (weighted)
        Trajectory.weight = getvelocity(pt)[0];
    else
        Trajectory.weight = 1;
    double t = 0;

    bool ex = false;
    while (pt.x < x_end && ex==false)
    {
        CVector V = getvelocity(pt);
        if (V.num == 2)
        {
            double dt = fabs(dx / sqrt(pow(V[0],2)+pow(V[1],2)));
            if (V.getsize() == 0) return Trajectory;
            CPosition ps;
            CVector V = getvelocity(pt);
            ps.x = pt.x;
            ps.y = pt.y;
            ps.t = t;
            ps.v = V;
            Trajectory.append(ps);
            pt.x += V[0] * dt;
            pt.y += V[1] * dt;
            t += dt;
        }
        else ex = true;
    }

    return Trajectory;
}


CPathwaySet CGrid::gettrajectories(double dt, double t_end)
{
    CPathwaySet X;
    for (int i = 0; i < int(pts.size()); i++)
    {
        cout << i << endl;
        CPathway X1 = gettrajectory(pts[i], dt, t_end);
        if (weighted)
        {   X.weighted = true;
            X.append(X1);
        }
        else
        {   X.weighted = false;
            X.append(X1);
        }
        set_progress_value(double(i) / double(pts.size()));
    }
    cout<<endl;
    return X;
}

CPathwaySet CGrid::gettrajectories_fixed_dx(double dx, double x_end)
{
    //qDebug() << "Simulating trajectories"<<endl;

    CPathwaySet X;
    for (int i = 0; i < int(pts.size()); i++)
    {
//	qDebug() << i << endl;

        CPathway X1 = gettrajectory_fix_dx(pts[i], dx, x_end);
        //cout << "\r" << "Trajectory #"<<  i << " Weight: " << X1.weight<< std::flush;
        if (weighted)
            {   X.weighted = true;
                X.append(X1);
            }
            else
            {   X.weighted = false;
                X.append(X1);
            }
        set_progress_value(double(i) / double(pts.size()));
    }
    cout<<endl;
    return X;
}


CBTC CGrid::initialize(int numpoints, double x_0, double smoothing_factor, bool _weighted)
{
    weighted = _weighted;
    cout << "weighted = " << weighted<<endl;
    CBTC vels;
    if (!_weighted)
    {
        cout << "weighted = " << weighted<<endl;
        int burnout = 0;
        double v_max = get_v_btc(x_0,0).maxC();
        double y_0 = unitrandom()*GP.dy*(GP.ny-1);
        point pt_0; pt_0.x = x_0; pt_0.y = y_0;
        double v_x = getvelocity(pt_0)[0];
        pts.push_back(pt_0);
        for (int i = 1; i < numpoints; i++)
        {
            bool accepted = false;
            while (!accepted)
            {
                y_0 = unitrandom()*GP.dy*(GP.ny-1);
                pt_0.x = x_0; pt_0.y = y_0;
                v_x = getvelocity(pt_0)[0];
                double u = unitrandom();
                if (u < (v_x / v_max/5)) accepted = true;
            }
            pts.push_back(pt_0);
            vels.append(i,v_x);

            set_progress_value(double(i) / double(numpoints));
        }

    }
    else
    {   double y_0 = unitrandom()*GP.dy*(GP.ny-1);
        point pt_0; pt_0.x = x_0; pt_0.y = y_0;
        double v_x = getvelocity(pt_0)[0];
        pts.push_back(pt_0);
        for (int i = 1; i < numpoints; i++)
        {
            y_0 = unitrandom()*GP.dy*(GP.ny-1);
            pt_0.x = x_0; pt_0.y = y_0;
            double v_xp = v_x;
            v_x = getvelocity(pt_0)[0];
            double u = unitrandom();
            pts.push_back(pt_0);
            vels.append(i,v_x,v_x);
            set_progress_value(double(i) / double(numpoints));
        }
        cout<<endl;
    }
    return vels.distribution(40,(vels.maxC()-vels.minC())*smoothing_factor);
}

CMatrix_arma_sp CGrid::create_stiffness_matrix_arma()
{
	string averaging = "arithmetic";
//	qDebug() << "Creating stiffness matrix" << endl;
        CMatrix_arma_sp K((GP.nx+1)*(GP.ny+1), (GP.nx + 1)*(GP.ny + 1));
//	qDebug() << "Stiffness matrix created" << endl;
	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			K.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = avg(p[i-1][j-1].K[0], p[i-1][j].K[0], averaging) / (GP.dx*GP.dx);
			K.matr(get_cell_no(i, j), get_cell_no(i, j)) = -(avg(p[i-1][j-1].K[0], p[i-1][j].K[0], averaging) + avg(p[i][j-1].K[0], p[i][j].K[0], averaging)) / (GP.dx*GP.dx);
			K.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) = avg(p[i][j-1].K[0], p[i][j].K[0], averaging) / (GP.dx*GP.dx);

			K.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = avg(p[i-1][j-1].K[0], p[i][j-1].K[0], averaging) / (GP.dy*GP.dy);
			K.matr(get_cell_no(i, j), get_cell_no(i, j)) += -(avg(p[i-1][j-1].K[0], p[i][j-1].K[0], averaging) + avg(p[i-1][j].K[0], p[i][j].K[0], averaging)) / (GP.dy*GP.dy);
			K.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = avg(p[i-1][j].K[0], p[i][j].K[0], averaging) / (GP.dy*GP.dy);
		}
		// top boundary
		int j = GP.ny;
		K.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = 1;
		K.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;

		// bottom boundary
		j = 0;
		K.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = 1;
		K.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;
	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		K.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		K.matr(get_cell_no(i, j), get_cell_no(i+1, j)) = 1;

	}

	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		K.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		K.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = 1;
		}

	return K;

}

CMatrix CGrid::create_stiffness_matrix()
{
	string averaging = "arithmetic";
	CMatrix K((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	CVector V((GP.nx + 1)*(GP.ny + 1));
	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			K[get_cell_no(i, j)][get_cell_no(i - 1, j)] = avg(p[i - 1][j - 1].K[0], p[i - 1][j].K[0], averaging) / (GP.dx*GP.dx);
			K[get_cell_no(i, j)][get_cell_no(i, j)] = -(avg(p[i - 1][j - 1].K[0], p[i - 1][j].K[0], averaging) + avg(p[i][j - 1].K[0], p[i][j].K[0], averaging)) / (GP.dx*GP.dx);
			K[get_cell_no(i, j)][get_cell_no(i + 1, j)] = avg(p[i][j - 1].K[0], p[i][j].K[0], averaging) / (GP.dx*GP.dx);

			K[get_cell_no(i, j)][get_cell_no(i, j - 1)] = avg(p[i - 1][j - 1].K[0], p[i][j - 1].K[0], averaging) / (GP.dy*GP.dy);
			K[get_cell_no(i, j)][get_cell_no(i, j)] += -(avg(p[i - 1][j - 1].K[0], p[i][j - 1].K[0], averaging) + avg(p[i - 1][j].K[0], p[i][j].K[0], averaging)) / (GP.dy*GP.dy);
			K[get_cell_no(i, j)][get_cell_no(i, j + 1)] = avg(p[i - 1][j].K[0], p[i][j].K[0], averaging) / (GP.dy*GP.dy);
		}
		// top boundary
		int j = GP.ny;
		K[get_cell_no(i, j)][get_cell_no(i, j - 1)] = 1;
		K[get_cell_no(i, j)][get_cell_no(i, j)] = -1;

		// bottom boundary
		j = 0;
		K[get_cell_no(i, j)][get_cell_no(i, j + 1)] = 1;
		K[get_cell_no(i, j)][get_cell_no(i, j)] = -1;
	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		K[get_cell_no(i, j)][get_cell_no(i, j)] = 1;
		K[get_cell_no(i, j)][get_cell_no(i + 1, j)] = 1;
	}

	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		K[get_cell_no(i, j)][get_cell_no(i, j)] = 1;
		K[get_cell_no(i, j)][get_cell_no(i - 1, j)] = 1;
	}

	return K;

}


CVector_arma CGrid::create_RHS_arma()
{
	CVector_arma V((GP.nx + 1)*(GP.ny + 1));

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
		V[get_cell_no(i, j)] = 2*leftboundary_h;


	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
		V[get_cell_no(i, j)] = 2*rightboundary_h;


	return V;
}

CVector CGrid::create_RHS()
{
	CVector V((GP.nx + 1)*(GP.ny + 1));

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
		V[get_cell_no(i, j)] = 2 * leftboundary_h;


	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
		V[get_cell_no(i, j)] = 2 * rightboundary_h;


	return V;

}


int CGrid::get_cell_no(int i, int j)
{
	return i*(GP.ny+1) + j;
}

int CGrid::get_cell_no_OU(int i, int j)
{
	return i*(GP.ny) + j;
}

CMatrix CGrid::solve()
{
    set_progress_value(0);
    CMatrix_arma_sp K = create_stiffness_matrix_arma();
    CVector_arma V = create_RHS_arma();
 //   qDebug() << "Solving the system of equations"<< endl;
    CVector_arma S = solve_ar(K, V);
//    qDebug() << "Solved"<< endl;
    H = CMatrix(GP.nx+1,GP.ny+1);
    vx = CMatrix(GP.nx, GP.ny - 1);
    vy = CMatrix(GP.nx - 1, GP.ny);
    for (int i = 0; i < GP.nx+1; i++)
        for (int j = 0; j < GP.ny+1; j++)
            H[i][j] = S[get_cell_no(i, j)];

    for (int i = 0; i < GP.nx; i++)
	for (int j = 0; j < GP.ny; j++)
	{
            double Kx1 = 0.5*(p[i][max(j - 1,0)].K[0] + p[i][j].K[0]);
            double Kx2 = 0.5*(p[min(i + 1, GP.nx - 1)][max(j - 1, 0)].K[0] + p[min(i+1,GP.nx-1)][j].K[0]);
            double Ky1 = 0.5*(p[max(i - 1, 0)][j].K[0] + p[i][j].K[0]);
            double Ky2 = 0.5*(p[i][min(j+1,GP.ny-1)].K[0] + p[min(i+1,GP.nx-1)][min(j + 1, GP.ny - 1)].K[0]);
            p[i][j].V[0] = -(0.5*Kx1 * (H[i + 1][j] - H[i][j]) / GP.dx + 0.5*Kx2 * (H[i + 1][j+1] - H[i][j+1]) / GP.dx);
            p[i][j].V[1] = -(0.5*Ky1 * (H[i][j+1] - H[i][j]) / GP.dx + 0.5*Ky2 * (H[i + 1][j + 1] - H[i+1][j]) / GP.dx);
	}

	for (int i=0; i<GP.nx; i++)
		for (int j = 0; j < GP.ny - 1; j++)
		{
			double Kx = 0.5*(p[i][j].K[0] + p[i][j + 1].K[0]);
			vx[i][j] = -Kx * (H[i + 1][j + 1] - H[i][j + 1]) / GP.dx;
		}

        set_progress_value(0.5);
	for (int i = 0; i<GP.nx-1; i++)
		for (int j = 0; j < GP.ny; j++)
		{
			double Ky = 0.5*(p[i][j].K[0] + p[i+1][j].K[0]);
			vy[i][j] = -Ky * (H[i + 1][j+1] - H[i+1][j]) / GP.dy;
		}

        set_progress_value(1);

	max_v_x = max_vx();
	min_v_x = min_vx();
	return H;
    cout<<endl;
}

vector<int> CGrid::get_ij(int k)
{
	vector<int> out(2);
	out[0] = k / (GP.ny + 1);
	out[1] = k % (GP.ny + 1);
	return out;

}

CBTC CGrid::get_v_btc(int k)
{
	CBTC out;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			out.append(i + j * 10000, p[i][j].V[k]);

	return out;
}

CBTC CGrid::get_v_mag_btc()
{
	CBTC out;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			out.append(i + j * 10000, sqrt(pow(p[i][j].V[0],2)+pow(p[i][j].V[1],2)));

	return out;
}

CBTC CGrid::get_kg_btc(int k)
{
	CBTC out;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			out.append(i + j * 10000, p[i][j].K_gauss[k]);

	return out;
}

void CGrid::remap_K(int k)
{
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			p[i][j].K[k] = map_to_KCDF(getnormalcdf(p[i][j].K_gauss[k]));

}

CBTC CGrid::get_v_btc(double x, int k)
{
	CBTC out;
	int i = int(x / GP.dx);
		for (int j = 0; j < GP.ny; j++)
			out.append(i + j * 10000, p[i][j].V[k]);
	return out;

}

CBTC CGrid::get_v_dist(double x, int k, int nbins, double smoothing_factor)
{
    CBTC v = get_v_btc(x,k);
    return v.distribution(nbins,(v.maxC()-v.minC())*smoothing_factor);
}

void CGrid::set_inv_K_dist(int ninc)
{
	inv_K_dist.clear();
	double epsilon = 1e-8;

	for (int i = 0; i < ninc+1; i++)
	{
		double u = double(i) / double(ninc)*(1 - 2 * epsilon) + epsilon;
		if (i==0)
			inv_K_dist.append(u, inv_function_s(u));
		else
			inv_K_dist.append(u, inv_function_s(u));
	}

}

double CGrid::map_to_KCDF(double u)
{
	return inv_K_dist.interpol(u);
}

double CGrid::K_CDF(double x)
{
	double out = 0;
	if (marginal_K_dist_type == "lognormal")
	{
		int n = marginal_K_dist_params.size() / 3.0;

		for (int i = 0; i < n; i++)
			out += marginal_K_dist_params[3 * i] * 0.5*(1.0 + erf((log(x) - log(marginal_K_dist_params[3 * i + 1])) / (sqrt(2.0)*marginal_K_dist_params[3 * i + 2])));
	}
	return out;

}


double CGrid::inv_function(double u, double guess)
{
	double x = guess;
	double dx = 1e-8;
	double tol = 1e-6;
	double err = K_CDF(x) - u;
	double err_p;
	double lambda = 1;
	while (fabs(err) > tol)
	{
		double dfdx = (K_CDF(x + dx) - K_CDF(x)) / dx;
		x -= lambda*err / dfdx;
		err_p = err;
		err = K_CDF(x) - u;
		if (x < 0 || !(x == x) || fabs(err)>fabs(err_p))
		{
			x += lambda*err_p / dfdx;
			lambda /= 2;
			err = err_p;
		}
		else lambda = min(lambda*1.3, 1.0);

	}
	return x;
}

double CGrid::inv_function_s(double u)
{
	double x2 = 1000;
	double x1 = 0.001;
	double tol = 1e-8;
	double err1 = K_CDF(x1) - u;
	double err2 = K_CDF(x2) - u;
	while (err1 > 0)
	{
		x1 /= 2;
		err1 = K_CDF(x1) - u;
	}

	while (err2 < 0)
	{
		x2 *= 2;
		err2 = K_CDF(x2) - u;
	}

	while (min(fabs(err1),fabs(err2)) > tol && fabs(x1-x2)>tol)
	{
		double slope = (err2 - err1) / (log(x2) - log(x1));
		double x_p = exp(log(x1) - err1 / slope);
		double err_p = K_CDF(x_p) - u;
		if (err_p > 0)
		{
			x2 = x_p;
			err2 = K_CDF(x2) - u;
		}
		else
		{
			x1 = x_p;
			err1 = K_CDF(x1) - u;
		}

	}
	if (fabs(err1) > fabs(err2))
		return x2;
	else
		return x1;
}

CBTC CGrid::get_K_CDF(double x0, double x1, double log_inc)
{
	CBTC out;
	for (double logx = log(x0); logx <= log(x1); logx += log_inc )
	{
		double x = exp(logx);
		out.append(x, K_CDF(x));
	}

	return out;
}

CBTC CGrid::get_margina_traj_v_dist(double vmin, double vmax, double nbins, string val)
{
	CBTC out(nbins + 1);
	int count = 0;
	for (int i = 0; i < int(Traj.paths.size()); i++)
	{
		for (int j = 0; j < int(Traj.paths[i].positions.size()); j++)
			out.C[int((Traj.paths[i].positions[j].getvar(val) - vmin) / nbins)] += 1;
		count += 1;
	}
	out = out/count;
	return out;
}


double avg(double x, double y, string type)
{
	if (type == "arithmetic")
            return 0.5*(x + y);
	if (type == "geometric")
            return sqrt(x*y);
	if (type == "harmonic")
            return (2.0*x*y / (x + y));
        else
            return 0.5*(x + y);

}

vector<ijval> get_top_n(vector<ijval> vec, int n)
{
	vector<ijval> out;
	vector<bool> extracted(vec.size());
	for (int i = 0; i < int(vec.size()); i++) extracted[i] = false;
	int smallest_dist = -1;

	for (int i = 0; i < n; i++)
	{
		double min_dist = 1e12;
		for (int j = 0; j < int(vec.size()); j++)
			if ((vec[j].val < min_dist) && (extracted[j]==false))
			{
				smallest_dist = j;
				min_dist = vec[j].val;
			}
		out.push_back(vec[smallest_dist]);
		extracted[smallest_dist] = true;
	}

	return out;
}

void CGrid::runcommands_qt()
{

	set_progress_value(0);
	vector<vtkSmartPointer<vtkActor>> actors;
	CBTCSet All_Breakthroughpoints;
	CBTCSet Allpoints_velocities_eulerian(3);
	CBTCSet velocity_samples_from_trajs(3);
	CBTCSet extracted_OU_parameters;

	for (int i = 0; i < commands.size(); i++)
	{
            #ifdef QT_version
            main_window->setCursor(Qt::WaitCursor);
            QApplication::processEvents();
            qDebug() << QString::fromStdString(commands[i].command) << endl;
            #endif // QT_version
            if (commands[i].command == "get_all_btc_points")
            {
                show_in_window("Getting all btc points");
                All_Breakthroughpoints.getfromfile(pathout + commands[i].parameters["filename"],1);
            }

            if (commands[i].command == "get_points_v_eulerian")
            {
                show_in_window("Getting all Eulerian velocities");
                Allpoints_velocities_eulerian.getfromfile(pathout + commands[i].parameters["filename"],1);
            }
            if (commands[i].command == "get_points_v_from_traj")
            {
                show_in_window("Getting all Lagrangian velocities");
                velocity_samples_from_trajs.getfromfile(pathout + commands[i].parameters["filename"],1);
            }
            if (commands[i].command == "get_ou_params")
            {
                show_in_window("Getting all OU parameters");
                extracted_OU_parameters.getfromfile(pathout + commands[i].parameters["filename"],1);
            }


            if (commands[i].command == "write_all_btc_points")
            {
                show_in_window("Writing all btc points");
                All_Breakthroughpoints.writetofile(pathout + commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_points_v_eulerian")
            {
                show_in_window("Writing all Eulerian velocities");
                Allpoints_velocities_eulerian.writetofile(pathout + commands[i].parameters["filename"]);
            }
            if (commands[i].command == "write_points_v_from_traj")
            {
                show_in_window("Writing all Lagrangian velocities");
                velocity_samples_from_trajs.writetofile(pathout + commands[i].parameters["filename"]);
            }
            if (commands[i].command == "write_ou_params")
            {
                 show_in_window("Writing all OU parameters");
                extracted_OU_parameters.writetofile(pathout + commands[i].parameters["filename"]);
            }


            if (commands[i].command == "generate_k_field")
            {
                clear_contents();
                show_in_window("Assigning K...,");

                //cout << "Assigning K..." << endl;
                field_gen.max_correl_n = atoi(commands[i].parameters["n_neighbors"].c_str());
                field_gen.k_correlation_lenght_scale_x = atof(commands[i].parameters["corr_length_scale_x"].c_str());
                field_gen.k_correlation_lenght_scale_y = atof(commands[i].parameters["corr_length_scale_y"].c_str());
                assign_K_gauss();
            }

            if (commands[i].command == "write_k_field")
            {
                show_in_window("Writing K field...");
                cout << "Writing K field..." << endl;
                writeasmatrixK(pathout+commands[i].parameters["filename"], 0);
            }

            if (commands[i].command == "read_k_field")
            {
                    show_in_window("Reading K field...");
                    cout << "Reading K field..." << endl;
                    GP.nx = atoi(commands[i].parameters["nx"].c_str());
                    GP.ny = atoi(commands[i].parameters["ny"].c_str());
                    GP.nx_data = atoi(commands[i].parameters["nx_data"].c_str());
                    GP.ny_data = atoi(commands[i].parameters["ny_data"].c_str());
                    getKfromfile(pathout+commands[i].parameters["filename"], GP.nx, GP.ny, GP.nx_data, GP.ny_data);
            }
            if (commands[i].command == "solve_hydro")
            {
                    show_in_window("Solving hydro ...");
                    leftboundary_h = atof(commands[i].parameters["l_boundary"].c_str());
                    rightboundary_h = atof(commands[i].parameters["r_boundary"].c_str());
                    //cout << "Solving hydro ..." << endl;
                    CMatrix H = solve();
            }

            if (commands[i].command == "solve_transport")
            {
                    show_in_window("Solving transport ...");
                    time_weight = atof(commands[i].parameters["weight"].c_str());
                    leftboundary_C = atof(commands[i].parameters["l_boundary"].c_str());
                    D = atof(commands[i].parameters["diffusion"].c_str());
                    dt = atof(commands[i].parameters["dt"].c_str());

                    cout << "Solving transport ..." << endl;
                    solve_transport(atof(commands[i].parameters["t_end"].c_str()));
            }

            if (commands[i].command == "solve_transport_ou")
            {
                show_in_window("Solving transport (Ornstein-Uhlenbeck)...");
                time_weight = atof(commands[i].parameters["weight"].c_str());
                D = atof(commands[i].parameters["diffusion"].c_str());
                dt = atof(commands[i].parameters["dt"].c_str());
                OU.kappa = atof(commands[i].parameters["kappa"].c_str());
                solve_transport_OU(atof(commands[i].parameters["t_end"].c_str()));
                if (commands[i].parameters.count("filename") > 0) OU.BTCs.writetofile(pathout + commands[i].parameters["filename"]);
                if (commands[i].parameters.count("filename_d") > 0) OU.BTCs.detivative().writetofile(pathout + commands[i].parameters["filename_d"]);
                if (commands[i].parameters.count("filename_n") > 0) OU.BTC_normal.writetofile(pathout + commands[i].parameters["filename_n"]);
                if (commands[i].parameters.count("filename_nd") > 0) OU.BTC_normal.detivative().writetofile(pathout + commands[i].parameters["filename_nd"]);
                if (commands[i].parameters.count("filename_fw_nd") > 0) OU.BTC_normal_fw.detivative().writetofile(pathout + commands[i].parameters["filename_fw_nd"]);
            }


            if (commands[i].command == "solve_transport_laplace")
            {
                show_in_window("Solving transport (laplace)...");
                leftboundary_C = 1;
                D = atof(commands[i].parameters["diffusion"].c_str());

                solve_transport_laplace(atof(commands[i].parameters["s"].c_str()));
            }

            if (commands[i].command == "write_h_field")
            {
                show_in_window("Writing H field");
               #ifdef QT_version
                QApplication::processEvents();
               #endif // QT_version
                H.writetofile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_v_field")
            {
                show_in_window("Writing velocities ...");

                //cout << "Writing velocities ..." << endl;
                writeasmatrix(pathout+commands[i].parameters["filename_x"], 0);
                writeasmatrix(pathout+commands[i].parameters["filename_y"], 1);
            }

            if (commands[i].command == "initialize_trajectories")
            {
                show_in_window("Initializing trajectories ...");
                cout << "Initializing trajectories ..." << endl;
                CBTC ini_vel=initialize(atoi(commands[i].parameters["n"].c_str()), atof(commands[i].parameters["x_0"].c_str()),atoi(commands[i].parameters["weighted"].c_str()));
                weighted = atoi(commands[i].parameters["weighted"].c_str());
                cout << "weighted = " << weighted;
                if (commands[i].parameters.count("filename"))
                    ini_vel.writefile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "create_trajectories")
            {
                show_in_window("Simulating trajectories ...");
                #ifdef QT_version
                QApplication::processEvents();
                #endif // QT_version
                cout << "Simulating trajectories ..." << endl;
                Traj = gettrajectories(atof(commands[i].parameters["dt"].c_str()), atof(commands[i].parameters["t_end"].c_str()));
            }

            if (commands[i].command == "create_trajectories_fix_dx")
            {
                show_in_window("Simulating trajectories with fixed dx...");
                #ifdef QT_version
                QApplication::processEvents();
                #endif // QT_version
                cout << "Simulating trajectories with fixed dx..." << endl;
                Traj = gettrajectories_fixed_dx(atof(commands[i].parameters["dx"].c_str()), atof(commands[i].parameters["x_end"].c_str()));
            }

            if (commands[i].command == "write_breakthrough_curve")
            {
                show_in_window("Get breakthrough curve at x = " + commands[i].parameters["x"]);
                CBTC Breakthroughcurve_from_trajs = Traj.get_BTC(atof(commands[i].parameters["x"].c_str()), atoi(commands[i].parameters["nbins"].c_str()), atof(commands[i].parameters["smoothing_factor"].c_str()));
                if (All_Breakthroughpoints.lookup("x=" + commands[i].parameters["x"]) == -1)
                {
                    CBTC BTC = Traj.get_BTC_points(atof(commands[i].parameters["x"].c_str()));
                    All_Breakthroughpoints.append(BTC, "x=" + commands[i].parameters["x"]);
                }
                else
                {   CBTC BTC = Traj.get_BTC_points(atof(commands[i].parameters["x"].c_str()));
                    All_Breakthroughpoints.BTC[All_Breakthroughpoints.lookup("x=" + commands[i].parameters["x"])].append(BTC);
                }
                Breakthroughcurve_from_trajs.writefile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_breakthrough_curves_all")
            {
                show_in_window("Writing break through curves for all realizations");
                All_Breakthroughpoints.distribution(atoi(commands[i].parameters["nbins"].c_str()), All_Breakthroughpoints.BTC.size(),0, atof(commands[i].parameters["smoothing_factor"].c_str())).writetofile(pathout + commands[i].parameters["filename"],true);
            }



            if (commands[i].command == "write_trajectories")
            {
                show_in_window("Writing trajectories ...");
                //cout << "Writing trajectories ..." << endl;
                if (commands[i].parameters.count("interval")>0)
                    Traj.write(pathout+commands[i].parameters["filename"]);
                else
                    Traj.write(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "get_velocity_dist")
            {
                show_in_window("Get Velocities into vectors...");

                CBTC vx = get_v_btc(0);
                CBTC vy = get_v_btc(1);
                CBTC v_mag = get_v_mag_btc();
                Allpoints_velocities_eulerian.BTC[0].append(vx);
                Allpoints_velocities_eulerian.BTC[1].append(vy);
                Allpoints_velocities_eulerian.BTC[2].append(v_mag);

                show_in_window("Get Velocities distributions");
                #ifdef QT_version
                QApplication::processEvents();
                #endif // QT_version
                //cout << "Get Velocities distributions" << endl;
                vx_dist = vx.distribution(atoi(commands[i].parameters["nbins"].c_str()), (vx.maxC()-vx.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                vy_dist = vy.distribution(atoi(commands[i].parameters["nbins"].c_str()), (vy.maxC()-vy.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                v_dist = v_mag.distribution(atoi(commands[i].parameters["nbins"].c_str()), (v_mag.maxC()-v_mag.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
            }

            if (commands[i].command == "write_velocity_dist")
            {
                show_in_window("Writing velocities distributions...");

                //cout << "Writing velocities into vectors" << endl;
                vx_dist.writefile(pathout+commands[i].parameters["filename_x"]);
                vy_dist.writefile(pathout+commands[i].parameters["filename_y"]);
                v_dist.writefile(pathout+commands[i].parameters["filename_mag"]);
            }

            if (commands[i].command == "write_velocity_dist_all")
            {
                show_in_window("Writing velocities distributions for all realizations...");

                CBTC vx_dist_all = Allpoints_velocities_eulerian.BTC[0].distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[0].maxC()-Allpoints_velocities_eulerian.BTC[0].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                CBTC vy_dist_all = Allpoints_velocities_eulerian.BTC[1].distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[1].maxC()-Allpoints_velocities_eulerian.BTC[1].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                CBTC v_mag_dist_all = Allpoints_velocities_eulerian.BTC[2].distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[2].maxC()-Allpoints_velocities_eulerian.BTC[2].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));

                vx_dist_all.writefile(pathout + commands[i].parameters["filename_x"]);
                vy_dist_all.writefile(pathout + commands[i].parameters["filename_y"]);
                v_mag_dist_all.writefile(pathout + commands[i].parameters["filename_mag"]);
                if (commands[i].parameters.count("filename_x_log") > 0)
                {
                    #ifdef QT_version
                    qDebug() << "Calculating log-transformed distributions" << endl;
                    #else
                    cout << "Calculating log-transformed distributions" << endl;
                    #endif // QT_version
                    CBTC vx_dist_all_log = Allpoints_velocities_eulerian.BTC[0].Log().distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[0].Log().maxC()-Allpoints_velocities_eulerian.BTC[0].Log().minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                    #ifdef QT_version
                    qDebug() << "Calculating log-transformed distributions, Done!" << endl;
                    #else
                    cout<< "Calculating log-transformed distributions, Done!" << endl;
                    #endif // QT_version
                    vx_dist_all_log.writefile(pathout + commands[i].parameters["filename_x_log"]);

                    if (commands[i].parameters.count("filename_log_int") > 0)
                    {
                        vx_dist_all_log.unlog().writefile(pathout + commands[i].parameters["filename_log_int"]);
                    }
                    #ifdef QT_version
                    qDebug() << "Calculating log-transformed distributions for velocity magnitude" << endl;
                    #else
                    cout<< "Calculating log-transformed distributions for velocity magnitude" << endl;
                    #endif // QT_version
                    CBTC v_mag_dist_all_log = Allpoints_velocities_eulerian.BTC[2].Log().distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[2].Log().maxC()-Allpoints_velocities_eulerian.BTC[2].Log().minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                    #ifdef QT_version
                    qDebug() << "Calculating log-transformed distributions for velocity magnitude, Done!" << endl;
                    #else
                    cout << "Calculating log-transformed distributions for velocity magnitude, Done!" << endl;
                    #endif // QT_version
                    v_mag_dist_all_log.writefile(pathout + commands[i].parameters["filename_mag_log"]);
                    if (commands[i].parameters.count("filename_log_int_mag") > 0)
                    {
                        v_mag_dist_all_log.unlog().writefile(pathout + commands[i].parameters["filename_log_int_mag"]);
                    }
                }
            }

            if (commands[i].command == "write_velocity_dist_all_fw")
            {
                show_in_window("Writing velocities distributions for all realizations (flux weighted)...");

                CBTC vx_dist_all = Allpoints_velocities_eulerian.BTC[0].distribution_fw(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[0].maxC()-Allpoints_velocities_eulerian.BTC[0].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                CBTC vy_dist_all = Allpoints_velocities_eulerian.BTC[1].distribution_fw(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[1].maxC()-Allpoints_velocities_eulerian.BTC[1].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                CBTC v_mag_dist_all = Allpoints_velocities_eulerian.BTC[2].distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[2].maxC()-Allpoints_velocities_eulerian.BTC[2].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                vx_dist_all.writefile(pathout + commands[i].parameters["filename_x"]);
                vy_dist_all.writefile(pathout + commands[i].parameters["filename_y"]);
                v_mag_dist_all.writefile(pathout + commands[i].parameters["filename_mag"]);
                if (commands[i].parameters.count("filename_x_log") > 0)
                {

                    CBTC vx_dist_all_log = Allpoints_velocities_eulerian.BTC[0].Log().distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[0].Log().maxC()-Allpoints_velocities_eulerian.BTC[0].Log().minC())*atof(commands[i].parameters["smoothing_factor"].c_str())).make_flux_weighted("log");
                    vx_dist_all_log.writefile(pathout + commands[i].parameters["filename_x_log"]);
                    if (commands[i].parameters.count("filename_log_int") > 0)
                    {
                        vx_dist_all_log.unlog().writefile(pathout + commands[i].parameters["filename_log_int"]);
                    }

                    CBTC v_mag_dist_all_log = Allpoints_velocities_eulerian.BTC[2].Log().distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[2].Log().maxC()-Allpoints_velocities_eulerian.BTC[2].Log().minC())*atof(commands[i].parameters["smoothing_factor"].c_str())).make_flux_weighted("log");
                    v_mag_dist_all_log.writefile(pathout + commands[i].parameters["filename_mag_log"]);
                    if (commands[i].parameters.count("filename_log_int_mag") > 0)
                    {
                        v_mag_dist_all_log.unlog().writefile(pathout + commands[i].parameters["filename_log_int_mag"]);
                    }
                }
            }

            if (commands[i].command == "get_velocity_dist_at_sects")
            {
                show_in_window("Get Velocities distributions at cross-sections ...");

                //cout << "Get Velocities distributions at cross-sections ..." << endl;

                CBTC vx = get_v_btc(atof(commands[i].parameters["x"].c_str()), 0);
                CBTC vy = get_v_btc(atof(commands[i].parameters["x"].c_str()), 1);
                vx_dist = vx.distribution(atoi(commands[i].parameters["nbins"].c_str()), (vx.maxC()-vx.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                vy_dist = vy.distribution(atoi(commands[i].parameters["nbins"].c_str()), (vy.maxC()-vy.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                sect_dist.append(vx_dist, "x=" + commands[i].parameters["x"]);
                sect_dist.append(vy_dist, "x=" + commands[i].parameters["x"]);

            }

            if (commands[i].command == "write_velocities_at_sects")
            {
                show_in_window("Write velocities distributions at cross-sections ...");

                //cout << "Write velocities distributions at cross-sections ..." << endl;
                sect_dist.writetofile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "create_inverse_marginal_k")
            {
                show_in_window("Create inverse marginal hydraulic conductivity ...");

                //cout << "Create inverse marginal hydraulic conductivity ..." << endl;
                set_inv_K_dist(atoi(commands[i].parameters["ninc"].c_str()));
            }

            if (commands[i].command == "write_inverse_marginal_k")
            {
                show_in_window("Write inverse marginal hydraulic conductivity ...");

                //cout << "Write inverse marginal hydraulic conductivity ..." << endl;
                inv_K_dist.writefile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_marginal_k")
            {
                show_in_window("Write marginal hydraulic conductivity ...");

                //cout << "Write marginal hydraulic conductivity ..." << endl;
                CBTC K = get_K_CDF(atof(commands[i].parameters["x0"].c_str()), atof(commands[i].parameters["x1"].c_str()), atof(commands[i].parameters["log_inc"].c_str()));
                K.writefile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "renormalize_k")
            {
                show_in_window("Renormalizing hydraulic conductivity ...");

                //cout << "Renormalizing hydraulic conductivity ..." << endl;
                renormalize_k();
            }

            #ifdef QT_version
            if (commands[i].command == "show_k_field")
            {
                show_in_window("Showing K field...");
                show_K_field();
            }
            #endif // QT_version

            if (commands[i].command == "save_k_field_into_vtp")
            {
                show_in_window("Saving K field (VTK)...");
                if (commands[i].parameters.count("filename") == 0) commands[i].parameters["filename"] == "surface.vtp";
                write_K_field_to_vtp(pathout+commands[i].parameters["filename"], atof(commands[i].parameters["z_factor"].c_str()),atoi(commands[i].parameters["offset"].c_str()));
            }

            if (commands[i].command == "save_solution_into_vtp")
            {

                show_in_window("Saving solution (VTK)...");
                if (commands[i].parameters.count("filename") == 0) commands[i].parameters["filename"] == "surface.vtp";
                write_K_solution_to_vtp(pathout + commands[i].parameters["filename"], atof(commands[i].parameters["z_factor"].c_str()), atoi(commands[i].parameters["log"].c_str()));
            }

            if (commands[i].command == "save_concentration_into_vtp")
            {
                double concentration_interval = -1;
                show_in_window("Saving Concentration (VTK)...");
                if (commands[i].parameters.count("interval") > 0) concentration_interval = atof(commands[i].parameters["[c_interval"].c_str());
                vector<double> intervals;
                for (int i = 0; i < p[0][0].C.size(); i++)
                {
                        intervals.push_back(i*dt);
                }
                if (commands[i].parameters.count("filename") == 0) commands[i].parameters["filename"] == "surface.vtp";
                write_C_to_vtp(pathout + commands[i].parameters["filename"], atof(commands[i].parameters["z_factor"].c_str()), atoi(commands[i].parameters["log"].c_str()), intervals);
            }

            if (commands[i].command == "show")
            {
/*                show_in_window("Showing so far...");

                if (commands[i].parameters.count("filename") == 0)
                        showthings(actors);
                else if (commands[i].parameters.size() == 1)
                        showthings(actors,pathout+commands[i].parameters["filename"]);
                actors.clear();
                */
            }

            if (commands[i].command == "get_trajectory")
            {
                show_in_window("Showing trajectory " + commands[i].parameters[0]);
                actors.push_back(traj_vtk_pdt(atoi(commands[i].parameters["traj_no"].c_str()), atof(commands[i].parameters["z_factor"].c_str()), atof(commands[i].parameters["offset"].c_str())));
            }

            if (commands[i].command == "get_trajectories")
            {
                show_in_window("Converting trajectories into graphic objects... ");
                vector<vtkSmartPointer<vtkActor>> tmp_actor;
                tmp_actor = trajs_vtk_pdt(atof(commands[i].parameters["z_factor"].c_str()), atof(commands[i].parameters["offset"].c_str()));

                for (int i = 0; i < tmp_actor.size(); i++)
                        actors.push_back(tmp_actor[i]);

            }

            #if QT_version
            if (commands[i].command == "screen_shot_test")
            {
                show_in_window("Screen Shot test");
                screenshot_test();
            }
            #endif // QT_version

            if (commands[i].command == "save_trajs_into_vtp")
            {
                show_in_window("Writing trajectories into vtp... ");
                if (commands[i].parameters.count("filename") == 0) commands[i].parameters["filename"] == "paths.vtp";
                trajs_vtk_pdt_to_vtp(pathout+commands[i].parameters["filename"], atof(commands[i].parameters["z_factor"].c_str()), atof(commands[i].parameters["offset"].c_str()), atof(commands[i].parameters["log"].c_str()), atof(commands[i].parameters["color"].c_str()));
            }

            if (commands[i].command == "initialize_marginal_v_dist")
            {
                show_in_window("initializing " +commands[i].parameters["type"] + " distribution ..." );
                dist.name = commands[i].parameters["type"];
                for (int j = 0; j < split(commands[i].parameters["values"],',').size(); j++)
                        dist.params.push_back(atof(split(commands[i].parameters["values"],',')[j].c_str()));

            }

            if (commands[i].command == "create_correlated_v_pathways")
            {
                show_in_window("Creating pathways ...");
                pset.create(atoi(commands[i].parameters["n"].c_str()), &dist, atof(commands[i].parameters["x_min"].c_str()), atof(commands[i].parameters["x_max"].c_str()), atof(commands[i].parameters["kappa"].c_str()), atof(commands[i].parameters["dx"].c_str()));

            }

            if (commands[i].command == "write_pathways")
            {
                show_in_window("Writing pathways ...");
                pset.write(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_pathways_vtp")
            {
                show_in_window("Writing pathways to vtp ...");
                vtkSmartPointer<vtkPolyDataMapper> mapper;
                mapper = pset.pathways_vtk_pdt_vtp(atof(commands[i].parameters["z_factor"].c_str()), atof(commands[i].parameters["offset"].c_str()));
                pset.write_vtk(mapper, pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "get_dist_at_time")
            {
                show_in_window("Getting distributions at time");
                dist_stores = pset.snapshotattime(atof(commands[i].parameters["t"].c_str()));
            }

            if (commands[i].command == "get_dist_at_location")
            {
                show_in_window("Getting distributions at location");
                dist_stores = pset.snapshotatlocation(atof(commands[i].parameters["x"].c_str()));
            }

            if (commands[i].command == "write_dist")
            {
                show_in_window("Writing distribution ...");
                dist_stores.write(pathout+commands[i].parameters["filename"].c_str());
            }

            if (commands[i].command == "write_dist_function")
            {
                show_in_window("writing distribution curve ...");
                if (commands[i].parameters.count("var") > 0)
                {
                    CBTC dist_curve = dist_stores.get_distribution(commands[i].parameters["var"], atoi(commands[i].parameters["log"].c_str()), atoi(commands[i].parameters["nbins"].c_str()));
                    dist_curve.writefile(pathout+commands[i].parameters["filename"]);
                }
                else
                {
                    CBTCSet dist_curve = dist_stores.get_distribution(atoi(commands[i].parameters["log"].c_str()), atoi(commands[i].parameters["nbins"].c_str()));
                    dist_curve.writetofile(pathout+commands[i].parameters["filename"]);
                }


            }

            if (commands[i].command == "create_grid")
            {
                creategrid(atoi(commands[i].parameters["nx"].c_str()), atoi(commands[i].parameters["ny"].c_str()), atof(commands[i].parameters["dx"].c_str()), atof(commands[i].parameters["dy"].c_str()));
            }

            if (commands[i].command == "set_marginal_k_dist")
            {

                marginal_K_dist_type = commands[i].parameters["dist"];
                for (int j = 0; j < split(commands[i].parameters["params"],',').size(); j++)
                        marginal_K_dist_params.push_back(atof(split(commands[i].parameters["params"], ',')[j].c_str()));

            }

            if (commands[i].command == "make_paths_uniform")
            {
                show_in_window("Making trajectories uniform ... ");
                Traj.make_uniform_at_x(atof(commands[i].parameters["dx"].c_str()));
            }

            if (commands[i].command == "extract_pairs")
            {
                show_in_window("Extracting pairs ... ");
                CBTCSet pairs;
                if (!commands[i].parameters.count("nsequence"))
                    commands[i].parameters["nsequence"] = "2";
                if (Traj.paths.size())
                    pairs = Traj.get_pair_v(atoi(commands[i].parameters["increment"].c_str()), atoi(commands[i].parameters["n"].c_str()),atoi(commands[i].parameters["nsequence"].c_str()));
                else if (pset.paths.size())
                    pairs = pset.get_pair_v(atoi(commands[i].parameters["increment"].c_str()), atoi(commands[i].parameters["n"].c_str()),atoi(commands[i].parameters["nsequence"].c_str()));
                else
                {
                    show_in_window("No trajectories has been initialized!");
                    return;
                }

                pairs.writetofile(pathout+commands[i].parameters["filename"]);
                if (commands[i].parameters.count("dist_filename") > 0)
                {
                        pairs.BTC[0].distribution(atoi(commands[i].parameters["nbins"].c_str()),(pairs.BTC[0].maxC()-pairs.BTC[0].minC())*atof(commands[i].parameters["smoothing_factor"].c_str())).writefile(pathout+commands[i].parameters["dist_filename"]);
                }
                if (commands[i].parameters.count("v_gnu_file"))
                {
                    TDMap GNU_out = pairs.get2DMap(atoi(commands[i].parameters["nbins"].c_str()));
                    GNU_out.writetofile_GNU(pathout + commands[i].parameters["v_gnu_file"],"", "v", "v'", "p(v,v')", true);
                }
                if (commands[i].parameters.count("ranks_filename") > 0)
                {
                    show_in_window("Writing ranks");
                    CBTCSet ranks(atoi(commands[i].parameters["nsequence"].c_str()));
                    for (int ii=0; ii<atoi(commands[i].parameters["nsequence"].c_str()); ii++)
                    {
                        cout<<pairs.BTC[ii].n<<endl;
                        pairs.BTC[ii].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout+commands[i].parameters["ranks_filename"]+"_"+numbertostring(ii));
                        cout<<pairs.BTC[ii].n<<endl;
                        ranks.BTC[ii] = pairs.BTC[ii].rank_bd(atoi(commands[i].parameters["nbins"].c_str()));
                        cout<<ranks.BTC[ii].n<<endl;
                        ranks.BTC[ii].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout+commands[i].parameters["ranks_filename"]+"_u"+numbertostring(ii));
                    }
                    ranks.writetofile(pathout+commands[i].parameters["ranks_filename"]);

                    if (commands[i].parameters.count("u_gnu_file"))
                    {
                        TDMap GNU_out = ranks.get2DMap(atoi(commands[i].parameters["nbins"].c_str()),0,1);
                        GNU_out.writetofile_GNU(pathout + commands[i].parameters["u_gnu_file"],"", "u", "u'", "p(u,u')");
                    }
                }
                if (commands[i].parameters.count("normal_filename") > 0)
                {
                    show_in_window("Writing normals");
                    CBTCSet normals(atoi(commands[i].parameters["nsequence"].c_str()));
                    for (int ii=0; ii<atoi(commands[i].parameters["nsequence"].c_str()); ii++)
                        normals.BTC[ii] = pairs.BTC[ii].map_to_standard_normal(atoi(commands[i].parameters["nbins"].c_str()));

                    if (commands[i].parameters.count("w_gnu_file"))
                    {
                        TDMap GNU_out = normals.get2DMap(atoi(commands[i].parameters["nbins"].c_str()),-4,4);
                        GNU_out.writetofile_GNU(pathout + commands[i].parameters["w_gnu_file"],"", "w", "w'", "p(w,w')");
                    }

                    normals.writetofile(pathout + commands[i].parameters["normal_filename"]);
                    if (commands[i].parameters.count("OU_parameters_filename") > 0)
                    {
                        show_in_window("Calculating OU params");

                        CVector X = normals.get_kappa_gamma(atof(commands[i].parameters["delta_x"].c_str()));
                        extracted_OU_parameters.append("p1_"+commands[i].parameters["increment"], atoi(commands[i].parameters["increment"].c_str()), X[0]);
                        extracted_OU_parameters.append("p2_"+commands[i].parameters["increment"], atoi(commands[i].parameters["increment"].c_str()), X[1]);
                        extracted_OU_parameters.append("p3_"+commands[i].parameters["increment"], atoi(commands[i].parameters["increment"].c_str()), X[2]);
                        show_in_window("Writing OU params");
                        X.writetofile(pathout + commands[i].parameters["OU_parameters_filename"]);
                    }
                }
            }

            if (commands[i].command == "write_all_ou_parameters")
            {
                show_in_window("Writing OU parameters ... ");
                extracted_OU_parameters.writetofile(pathout+commands[i].parameters["filename"]);
                if (commands[i].parameters.count("density_filename") > 0)
                {
                        CBTCSet OU_Density;
                        for (int j = 0; j < extracted_OU_parameters.nvars; j++)
                        {   CBTC BTC = extracted_OU_parameters.BTC[j].distribution(atoi(commands[i].parameters["nbins"].c_str()),0);
                            OU_Density.append(BTC,extracted_OU_parameters.BTC[j].name);

                        }

                        OU_Density.writetofile(pathout+commands[i].parameters["density_filename"]);
                }
            }

            if (commands[i].command == "clear_all")
            {
                show_in_window("Clearing variables ... ");
                C.matr.clear();
                H.matr.clear();
                KD.matr.clear();
                Kt.matr.clear();
                Kv.matr.clear();
                clear();
                trajectories.clear();
                Traj.paths.clear();
                sect_dist.clear();
                pts.clear();
            }

            if (commands[i].command == "write")
            {
                    show_in_window(commands[i].parameters["content"]);
            }

	}

	show_in_window("Finished!");
#ifdef QT_version
	main_window->setCursor(Qt::ArrowCursor);
#endif // QT_version
}




double CGrid::max_K()
{
	double max_k = -1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			max_k = max(max_k, p[i][j].K[0]);

	return max_k;

}
double CGrid::min_K()
{
	double min_k = 1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			min_k = min(min_k, p[i][j].K[0]);
	return min_k;
}
double CGrid::max_vx()
{
	double max_k = -1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			max_k = max(max_k, p[i][j].V[0]);
	return max_k;
}

double CGrid::min_vx()
{
	double min_k = 1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			min_k = min(min_k, p[i][j].V[0]);
	return min_k;
}
double CGrid::max_vy()
{
	double max_k = -1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			max_k = max(max_k, p[i][j].V[1]);
	return max_k;

}
double CGrid::min_vy()
{
	double min_k = 1e23;
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			min_k = min(min_k, p[i][j].V[1]);
	return min_k;
}

void CGrid::set_K_transport(double dt, double D, double weight)
{
	string averaging = "arithmetic";
	Kv = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	KD = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	Kt = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));

	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			Kv.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -weight*pos(vx[i-1][j-1]) / (GP.dx);
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += weight*(neg(vx[i-1][j-1])/GP.dx + pos(vx[i][j-1])/GP.dx);
			Kv.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -weight*neg(vx[i][j - 1]) / (GP.dx);

			Kv.matr(get_cell_no(i, j), get_cell_no(i, j-1)) += -weight*pos(vy[i - 1][j - 1]) / GP.dy;
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += weight*(neg(vy[i - 1][j - 1]) / GP.dy + pos(vy[i-1][j]) / GP.dy);
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j-1)) += -weight*neg(vy[i-1][j]) / GP.dy;

			KD.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -weight*D / (GP.dx*GP.dx);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2*weight*D / (GP.dx*GP.dx);
			KD.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -weight*D / (GP.dx*GP.dx);

			KD.matr(get_cell_no(i, j), get_cell_no(i, j-1)) += -weight*D / (GP.dy*GP.dy);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2*weight*D / (GP.dy*GP.dy);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j+1)) += -weight*D / (GP.dy*GP.dy);

			Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) += 1.0 / dt;

		}
		// top boundary
		int j = GP.ny;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;

		// bottom boundary
		j = 0;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;
	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) = 1;
	}

	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = -1;
	}


}

CVector_arma CGrid::create_RHS_transport(double dt, double weight,double D)
{
	CVector_arma RHS((GP.nx + 1)*(GP.ny + 1));
	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			double rhs = 0;
			rhs += (1-weight)*pos(vx[i - 1][j - 1]) / (GP.dx)*C[i-1][j];
			rhs += -(1-weight)*(neg(vx[i - 1][j - 1]) / GP.dx + pos(vx[i][j-1]) / GP.dx)*C[i][j];
			rhs += (1-weight)*neg(vx[i][j - 1]) / (GP.dx)*C[i+1][j];

			rhs += (1-weight)*pos(vy[i - 1][j - 1]) / GP.dy*C[i][j-1];
			rhs += -(1 - weight)*(neg(vy[i - 1][j - 1]) / GP.dy + pos(vy[i - 1][j]) / GP.dy)*C[i][j];
			rhs += (1-weight)*neg(vy[i - 1][j]) / GP.dy*C[i][j+1];

			rhs += (1-weight)*D / (GP.dx*GP.dx)*C[i-1][j];
			rhs += -2 * (1-weight)*D / (GP.dx*GP.dx)*C[i][j];
			rhs += (1-weight)*D / (GP.dx*GP.dx)*C[i + 1][j];

			rhs += (1-weight)*D / (GP.dy*GP.dy)*C[i][j-1];
			rhs += -2 * (1-weight)*D / (GP.dy*GP.dy)*C[i][j];
			rhs += (1-weight)*D / (GP.dy*GP.dy)*C[i][j+1];

			rhs += 1.0 / dt*C[i][j];
			RHS[get_cell_no(i, j)] = rhs;
		}
		// top boundary
		int j = GP.ny;
		RHS[get_cell_no(i, j)] = 0;


		// bottom boundary
		j = 0;
		RHS[get_cell_no(i, j)] = 0;

	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
		RHS[get_cell_no(i, j)] = 2*leftboundary_C;


	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
		RHS[get_cell_no(i, j)] = 0;

	return RHS;

}

CVector_arma CGrid::create_RHS_transport_laplace(double weight, double D, double s)
{
	CVector_arma RHS((GP.nx + 1)*(GP.ny + 1));
	for (int i = 1; i < GP.nx; i++)
	{
		// top boundary
		int j = GP.ny;
		RHS[get_cell_no(i, j)] = 0;


		// bottom boundary
		j = 0;
		RHS[get_cell_no(i, j)] = 0;

	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
		RHS[get_cell_no(i, j)] = 2 * leftboundary_C;


	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
		RHS[get_cell_no(i, j)] = 0;

	return RHS;

}

void CGrid::solve_transport(double t_end)
{
	set_K_transport(dt, D, time_weight);
	C = CMatrix(GP.nx + 1, GP.ny + 1);// = leftboundary_C;
	CMatrix_arma_sp K = KD + Kt + Kv;
	Kt.writetofile(pathout + "Kt_matrix.txt");
	KD.writetofile(pathout + "KD_matrix.txt");
	Kv.writetofile(pathout + "Kv_matrix.txt");
	K.writetofile(pathout + "transport_matrix.txt");
	set_progress_value(0);
        for (double t = 0; t < t_end; t += dt)
        {
            CVector_arma RHS = create_RHS_transport(dt, time_weight, D);
            CVector_arma S = solve_ar(K, RHS);

            for (int i=0; i<GP.nx+1; i++)
                for (int j=0; j<GP.ny+1; j++)
                    C[i][j] = S[get_cell_no(i, j)];

            for (int i = 0; i < GP.nx; i++)
                for (int j = 0; j < GP.ny; j++)
                    p[i][j].C.push_back(0.25*(C[i][j] + C[i + 1][j] + C[i][j + 1] + C[i + 1][j + 1]));

            set_progress_value(t / t_end);
    //                tbrowse->append("t = " + QString::number(t));

        }
        cout<<endl;
}

void CGrid::solve_transport_laplace(double s)
{
	set_K_transport_laplace(D, s);
	C = CMatrix(GP.nx + 1, GP.ny + 1);// = leftboundary_C;
	CMatrix_arma_sp K = KD + Kt + Kv;


	CVector_arma RHS = create_RHS_transport_laplace(dt, time_weight, D);
	CVector_arma S = solve_ar(K, RHS);

	for (int i = 0; i<GP.nx + 1; i++)
		for (int j = 0; j<GP.ny + 1; j++)
			C[i][j] = S[get_cell_no(i, j)];

	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			p[i][j].C.push_back(0.25*(C[i][j] + C[i + 1][j] + C[i][j + 1] + C[i + 1][j + 1]));


}


void CGrid::set_K_transport_laplace(double D, double s)
{
	string averaging = "arithmetic";
	Kv = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	KD = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));
	Kt = CMatrix_arma_sp((GP.nx + 1)*(GP.ny + 1), (GP.nx + 1)*(GP.ny + 1));

	for (int i = 1; i < GP.nx; i++)
	{
		for (int j = 1; j < GP.ny; j++)
		{
			Kv.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -pos(vx[i - 1][j - 1]) / (GP.dx);
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += (neg(vx[i - 1][j - 1]) / GP.dx + pos(vx[i][j - 1]) / GP.dx);
			Kv.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -neg(vx[i][j - 1]) / (GP.dx);

			Kv.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) += -pos(vy[i - 1][j - 1]) / GP.dy;
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j)) += (neg(vy[i - 1][j - 1]) / GP.dy + pos(vy[i - 1][j]) / GP.dy);
			Kv.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) += -neg(vy[i - 1][j]) / GP.dy;

			KD.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) += -D / (GP.dx*GP.dx);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2 * D / (GP.dx*GP.dx);
			KD.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) += -D / (GP.dx*GP.dx);

			KD.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) += -D / (GP.dy*GP.dy);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j)) += 2 * D / (GP.dy*GP.dy);
			KD.matr(get_cell_no(i, j), get_cell_no(i, j+1)) += -D / (GP.dy*GP.dy);

			Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) += s;

		}
		// top boundary
		int j = GP.ny;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j - 1)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;

		// bottom boundary
		j = 0;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j + 1)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = -1;
	}

	//left boundary
	int i = 0;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i + 1, j)) = 1;
	}

	//right boundary
	i = GP.nx;
	for (int j = 0; j < GP.ny + 1; j++)
	{
		Kt.matr(get_cell_no(i, j), get_cell_no(i, j)) = 1;
		Kt.matr(get_cell_no(i, j), get_cell_no(i - 1, j)) = -1;
	}


}


void CGrid::create_inverse_K_OU(double dt)
{
	CMatrix_arma M(GP.ny*(GP.nx+2), GP.ny*(GP.nx+2));

	for (int i = 1; i < GP.nx+1; i++)
	{
		for (int j = 0; j < GP.ny; j++)
		{
			// Advection
			M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) = 1.0 / dt;
			M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*OU.FinvU[j] / GP.dx;
			M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i - 1, j)) += -time_weight*OU.FinvU[j] / GP.dx;

			/*if (i < GP.nx - 1)
			{
				M(get_cell_no(i, j), get_cell_no(i, j)) = 2 * time_weight*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
				M(get_cell_no(i, j), get_cell_no(i-1, j)) = -w*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
				M(get_cell_no(i, j), get_cell_no(i+1, j)) = -w*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
			}*/

			//Diffusion
			if (i < GP.nx+ 1)
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += 2 * time_weight*D / pow(GP.dx, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i - 1, j)) += -time_weight*D / pow(GP.dx, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i + 1, j)) += -time_weight*D / pow(GP.dx, 2);
			}
			else
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*D / pow(GP.dx, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i - 1, j)) += -time_weight*D / pow(GP.dx, 2);
			}
			//Exchange
			if (j > 0 && j < GP.ny - 1)
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*OU.FinvU[j] * OU.kappa*(pow(std_normal_phi_inv((double(j))*GP.dy), 2) + pow(std_normal_phi_inv((double(j) + 1)*GP.dy), 2))/pow(GP.dy,2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j - 1)) -= time_weight*OU.FinvU[j] * OU.kappa*pow(std_normal_phi_inv((double(j))*GP.dy), 2) / pow(GP.dy, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j + 1)) -= time_weight*OU.FinvU[j] * OU.kappa*pow(std_normal_phi_inv((double(j) + 1)*GP.dy), 2) / pow(GP.dy, 2);
			}
			else if (j == 0)
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*OU.FinvU[j] * OU.kappa*pow(std_normal_phi_inv((double(j) + 1)*GP.dy), 2) / pow(GP.dy, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j + 1)) -= time_weight*OU.FinvU[j] * OU.kappa*pow(std_normal_phi_inv((double(j) + 1)*GP.dy), 2) / pow(GP.dy, 2);
			}
			else if (j == GP.ny - 1)
			{
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) += time_weight*OU.FinvU[j] * OU.kappa*pow(std_normal_phi_inv((double(j))*GP.dy), 2) / pow(GP.dy, 2);
				M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j - 1)) -= time_weight*OU.FinvU[j] * OU.kappa*pow(std_normal_phi_inv((double(j))*GP.dy), 2) / pow(GP.dy, 2);

			}
		}

	}

	int i = 0;
	for (int j = 0; j < GP.ny; j++)
	{
		M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) = 1;
		M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i+1, j)) = 1;
	}

	i = GP.nx+1;
	for (int j = 0; j < GP.ny; j++)
	{
		M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i, j)) = 1;
		M.matr(get_cell_no_OU(i, j), get_cell_no_OU(i - 1, j)) = -1;
	}

	OU.Inv_M = inv(M);
	CMatrix MM_inv(OU.Inv_M);

}


void CGrid::create_f_inv_u()
{
	OU.FinvU = CVector(GP.ny);
	for (int j = 0; j < GP.ny; j++)
	{
		double u1 = double(j)*GP.dy + GP.dy / 2;
		OU.FinvU[j] = dist.inverseCDF(u1);
	}
}

CVector_arma CGrid::create_RHS_OU(double dt)
{
	CVector_arma RHS(GP.ny*(GP.nx+2));

	for (int i = 1; i < GP.nx+1; i++)
	{
		for (int j = 0; j < GP.ny; j++)
		{
			// Advection
			RHS[get_cell_no_OU(i, j)] = 1.0 / dt*C[i][j];
			RHS[get_cell_no_OU(i, j)] += -(1-time_weight)*OU.FinvU[j] / GP.dx*C[i][j];
			RHS[get_cell_no_OU(i, j)] += (1-time_weight)*OU.FinvU[j] / GP.dx*C[i-1][j];

			/*if (i < GP.nx - 1)
			{
			M(get_cell_no(i, j), get_cell_no(i, j)) = 2 * time_weight*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
			M(get_cell_no(i, j), get_cell_no(i-1, j)) = -w*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
			M(get_cell_no(i, j), get_cell_no(i+1, j)) = -w*pow(FinvU[j], 2)*dt / 2 / pow(dx, 2);
			}*/

			//Diffusion
			if (i < GP.nx - 1)
			{
				RHS[get_cell_no_OU(i, j)] -= 2 * (1-time_weight)*D / pow(GP.dx, 2)*C[i][j];
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*D / pow(GP.dx, 2)*C[i-1][j];
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*D / pow(GP.dx, 2)*C[i+1][j];
			}
			else
			{
				RHS[get_cell_no_OU(i, j)] += -(1-time_weight)*D / pow(GP.dx, 2)*C[i][j];
				RHS[get_cell_no_OU(i, j)] += (1-time_weight)*D / pow(GP.dx, 2)*C[i-1][j];
			}
			//Exchange
			if (j > 0 && j < GP.ny - 1)
			{
				RHS[get_cell_no_OU(i, j)] += -OU.FinvU[j]*(1-time_weight)*OU.kappa*(pow(std_normal_phi_inv((double(j))*GP.dy), 2) + pow(std_normal_phi_inv((double(j) + 1)*GP.dy), 2))*C[i][j] / pow(GP.dy, 2);
				RHS[get_cell_no_OU(i, j)] += OU.FinvU[j]*(1-time_weight)*OU.kappa*pow(std_normal_phi_inv((double(j))*GP.dy), 2)*C[i][j-1] / pow(GP.dy, 2);
				RHS[get_cell_no_OU(i, j)] += OU.FinvU[j]*(1-time_weight)*OU.kappa*pow(std_normal_phi_inv((double(j) + 1)*GP.dy), 2)*C[i][j+1] / pow(GP.dy, 2);
			}
			else if (j == 0)
			{
				RHS[get_cell_no_OU(i, j)] += -(1 - time_weight)*OU.FinvU[j]*OU.kappa*(pow(std_normal_phi_inv((double(j) + 1)*GP.dy), 2))*C[i][j] / pow(GP.dy, 2);
				RHS[get_cell_no_OU(i, j)] += (1 - time_weight)*OU.FinvU[j]*OU.kappa*pow(std_normal_phi_inv((double(j) + 1)*GP.dy), 2)*C[i][j + 1] / pow(GP.dy, 2);
			}
			else if (j == GP.ny - 1)
			{
				RHS[get_cell_no_OU(i, j)] += -(1 - time_weight)*OU.FinvU[j]*OU.kappa*(pow(std_normal_phi_inv((double(j))*GP.dy), 2))*C[i][j] / pow(GP.dy, 2);
				RHS[get_cell_no_OU(i, j)] += (1 - time_weight)*OU.FinvU[j] *OU.kappa*pow(std_normal_phi_inv((double(j))*GP.dy), 2)*C[i][j - 1] / pow(GP.dy, 2);
			}
		}

	}

	int i = 0;
	for (int j = 0; j < GP.ny; j++)
	{
		RHS[get_cell_no_OU(i, j)] = 2;
	}

	i = GP.nx+1;
	for (int j = 0; j < GP.ny; j++)
	{
		RHS[get_cell_no_OU(i, j)] = 0;
	}
	return RHS;

}

void CGrid::solve_transport_OU(double t_end)
{
	create_f_inv_u();
	create_inverse_K_OU(dt);
	C = CMatrix(GP.nx+2, GP.ny);
	OU.BTCs = CBTCSet(GP.nx+2);
	OU.BTC_normal = CBTCSet(GP.nx + 2);
	OU.BTCs_fw = CBTCSet(GP.nx + 2);
	OU.BTC_normal_fw = CBTCSet(GP.nx + 2);
	for (int i = 0; i < GP.nx+2; i++) OU.BTCs.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTCs_fw.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal_fw.names[i] = ("x=" + numbertostring((i - 0.5)*GP.dx));
	//K.writetofile(pathout + "transport_matrix.txt");
	for (int i = 0; i < GP.nx + 2; i++) OU.BTCs.BTC[i].append(0, 0);
	for (int i = 0; i < GP.nx + 2; i++) OU.BTC_normal.BTC[i].append(0, 0);
	set_progress_value(0);
    for (double t = 0; t < t_end; t += dt)
    {
        CVector_arma RHS = create_RHS_OU(dt);
        CVector_arma S = OU.Inv_M*RHS;

        for (int i = 0; i < GP.nx+2; i++)
        {
            double sum = 0;
            double sum_fw = 0;

            for (int j = 0; j < GP.ny; j++)
            {
                C[i][j] = S[get_cell_no_OU(i, j)];
                sum += C[i][j] * GP.dy;
                sum_fw += C[i][j] * OU.FinvU[j]/OU.FinvU.sum();
            }
            OU.BTCs.BTC[i].append(t+dt, sum);
            OU.BTCs_fw.BTC[i].append(t + dt, sum_fw);
            OU.BTC_normal.BTC[i].append((t + dt)/((i-0.5)*GP.dx), sum);
            OU.BTC_normal_fw.BTC[i].append((t + dt) / ((i - 0.5)*GP.dx), sum_fw);
        }

        #if QT_version
                set_progress_value(t / t_end);
        tbrowse->append("t = " + QString::number(t));
        #else
        cout<<"t = " <<t <<endl;
        #endif // QT_version

    }

}

void CGrid::renormalize_k()
{
	CBTC K = get_kg_btc(0);
	double mu = K.mean();
	double std = K.std();
	for (int i = 0; i < GP.nx; i++)
		for (int j = 0; j < GP.ny; j++)
			p[i][j].K_gauss[0] = (p[i][j].K_gauss[0] - mu) / std;
	remap_K(0);
}


void CGrid::show_in_window(string s)
{
    #ifdef QT_version
    qDebug()<<QString::fromStdString(s);
    main_window->get_ui()->ShowOutput->append(QString::fromStdString(s));
    QApplication::processEvents();
    #else
    cout<<s<<endl;
    #endif // Qt_version
}

void CGrid::set_progress_value(double s)
{
#ifdef QT_version
	main_window->get_ui()->progressBar->setValue(s*100);
	QApplication::processEvents();
#endif // QT_version
    cout << "\r Progress: " << s*100 << "%";
}

void CGrid::clear_contents()
{
	this->n_k_dets = 0;
	p.resize(GP.nx);
	for (int i = 0; i < GP.nx; i++)
	{
		p[i].resize(GP.ny);
		for (int j = 0; j < GP.ny; j++)
		{
			p[i][j].k_det = false;
			p[i][j].V = CVector(2);
			p[i][j].K = CVector(2);
			p[i][j].K_gauss = CVector(2);
		}
	}
}
