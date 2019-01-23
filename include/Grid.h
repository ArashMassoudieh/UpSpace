#pragma once
#include <vector>
#include "Matrix.h"
#include "Vector.h"
#include "Matrix_arma.h"
#include "Matrix_arma_sp.h"
#include "BTCSet.h"
#include "Distribution.h"
#include "NormalDist.h"
#include "vtk.h"
#include "PathwaySet.h"
#include "_command.h"
#ifdef Qt_version
#include "qtextbrowser.h"
#endif // Qt_version


#define arma
//#define Qt_version

class HETEVAL;

using namespace std;

struct point
{
	double x;
	double y;
	CVector v;
	double t;
};

struct prop
{
	CVector K;
	CVector K_gauss;
	CVector V;
	bool k_det = false;
	vector<double> C;
};

struct prop_cell
{
	double H;
};

struct ijval
{
	int i;
	int j;
	double val;
};

struct correl_mat_vec
{
#ifdef  arma

	CMatrix_arma M_22;
	CVector_arma V_21;
	CVector_arma V_RHS;
#else
	CMatrix M_22;
	CVector V_21;
	CVector V_RHS;
#endif //  arma
};

struct _field_gen_params
{
	int max_correl_n;
	double k_correlation_lenght_scale_x;
	double k_correlation_lenght_scale_y;
};

struct _trajectory_params
{
	int T_n_samples;
	double T_dt;
	double T_end;
	double x0_trajs;
};

struct _filenames
{
	string K_output;
	string vx_output;
	string vy_output;
	string H_output;
	string T_output;
	string vx_dist_out;
	string vy_dist_out;
	string v_sect_dist_out;
	string K_input;
};

struct _geom_params
{
	int nx, ny;
	double dx, dy;
	int nx_data, ny_data;
};

struct _OU_params
{
	double kappa;
	CVector FinvU;
	CMatrix_arma Inv_M;
	CBTCSet BTCs;
	CBTCSet BTC_normal;
	CBTCSet BTCs_fw;
	CBTCSet BTC_normal_fw;
};




class CGrid
{
public:
	CGrid();
	~CGrid();
	_geom_params GP;
	int n_k_dets = 0;
	CPathwaySet Traj;
	_field_gen_params field_gen;
	_filenames filenames;
	_trajectory_params trajectory_params;
	CMatrix H;
	CMatrix vx;
	CMatrix vy;
	vector<ijval> get_closest_K_dets(int i, int j, int n);
	correl_mat_vec get_correll_matrix_vec(int i, int j);
	void assign_K_gauss(int i, int j);
	void assign_K_gauss();
	vector<vector<prop> > p;
	bool getfromfile(string filename, int _nx, int _ny);
	bool getKfromfile(string filename, int _nx, int _ny, int _nx_act, int _ny_act);
	bool createuniform(double K, int _nx, int _ny);
	CGrid(string filename);
	void creategrid(int nx, int ny, double dx, double dy);
	vector<CBTCSet> trajectories;
	CGrid(string filename_V, string filename_K, int _nx, int _ny);
	void writeasmatrix(string filename, int);
	void writeasmatrixK(string filename, int component);
	CVector getvelocity(point pp);
	CPathway gettrajectory(point pp, double dt, double t_end);
	CPathway gettrajectory_fix_dx(point pp, double dt, double t_end);
	CPathwaySet gettrajectories(double dt, double t_end);
	CPathwaySet gettrajectories_fixed_dx(double dt, double x_end);
	vector<point> pts;
        bool weighted;
	CBTC initialize(int numpoints, double x_0, double smoothing_factor=0, bool weighted=false);
	CMatrix_arma_sp create_stiffness_matrix_arma();
	CVector_arma create_RHS_arma();
	CVector create_RHS();
	CMatrix create_stiffness_matrix();
	int get_cell_no(int i, int j);
	int get_cell_no_OU(int i, int j);
	double leftboundary_h, rightboundary_h;
	double dt=1;
	double leftboundary_C;
	double D;
	CMatrix solve();
	vector<int> get_ij(int k);
	CBTC get_v_btc(int k);
	CBTC get_kg_btc(int k);
	void remap_K(int k);
	CBTC get_v_btc(double x,int k=0);
	CBTC get_v_dist(double x, int k=0, int nbins=40, double smoothing_factor=0);
    CBTC get_v_mag_btc();
	vector<_command> commands;
	void runcommands();
	CBTC vx_dist;
	CBTC vy_dist;
        CBTC v_dist;
	CBTCSet sect_dist;
	vector<double> marginal_K_dist_params;
	void set_inv_K_dist(int ninc);
	double map_to_KCDF(double u);
	double K_CDF(double x);
	double inv_function(double u, double guess = 0.5);
	double inv_function_s(double u);
	CBTC inv_K_dist;
	string marginal_K_dist_type;
	CBTC get_K_CDF(double x0, double x1, double log_inc);
	CBTC get_margina_traj_v_dist(double vmin, double vmax, double nbins, string val);
	double interpolate_K(double x, double y);
	vector<double> interpolate_V(double x, double y);
	void runcommands_qt();
	double max_K();
	double min_K();
	double max_vx();
	double min_vx();
	double max_vy();
	double min_vy();
	void set_K_transport(double dt, double D=0, double weight=0.5);
	CVector_arma create_RHS_transport(double dt, double weight=0.5, double D=0);
	CVector_arma create_RHS_transport_laplace(double weight, double D, double s);
	void solve_transport(double t_end);
	void solve_transport_laplace(double s);
	void set_K_transport_laplace(double D, double s);
	void create_f_inv_u();
	CVector_arma create_RHS_OU(double dt);
	void solve_transport_OU(double t_end);
	double time_weight;
	double min_v_x = 0;
	double max_v_x=0;
	void create_inverse_K_OU(double dt);
	void write_K_solution_to_vtp(string filename, double z_factor, bool _log);
	void write_C_to_vtp(string filename, double z_factor, bool _log, vector<double> t);
	void write_C_to_vtp(string filename, double z_factor, bool _log, double t);
	void clear();
	void showthings(vector<vtkSmartPointer<vtkActor>> actors, string filename = "");
	void write_K_field_to_vtp(string filename="surface.vtp", double z_factor=0.5, bool _log = false);
	vtkSmartPointer<vtkActor> traj_vtk_pdt(int trajno, double z_factor=0.5,double offset = 0);
	vtkSmartPointer<vtkPolyData> traj_vtk_pdt_vtp(int trajno, double z_factor=0.5, double offset=0, bool _log = false, bool _color = true);
	vector<vtkSmartPointer<vtkActor>> trajs_vtk_pdt(double z_factor=0.5, double offset = 0);
	void trajs_vtk_pdt_to_vtp(string filename = "paths.vtp", double z_factor = 0.5, double offset = 0, bool _log = false, bool _color = true);
    vtkSmartPointer<vtkActor> get_K_field_vtk_pdt(double z_factor=0.5);
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
    void show_K_field();
	void show_K_field_vtk(double z_factor = 0.5);

	void screen_shot(string filename = "screen_shot.png");
	void screenshot_test();
#ifdef Qt_version



	HETEVAL *main_window;
#endif // Qt_version
	CDistribution dist;
	CPathwaySet pset;
	double get_min_traj_v(int k);
	double get_max_traj_v(int k);
	CPathway dist_stores;
	string pathout;
	CMatrix_arma_sp Kv;
	CMatrix_arma_sp KD;
	CMatrix_arma_sp Kt;
	CMatrix C;
	_OU_params OU;
        void renormalize_k();
        void show_in_window(string s);
        void set_progress_value(double s);
        void clear_contents();
};


double avg(double x, double y, string type);
vector<ijval> get_top_n(vector<ijval> vec, int n);
vtkSmartPointer<vtkPolyData> pathway_vtk_pdt_vtp(CPathway &pathway, double z_factor=1, double offset=0);
void write_vtk(vtkSmartPointer<vtkPolyDataMapper> mapper, string filename);
vtkSmartPointer<vtkPolyDataMapper> pathways_vtk_pdt_vtp(CPathwaySet *pathwayset, double z_factor, double offset);
