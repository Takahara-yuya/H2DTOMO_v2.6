#pragma once
#ifndef TOMO_FSM2D_H_
#define TOMO_FSM2D_H_

#include "CommonHeaders.h"
#include "Tomo_mesh2D.h"
#include "Tomo_srs2D.h"
#include "tool.h"
#include "Filter.h"

class Tomofsm2D
{
public:
	Tomofsm2D() {
	}

	~Tomofsm2D() {
	}
    //mpi parameters//
    int np, myid;
    //
	int nsweep;

	string refrac_outroot, reflec_outroot, filetype;
	string rectt_outroot, rectt_reflec_draw_outroot, rectt_refrac_draw_outroot;

	double** grad_refl, ** grad_refr, grad_refr_max, grad_refl_max;

	bool T_refrac_print, T_reflec_print, grad_refr_print, grad_refl_print;

	bool precondition;

	double damp;

	void Init_fsm(string filename, TomoMesh2D model);

	void Fsm_fwd(TomoMesh2D model, TomoSrs2D data_obs, TomoSrs2D& data_cal_all, bool inv);

	void output(TomoSrs2D data_cal);

	void gradient_print(double** grad, int raycode, int n, TomoMesh2D model);

protected:

	void _Fsm_fwd_1shot_ttfield(int nx, int ny, int nrecfr, int nrecfl, int nx_s, int nx_e,
		double soux, double souy, vector< vector<int> >  sou_index, vector<int> sou_type, vector<int> intf_ynum,
		double** T1, double* tt_interface, vector<double> dX, vector<double> dY, vector< vector<double> >  Slowness,
		vector<double> xCen, vector<double> yCen, vector< vector<int> >  rec_refr_index, vector< vector<int> >  rec_refl_index,
		int raycode, int ishot);
	void _upa(int* a, int* b, int* c, int* d, int ix, int iz, int nx, int nz);
	double _fin_min(double Tx_min, double Tz_min, double dx, double dz, double s);
	void _add_source(int nx, int ny, double** T, double soux, double souy, int sx_index, int sy_index,
		int stype, vector<double> xCen, vector<double> yCen, vector< vector<double> >  Slowness);
	void _add_interface(double** T, double* tt_interface, int nx, vector<int> intf_ynum);
	void _ttfield_print(double** T, int ishot, int raycode, TomoMesh2D model);
	double _cal_rec_tt(double** T, int ishot, int irec, int x_index, int y_index, double recx, double recy, int type,
		vector<double>xCen, vector<double>yCen);
	//inverse
	void _get_grad_1shot(int nrec, int nx, int ny, int nx_s, int nx_e, vector<double> err1, vector<double> tobs, vector<double> tcal, vector<int> type1,
		vector< vector<int> > index, vector<double> recx1, vector<double> recy1, vector<double> xCen, vector<double> yCen, vector<int> intf_ynum,
		vector<double> dX, vector<double> dY, vector<double> xTopo, vector< vector<double> > Slowness,
		int ishot, double** T, double** T2, double** grad, double** dens, int raycode);
	void _calc_lamda(int nrec, int nx, int ny, int nx_s, int nx_e, vector<double>dx, vector<double>dy, vector<int> intf_ynum,
		double** T1, double** lamda, double* lamda_tmp, int** lamda_type_nidex, int ishot, int minmax, int raycode, int n2s);
	void _upab(double dx1, double dx2, double dy1, double dy2, double** T1, double& a_p_1, double& a_p_2, double& a_r_1,
		double& a_r_2, double& b_p_1, double& b_p_2, double& b_r_1, double& b_r_2, int ix, int iy);
	void _uplamda(double dx1, double dx2, double dy1, double dy2, double** lamda, double a_p_1, double a_p_2, double a_r_1, double a_r_2,
		double b_p_1, double b_p_2, double b_r_1, double b_r_2, int ix, int iy, int minmax);
	void _precondition(double** dens, double* rn, TomoMesh2D model, TomoSrs2D data_obs, int raycode, int ishot);
};

#endif
