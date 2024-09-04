#pragma once
#ifndef TOMO_INVERSE2D_H_
#define TOMO_INVERSE2D_H_
#include "CommonHeaders.h"
#include "Tomo_fsm2D.h"

class Tomoinv2D
{
public:
	Tomoinv2D() {
	}

	~Tomoinv2D() {
	}
    //mpi parameter//
    int np, myid;
    //
	void Init_inv(string filename);

	void ASTinv(TomoMesh2D& model, TomoMesh2D mapr, TomoSrs2D data_obs, TomoSrs2D& data_cal, Tomofsm2D& fwd);
	double vmin, vmax;
protected:

	string LogFile_Root;
	string outroot, tmp_outroot, filetype;
	int steplen;
	SparseMatrix<double> We, Wex, Wey;

	bool outtmp;
	int Max_iteration;
	int smx, smy, smtimes;
	double fracx, fracy, vref, frequency;
	double dv_max;
	double alpha, target_misfit;
	string sm_type;

	bool _Iteration_Termination_Condition(int nit, double misfit);
	void _misfit_calculate(double* rn_refr, double* rn_refl, TomoSrs2D data_obs, TomoSrs2D data_cal, double& data_misfit,
		double& refr_misfit, double& refl_misfit, double& dt_refr_total, double& dt_refl_total);
	void _Golden_section_Linesearch(TomoMesh2D model, Tomofsm2D fwd, VectorXd gl, TomoSrs2D data_obs, TomoSrs2D data_cal,
		double& kn, double ks, double ke, double vref);
	void _Mesh_update(TomoMesh2D& model, double kn, VectorXd gl);
	void _Stabilizer_Function(TomoMesh2D model);
};

#endif
