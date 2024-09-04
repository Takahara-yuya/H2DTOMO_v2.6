/*
***********************************************************************

Tomo_inverse2D.cpp	(Tomo inversion)
This file is part of H2DTOMO.

***********************************************************************

Nov 22, 2023
Copyright 2023

Zuwei Huang
hzw1498218560@tongji.edu.cn
School of Ocean and Earth Science, Tongji University
Integrated Geophysics Group

Chongjin Zhao
zcjadc@126.com
School of Ocean and Earth Science, Tongji University
Integrated Geophysics Group

version 2.6.0


***********************************************************************
*/

#include "../include/Tomo_inverse2D.h"

void Tomoinv2D::Init_inv(string filename)
{
	string tmp, tmp2, tmp3;
	ifstream fin(filename);
	fin >> tmp >> tmp >> tmp >> Max_iteration;
	fin >> tmp >> tmp >> dv_max;
	dv_max *= 1000.;
	fin >> tmp >> tmp >> tmp3;
	for (char& c : tmp3) {
		c = std::tolower(c);
	}
	if (tmp3 == "yes")
	{
		outtmp = true;
	}
	else
	{
		outtmp = false;
	}
	fin >> tmp >> tmp >> target_misfit;
	fin >> tmp >> tmp >> vmin >> vmax;
	vmin *= 1000.;
	vmax *= 1000.;
	fin >> tmp >> tmp >> sm_type;
	for (char& c : sm_type) {
		c = std::tolower(c);
	}
	fin >> tmp >> tmp >> smx >> smy >> smtimes;
	fin >> tmp >> tmp >> fracx >> fracy >> vref >> frequency;
	vref *= 1000.;
	fin >> tmp >> tmp >> steplen;
	fin >> tmp >> filetype;
	for (char& c : filetype) {
		c = std::tolower(c);
	}
	fin >> tmp >> tmp >> LogFile_Root;
	fin >> tmp >> tmp >> outroot;
	fin >> tmp >> tmp >> tmp >> tmp_outroot;
	fin.close();
}

void Tomoinv2D::ASTinv(TomoMesh2D& model, TomoMesh2D mapr, TomoSrs2D data_obs, TomoSrs2D& data_cal, Tomofsm2D& fwd)
{
	///
	///initial model min-max
	for (int ix = 0; ix < model.Num_xCen; ix++)
	{
		for (int iy = 0; iy < model.Num_yCen; iy++)
		{
			int index = iy + ix * model.Num_yCen;
			if (model.Vel1D.coeffRef(index) < vmin && iy > model.topo_ynum[ix] - 1)model.Vel1D.coeffRef(index) = vmin;
			if (model.Vel1D.coeffRef(index) > vmax && iy > model.topo_ynum[ix] - 1)model.Vel1D.coeffRef(index) = vmax;
			model.Vel[ix][iy] = model.Vel1D.coeffRef(index);
			model.Slowness[ix][iy] = 1. / model.Vel[ix][iy];
			model.Slowness1D.coeffRef(index) = model.Slowness[ix][iy];
		}
	}
	///
	ofstream flog;
	if (myid == 0)
	{
		model.Outputmesh(outroot + "_init.dat", model.filetype);
		flog.open(LogFile_Root);
		if (!flog)
		{
			cerr << "Error writing LogFile" << endl;
			exit(1);
		}
		std::cerr << endl;
		std::cerr << "-------------INVERSION MODULE INITIALIZATING!--------------" << endl;
		std::cerr << endl;
		flog << "niter data_misfit refr_misfit  dt_total_refr(Refrac) refl_misfit dt_total_refl(Reflec) max_vel_perb roughness" << endl;
	}
	//_Stabilizer_Function(model);
	int nit, nx, ny;
	int nsearch;
	TomoMesh2D model_test;
	TomoSrs2D data_cal_test;
	data_cal_test = data_cal;
	nx = model.Num_xCen; ny = model.Num_yCen;
	VectorXd gl, gl_last, pl, pl_last;
	gl.resize(nx * ny); gl_last.resize(nx * ny); pl.resize(nx * ny);
	double data_misfit, data_misfit_last, refl_misfit, refr_misfit, dt_total_refr, dt_total_refl, ke, ks, kn, grad_max, frmax, flmax;
	double data_misfit_test, refl_misfit_test, refr_misfit_test, dt_total_refr_test, dt_total_refl_test;
	double dt_total_0;
	double R;
	double* v_ref;
	data_misfit = 99999.;
	double* rn_refr, * rn_refl;
	double* rn_refr_test, * rn_refl_test;
	rn_refr = alloc_double_1d(data_obs.Num_rec_refr);
	rn_refl = alloc_double_1d(data_obs.Num_rec_refl);
	rn_refr_test = alloc_double_1d(data_obs.Num_rec_refr);
	rn_refl_test = alloc_double_1d(data_obs.Num_rec_refl);
	v_ref = alloc_double_1d(ny);
	zero_1d(v_ref, ny);
	bool new_iter, conj_break;
	new_iter = true; nit = 0; data_misfit_last = 1.e+9; conj_break = false;
	while (new_iter)
	{
		nit++;
		fwd.Fsm_fwd(model, data_obs, data_cal, true);
		////smooth the grad
		if (sm_type == "gauss")
		{
			smoothgauss2dtopo(fwd.grad_refr, ny, nx, model.dY[0], frequency, model.topo_ynum, 2000, fracx, fracy, vref, 1);
			smoothgauss2dtopo(fwd.grad_refl, ny, nx, model.dY[0], frequency, model.topo_ynum, 2000, fracx, fracy, vref, 1);
		}
		else
		{
			movingAverage2D(fwd.grad_refl, nx, ny, smx, smy, smtimes);
			movingAverage2D(fwd.grad_refr, nx, ny, smx, smy, smtimes);
		}
		//add constraint
		frmax = -99999.; flmax = -99999.;
		for (int ix = 0; ix < nx; ix++)
			for (int iy = 0; iy < ny; iy++)
			{
				if (iy < model.topo_ynum[ix])
				{
					fwd.grad_refr[ix][iy] = 0;
					fwd.grad_refl[ix][iy] = 0;
				}
				if (iy > model.intf_ynum[ix] - 1)
				{
					fwd.grad_refr[ix][iy] = 0;
					fwd.grad_refl[ix][iy] = 0;
				}
				if (abs(fwd.grad_refr[ix][iy]) > frmax)frmax = abs(fwd.grad_refr[ix][iy]);
				if (abs(fwd.grad_refl[ix][iy]) > flmax)flmax = abs(fwd.grad_refl[ix][iy]);
			}
		//min-max normoliazation
		normolization(fwd.grad_refr, nx, ny);
		normolization(fwd.grad_refl, nx, ny);
		//output
		fwd.output(data_cal);
		if (fwd.grad_refr_print && outtmp)fwd.gradient_print(fwd.grad_refr, 0, nit, model);
		if (fwd.grad_refl_print && outtmp)fwd.gradient_print(fwd.grad_refl, 1, nit, model);
		////calculate roughness
		for (int iy = 0; iy < ny; iy++)
		{
			for (int ix = 0; ix < nx; ix++)
			{
				if (frmax > 0)fwd.grad_refr[ix][iy] = fwd.grad_refr[ix][iy] / frmax;
				else fwd.grad_refr[ix][iy] = 0;
				if (flmax > 0)fwd.grad_refl[ix][iy] = fwd.grad_refl[ix][iy] / flmax;
				else fwd.grad_refl[ix][iy] = 0;
				v_ref[iy] += model.Vel[ix][iy];
			}
			v_ref[iy] /= double(ny);
		}
		for (int iy = 0; iy < ny; iy++)
		{
			for (int ix = 0; ix < nx; ix++)
			{
				R += pow((model.Vel[ix][iy] - v_ref[iy]) / v_ref[iy], 2);
			}
		}
		R = sqrt(R / (nx * ny)) * 100.;
		////misfit
		_misfit_calculate(rn_refr, rn_refl, data_obs, data_cal, data_misfit, refr_misfit, refl_misfit, dt_total_refr, dt_total_refl);
		if (nit == 1)dt_total_0 = dt_total_refr + dt_total_refl;
		if (myid == 0)
		{
			cerr << endl;
			cerr << "********************************INVERSION********************************" << endl;
			cerr << "NUMBER OF ITERATION: " << nit << endl;
			cerr << "Data Misfit: " << data_misfit << "	Refraction Misfit: " << refr_misfit << " " << "	Reflection Misfit: " << refl_misfit << endl;
			cerr << "Refraction traveltime residual: " << dt_total_refr << "		" << "Reflection traveltime residual: " << dt_total_refl << endl;
			cerr << "Roughness " << R << endl;
			cerr << "***************************************************************************" << endl;
			cerr << endl;
		}
		if (data_misfit_last - data_misfit < 0.2)conj_break = true;
		for (int ix = 0; ix < nx; ix++)
			for (int iy = 0; iy < ny; iy++)
			{
				gl(iy + ix * ny) = dt_total_refr / (dt_total_refr + dt_total_refl) * fwd.grad_refr[ix][iy] + 0.3 * dt_total_refl / (dt_total_refr + dt_total_refl) * fwd.grad_refl[ix][iy];
			}
		if (nit == 1 || conj_break)
		{
			pl = gl;
			conj_break = false;
		}
		else
		{
			if (gl_last.norm() != 0)
				pl = gl + (gl.norm() * gl.norm() / (gl_last.norm() * gl_last.norm())) * pl_last;
			else
				pl.setZero();
		}
		grad_max = 0;
		for (int ix = 0; ix < nx; ix++)
			for (int iy = 0; iy < ny; iy++)
			{
				if (abs(pl(iy + ix * ny)) > grad_max)grad_max = abs(pl(iy + ix * ny));
			}
		gl_last = gl;
		pl_last = pl;
		if (_Iteration_Termination_Condition(nit, data_misfit))break;
		if (grad_max == 0)break;
		if (steplen == 1)
		{
			//golden research!
			if ((dv_max - 100 * double(nit)) > 200.)
			{
				ke = (dv_max - 100 * double(nit)) / grad_max;
			}
			else
			{
				ke = 200. / grad_max;
			}
			ks = 0;
			_Golden_section_Linesearch(model, fwd, pl, data_obs, data_cal, kn, ks, ke, 1. / grad_max);
		}
		else
		{
			if (steplen == 2)
			{
				kn = 1. / grad_max * dv_max * sqrt((dt_total_refr + dt_total_refl) / dt_total_0);
			}
			///test///
			nsearch = 0;
			while (true)
			{
				model_test = model;
				nsearch++;
				if (nsearch > 10)
				{
					kn = 0;
					smtimes--;
					if (smtimes < 1)
					{
						new_iter = false;
						if (myid == 0)
						{
							cerr << endl;
							cerr << "Convergence problem, Iteration terminate!" << endl;
						}
					}
					else
					{
						if (myid == 0)
						{
							cerr << endl;
							cerr << "Convergence problem, decrease smooth times!" << endl;
						}
					}
					break;
				}
				_Mesh_update(model_test, kn, pl);
				fwd.Fsm_fwd(model_test, data_obs, data_cal_test, false);
				_misfit_calculate(rn_refr_test, rn_refl_test, data_obs, data_cal, data_misfit_test, refr_misfit_test, refl_misfit_test, dt_total_refr_test, dt_total_refl_test);
				if (data_misfit_test < data_misfit)break;
				kn *= 0.8;
			}
			if (new_iter)
			{
				if (myid == 0)
				{
					cerr << "Misfit decrease, iterate continue!" << endl;
					cerr << endl;
				}
			}
		}
		_Mesh_update(model, kn, pl);
		if (myid == 0)
		{
			flog << nit << " " << data_misfit << " " << refr_misfit << " " << dt_total_refr << " " << refl_misfit << " " << dt_total_refl << " " << kn * grad_max << " " << R << endl;
			if (outtmp)model.Outputmesh(tmp_outroot + to_string(nit) + ".dat", filetype);
			model.Outputmesh(outroot, filetype);
		}
		data_misfit_last = data_misfit;
	}
	flog.close();
}

bool Tomoinv2D::_Iteration_Termination_Condition(int nit, double misfit)
{
	if (nit > Max_iteration)
	{
		if (myid == 0)
		{
			cerr << endl;
			cerr << "Maximum number of iterations reached, iteration terminated!" << endl;
		}
		return true;
	}
	else if (misfit < target_misfit)
	{
		if (myid == 0)
		{
			cerr << endl;
			cerr << "Target misfit achieved, iteration terminated!" << endl;
		}
		return true;
	}
	/*else if (alpha < 0.00001)
	{
		cerr << "Regularization factor less than 0.00001, iteration terminated!" << endl;
		return true;
	}*/
	else
	{
		return false;
	}
}

void Tomoinv2D::_misfit_calculate(double* rn_refr, double* rn_refl, TomoSrs2D data_obs, TomoSrs2D data_cal, double& data_misfit,
	double& refr_misfit, double& refl_misfit, double& dt_refr_total, double& dt_refl_total)
{
	int count1, count2;
	count1 = 0; count2 = 0;
	int nrec_refl, nrec_refr;
	double t_cal, t_obs;
	double nshot = data_obs.Num_source;
	refr_misfit = 0; refl_misfit = 0; dt_refr_total = 0; dt_refl_total = 0;
	for (int ishot = 0; ishot < nshot; ishot++)
	{
		nrec_refr = data_obs.data[ishot].nrec_refrac;
		nrec_refl = data_obs.data[ishot].nrec_reflec;
		for (int irec = 0; irec < nrec_refr; irec++)
		{
			t_cal = data_cal.data[ishot].time_refrac[irec]; t_obs = data_obs.data[ishot].time_refrac[irec];
			rn_refr[count1] = t_cal - t_obs;
			refr_misfit += pow((t_cal - t_obs) / data_cal.data[ishot].err_refrac[irec], 2);
			dt_refr_total += pow(rn_refr[count1], 2);
			count1++;
		}
		for (int irec = 0; irec < nrec_refl; irec++)
		{
			t_cal = data_cal.data[ishot].time_reflec[irec]; t_obs = data_obs.data[ishot].time_reflec[irec];
			rn_refl[count2] = t_cal - t_obs;
			refl_misfit += pow((t_cal - t_obs) / data_cal.data[ishot].err_reflec[irec], 2);
			dt_refl_total += pow(rn_refl[count2], 2);
			count2++;
		}
	}
	if (data_obs.Num_rec_refr > 0)
		refr_misfit = sqrt(refr_misfit / data_obs.Num_rec_refr);
	else
		refr_misfit = 0;
	if (data_obs.Num_rec_refl > 0)
		refl_misfit = sqrt(refl_misfit / data_obs.Num_rec_refl);
	else
		refl_misfit = 0;
	data_misfit = double(data_obs.Num_rec_refr) / double(data_obs.Num_rec_refl + data_obs.Num_rec_refr) * refr_misfit
		+ double(data_obs.Num_rec_refl) / double(data_obs.Num_rec_refl + data_obs.Num_rec_refr) * refl_misfit;
}

void Tomoinv2D::_Golden_section_Linesearch(TomoMesh2D model, Tomofsm2D fwd, VectorXd gl, TomoSrs2D data_obs, TomoSrs2D data_cal,
	double& kn, double ks, double ke, double vref)
{
	if (myid == 0)cerr << "Golden section linesearch start!" << endl;
	double data_misfit, refl_misfit, refr_misfit, dt_total_refr, dt_total_refl;
	double* rn_refr_test, * rn_refl_test;
	rn_refr_test = alloc_double_1d(data_obs.Num_rec_refr);
	rn_refl_test = alloc_double_1d(data_obs.Num_rec_refl);
	TomoMesh2D model_test;
	VectorXd t, s, upper, lower, theta_t, theta_s;
	t.resize(50); s.resize(50);
	upper.resize(50); lower.resize(50);
	theta_t.resize(50); theta_s.resize(50);
	///
	lower(0) = ks; upper(0) = ke;//log
	t(0) = lower(0) + 0.382 * (upper(0) - lower(0));
	s(0) = lower(0) + 0.618 * (upper(0) - lower(0));
	//p0//Left
	model_test = model;
	_Mesh_update(model_test, t(0), gl);
	fwd.Fsm_fwd(model_test, data_obs, data_cal, false);
	_misfit_calculate(rn_refr_test, rn_refl_test, data_obs, data_cal, theta_t(0), refr_misfit, refl_misfit, dt_total_refr, dt_total_refl);
	if (myid == 0)
	{
		cerr << "Data Misfit: " << theta_t(0) << " " << "at kn=: " << t(0) / vref << endl;
		cerr << endl;
	}
	//q0//Right
	model_test = model;
	_Mesh_update(model_test, s(0), gl);
	fwd.Fsm_fwd(model_test, data_obs, data_cal, false);
	_misfit_calculate(rn_refr_test, rn_refl_test, data_obs, data_cal, theta_s(0), refr_misfit, refl_misfit, dt_total_refr, dt_total_refl);
	if (myid == 0)
	{
		cerr << "Data Misfit: " << theta_s(0) << " " << "at kn=: " << s(0) / vref << endl;
		cerr << endl;
	}
	//
	int kk = 0;
	while (1)
	{
		if (theta_t(kk) >= theta_s(kk))//left
		{
			if ((upper(kk) - t(kk)) < 40 * vref)
			{
				free_1d(rn_refr_test);
				free_1d(rn_refl_test);
				kn = s(kk);
				break;
			}
			else
			{
				lower(kk + 1) = t(kk);
				upper(kk + 1) = upper(kk);
				t(kk + 1) = s(kk);
				s(kk + 1) = lower(kk + 1) + 0.618 * (upper(kk + 1) - lower(kk + 1));
				//p0//Left
				model_test = model;
				_Mesh_update(model_test, t(kk + 1), gl);
				fwd.Fsm_fwd(model_test, data_obs, data_cal, false);
				_misfit_calculate(rn_refr_test, rn_refl_test, data_obs, data_cal, theta_t(kk + 1), refr_misfit, refl_misfit, dt_total_refr, dt_total_refl);
				if (myid == 0)
				{
					cerr << "Data Misfit: " << theta_t(kk + 1) << " " << "at kn=: " << t(kk + 1) / vref << endl;
					cerr << endl;
				}
				//q0//Right
				model_test = model;
				_Mesh_update(model_test, s(kk + 1), gl);
				fwd.Fsm_fwd(model_test, data_obs, data_cal, false);
				_misfit_calculate(rn_refr_test, rn_refl_test, data_obs, data_cal, theta_s(kk + 1), refr_misfit, refl_misfit, dt_total_refr, dt_total_refl);
				if (myid == 0)
				{
					cerr << "Data Misfit: " << theta_s(kk + 1) << " " << "at kn=: " << s(kk + 1) / vref << endl;
					cerr << endl;
				}
				kk++;
			}
		}
		else
		{
			if (s(kk) - lower(kk) < 40 * vref)
			{
				kn = t(kk);
				free_1d(rn_refr_test);
				free_1d(rn_refl_test);
				break;
			}
			else
			{
				lower(kk + 1) = lower(kk);
				upper(kk + 1) = s(kk);
				s(kk + 1) = t(kk);
				t(kk + 1) = lower(kk + 1) + 0.382 * (upper(kk + 1) - lower(kk + 1));
				//p0//Left
				model_test = model;
				_Mesh_update(model_test, t(kk + 1), gl);
				fwd.Fsm_fwd(model_test, data_obs, data_cal, false);
				_misfit_calculate(rn_refr_test, rn_refl_test, data_obs, data_cal, theta_t(kk + 1), refr_misfit, refl_misfit, dt_total_refr, dt_total_refl);
				if (myid == 0)
				{
					cerr << "Data Misfit: " << theta_t(kk + 1) << " " << "at kn=: " << t(kk + 1) / vref << endl;
					cerr << endl;
				}
				//q0//Right
				model_test = model;
				_Mesh_update(model_test, s(kk + 1), gl);
				fwd.Fsm_fwd(model_test, data_obs, data_cal, false);
				_misfit_calculate(rn_refr_test, rn_refl_test, data_obs, data_cal, theta_s(kk + 1), refr_misfit, refl_misfit, dt_total_refr, dt_total_refl);
				if (myid == 0)
				{
					cerr << "Data Misfit: " << theta_s(kk + 1) << " " << "at kn=: " << s(kk + 1) / vref << endl;
					cerr << endl;
				}
				kk++;
			}
		}
		if (kk > 20)
		{
			free_1d(rn_refr_test);
			free_1d(rn_refl_test);
			kn = s(kk);
			if (myid == 0)cerr << "Golden Research Failed!" << endl;
			break;
		}
	}
}

void Tomoinv2D::_Mesh_update(TomoMesh2D& model, double kn, VectorXd gl)
{
	model.Vel1D += kn * gl;
	for (int ix = 0; ix < model.Num_xCen; ix++)
	{
		for (int iy = 0; iy < model.Num_yCen; iy++)
		{
			int index = iy + ix * model.Num_yCen;
			if (model.Vel1D.coeffRef(index) < vmin && iy > model.topo_ynum[ix] - 1)model.Vel1D.coeffRef(index) = vmin;
			if (model.Vel1D.coeffRef(index) > vmax && iy > model.topo_ynum[ix] - 1)model.Vel1D.coeffRef(index) = vmax;
			model.Vel[ix][iy] = model.Vel1D.coeffRef(index);
			model.Slowness[ix][iy] = 1. / model.Vel[ix][iy];
			model.Slowness1D.coeffRef(index) = model.Slowness[ix][iy];
		}
	}
}

void Tomoinv2D::_Stabilizer_Function(TomoMesh2D model)
{
	vector<Triplet<double> > triplets, triplets2, triplets3;
	//vertical smooth
	for (int i = 0; i < model.Num_xCen; i++)
	{
		for (int j = 1; j < model.Num_yCen; j++)
		{
			int np = j + i * model.Num_yCen;
			triplets.push_back(Triplet<double>(np, np, 1));//EVENLY WEIGHT
			triplets.push_back(Triplet<double>(np, np - 1, -1));//EVENLY WEIGHT
		}
	}
	Wey.resize(model.Num_Cen, model.Num_Cen);
	Wey.setFromTriplets(triplets.begin(), triplets.end());
	Wey.makeCompressed();
	//horizontal smooth
	for (int np = 0; np < model.Num_xCen - 1; np++)//L layer
	{
		for (int nl = 0; nl < model.Num_yCen; nl++)//X
		{
			int tmp = nl + np * model.Num_yCen;
			triplets2.push_back(Triplet<double>(tmp, tmp, -1));
			triplets2.push_back(Triplet<double>(tmp, tmp + model.Num_yCen, 1));
		}
	}
	Wex.resize(model.Num_Cen, model.Num_Cen);
	Wex.setFromTriplets(triplets2.begin(), triplets2.end());
	Wex.makeCompressed();
	//
	for (int np = 0; np < model.Num_xCen; np++)//L layer
	{
		for (int nl = 0; nl < model.Num_yCen; nl++)//X
		{
			int tmp = nl + np * model.Num_yCen;
			triplets3.push_back(Triplet<double>(tmp, tmp, 1));
		}
	}
	We.resize(model.Num_Cen, model.Num_Cen);
	We.setFromTriplets(triplets3.begin(), triplets3.end());
	We.makeCompressed();
}
