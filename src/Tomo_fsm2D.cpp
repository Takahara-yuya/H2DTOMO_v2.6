/*
***********************************************************************

Tomo_fsm2D.cpp	(Forward modeling)
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

#include "../include/Tomo_fsm2D.h"

void Tomofsm2D::Init_fsm(string filename, TomoMesh2D model)
{
	string tmp, tmp2, tmp3;
	ifstream fin(filename);
	fin >> tmp >> tmp >> tmp >> nsweep;
	fin >> tmp >> tmp >> tmp >> tmp2;
	fin >> tmp >> tmp >> tmp >> tmp3;
	for (char& c : tmp2) {
		c = std::tolower(c);
	}
	for (char& c : tmp3) {
		c = std::tolower(c);
	}
	if (tmp2 == "yes")T_refrac_print = true;
	else T_refrac_print = false;
	if (tmp3 == "yes")T_reflec_print = true;
	else T_reflec_print = false;
	fin >> tmp >> tmp >> tmp2;
	for (char& c : tmp2) {
		c = std::tolower(c);
	}
	if (tmp2 == "yes")precondition = true;
	else precondition = false;
	fin >> tmp >> tmp >> damp;
	fin >> tmp >> tmp >> tmp >> tmp2;
	fin >> tmp >> tmp >> tmp >> tmp3;
	for (char& c : tmp2) {
		c = std::tolower(c);
	}
	for (char& c : tmp3) {
		c = std::tolower(c);
	}
	if (tmp2 == "yes")grad_refr_print = true;
	else grad_refr_print = false;
	if (tmp3 == "yes")grad_refl_print = true;
	else grad_refl_print = false;
	fin >> tmp >> tmp >> refrac_outroot;
	fin >> tmp >> tmp >> reflec_outroot;
	fin >> tmp >> filetype;
	fin >> tmp >> tmp >> rectt_outroot;
	fin >> tmp >> tmp >> tmp >> rectt_refrac_draw_outroot;
	fin >> tmp >> tmp >> tmp >> rectt_reflec_draw_outroot;
	fin.close();
	grad_refl = alloc_double_2d(model.Num_xCen, model.Num_yCen);
	grad_refr = alloc_double_2d(model.Num_xCen, model.Num_yCen);
}

void Tomofsm2D::output(TomoSrs2D data_cal)
{
	ofstream fout(rectt_outroot + ".dat");
	ofstream fout1(rectt_refrac_draw_outroot + ".dat");
	ofstream fout2(rectt_reflec_draw_outroot + ".dat");
	if (!fout || !fout1 || !fout2)
	{
		cerr << "Error writing travel time file" << endl;
		exit(1);
	}
	fout << data_cal.Num_source << endl;
	for (int is = 0; is < data_cal.Num_source; is++)
	{
		fout << "s" << " " << data_cal.soux[is] / 1000. << " " << data_cal.souy[is] / 1000. << " "
			<< (data_cal.data[is].nrec_reflec + data_cal.data[is].nrec_refrac) << endl;
		for (int ifrac = 0; ifrac < data_cal.data[is].nrec_refrac; ifrac++)
		{
			fout << "r" << " " << data_cal.data[is].recx_refrac[ifrac] / 1000. << " " << data_cal.data[is].recy_refrac[ifrac] / 1000. << " "
				<< "0" << " " << data_cal.data[is].time_refrac[ifrac] << " " << data_cal.data[is].err_refrac[ifrac] << endl;
			fout1 << data_cal.data[is].recx_refrac[ifrac] / 1000. << " " << data_cal.data[is].recy_refrac[ifrac] / 1000. << " "
				<< data_cal.data[is].time_refrac[ifrac] << endl;
		}
		for (int iflec = 0; iflec < data_cal.data[is].nrec_reflec; iflec++)
		{
			fout << "r" << " " << data_cal.data[is].recx_reflec[iflec] / 1000. << " " << data_cal.data[is].recy_reflec[iflec] / 1000. << " "
				<< "1" << " " << data_cal.data[is].time_reflec[iflec] << " " << data_cal.data[is].err_reflec[iflec] << endl;
			fout2 << data_cal.data[is].recx_reflec[iflec] / 1000. << " " << data_cal.data[is].recy_reflec[iflec] / 1000. << " "
				<< data_cal.data[is].time_reflec[iflec] << endl;
		}
	}
	fout.close();
	fout1.close();
	fout2.close();
}

void Tomofsm2D::Fsm_fwd(TomoMesh2D model, TomoSrs2D data_obs, TomoSrs2D& data_cal, bool inv)
{
	for (int ishot = 0; ishot < data_cal.Num_source; ishot++)
	{
		for (int irec = 0; irec < data_cal.data[ishot].nrec_refrac; irec++)data_cal.data[ishot].time_refrac[irec] = 0;
		for (int irec = 0; irec < data_cal.data[ishot].nrec_reflec; irec++)data_cal.data[ishot].time_reflec[irec] = 0;
	}
	auto start_time = std::chrono::high_resolution_clock::now();
	//
	//ishot
	int nx, ny;
	nx = model.Num_xCen;
	ny = model.Num_yCen;
	if (inv)
	{
		zero_2d(grad_refl, nx, ny);
		zero_2d(grad_refr, nx, ny);
	}
	int nx_s, nx_e;
	double** grad_refl_tmp, ** grad_refr_tmp, ** dens_refr, ** dens_refl, ** dens_refr_tmp, ** dens_refl_tmp;
	double** T_refrac, ** T_reflec, * tt_interface;
	grad_refl_tmp = alloc_double_2d(nx, ny);
	grad_refr_tmp = alloc_double_2d(nx, ny);
	dens_refr_tmp = alloc_double_2d(nx, ny);
	dens_refl_tmp = alloc_double_2d(nx, ny);
	dens_refr = alloc_double_2d(nx, ny);
	dens_refl = alloc_double_2d(nx, ny);
	T_refrac = alloc_double_2d(nx, ny);
	T_reflec = alloc_double_2d(nx, ny);
	tt_interface = alloc_double_1d(nx);
	zero_2d(dens_refr, nx, ny);
	zero_2d(dens_refl, nx, ny);
	zero_2d(T_refrac, nx, ny);
	zero_2d(T_reflec, nx, ny);
	if (myid == 0)
	{
		cerr << "Fast sweeping forward modeling start!" << endl;
		cerr << "Shot number:" << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for (int ishot = myid; ishot < data_cal.Num_source; ishot += np)
	{
		cerr << ishot + 1 << " ";
		///block
		nx_s = 999999999;
		nx_e = -999999999;
		for (int i = 0; i < data_obs.data[ishot].nrec_refrac; i++)
		{
			if (nx_s > data_obs.data[ishot].rec_refr_index[0][i])nx_s = data_obs.data[ishot].rec_refr_index[0][i];
			if (nx_e < data_obs.data[ishot].rec_refr_index[0][i])nx_e = data_obs.data[ishot].rec_refr_index[0][i];
		}
		for (int i = 0; i < data_obs.data[ishot].nrec_reflec; i++)
		{
			if (nx_s > data_obs.data[ishot].rec_refl_index[0][i])nx_s = data_obs.data[ishot].rec_refl_index[0][i];
			if (nx_e < data_obs.data[ishot].rec_refl_index[0][i])nx_e = data_obs.data[ishot].rec_refl_index[0][i];
		}
		if (nx_s > data_obs.sou_index[0][ishot])nx_s = data_obs.sou_index[0][ishot];
		if (nx_e < data_obs.sou_index[0][ishot])nx_e = data_obs.sou_index[0][ishot];
		nx_s -= 10;
		if (nx_s < 0)nx_s = 0;
		nx_e += 10;
		if (nx_e > nx)nx_e = nx;
		/////
		if (data_cal.data[ishot].nrec_refrac > 0)
		{
			_Fsm_fwd_1shot_ttfield(nx, ny, data_obs.data[ishot].nrec_refrac, data_obs.data[ishot].nrec_reflec, nx_s, nx_e,
				data_obs.soux[ishot], data_obs.souy[ishot], data_obs.sou_index, data_obs.sou_type, model.intf_ynum, T_refrac, tt_interface,
				model.dX, model.dY, model.Slowness, model.xCen, model.yCen, data_obs.data[ishot].rec_refr_index, data_obs.data[ishot].rec_refl_index,
				0, ishot);
		}
		if (data_cal.data[ishot].nrec_reflec > 0)
		{
			if (data_cal.data[ishot].nrec_refrac == 0)
			{
				_Fsm_fwd_1shot_ttfield(nx, ny, data_obs.data[ishot].nrec_refrac, data_obs.data[ishot].nrec_reflec, nx_s, nx_e,
					data_obs.soux[ishot], data_obs.souy[ishot], data_obs.sou_index, data_obs.sou_type, model.intf_ynum, T_refrac, tt_interface,
					model.dX, model.dY, model.Slowness, model.xCen, model.yCen, data_obs.data[ishot].rec_refr_index, data_obs.data[ishot].rec_refl_index,
					0, ishot);
			}
			for (int ix = 0; ix < nx; ix++)
				tt_interface[ix] = T_refrac[ix][model.intf_ynum[ix] - 1];
			_Fsm_fwd_1shot_ttfield(nx, ny, data_obs.data[ishot].nrec_refrac, data_obs.data[ishot].nrec_reflec, nx_s, nx_e,
				data_obs.soux[ishot], data_obs.souy[ishot], data_obs.sou_index, data_obs.sou_type, model.intf_ynum, T_reflec, tt_interface,
				model.dX, model.dY, model.Slowness, model.xCen, model.yCen, data_obs.data[ishot].rec_refr_index, data_obs.data[ishot].rec_refl_index,
				1, ishot);
		}
		//get Traveltime for nrec
		for (int irec = 0; irec < data_cal.data[ishot].nrec_refrac; irec++)
		{
			data_cal.data[ishot].time_refrac[irec] = _cal_rec_tt(T_refrac, ishot, irec, data_cal.data[ishot].rec_refr_index[0][irec], data_cal.data[ishot].rec_refr_index[1][irec],
				data_cal.data[ishot].recx_refrac[irec], data_cal.data[ishot].recy_refrac[irec], data_cal.data[ishot].rec_refr_type[irec], model.xCen, model.yCen);
		}
		//
		for (int irec = 0; irec < data_cal.data[ishot].nrec_reflec; irec++)
		{
			data_cal.data[ishot].time_reflec[irec] = _cal_rec_tt(T_reflec, ishot, irec, data_cal.data[ishot].rec_refl_index[0][irec], data_cal.data[ishot].rec_refl_index[1][irec],
				data_cal.data[ishot].recx_reflec[irec], data_cal.data[ishot].recy_reflec[irec], data_cal.data[ishot].rec_refl_type[irec], model.xCen, model.yCen);
		}
		if (T_refrac_print)_ttfield_print(T_refrac, int(ishot + 1), 0, model);
		if (T_reflec_print)_ttfield_print(T_reflec, int(ishot + 1), 1, model);
		//inverse!!
		if (inv)
		{
			zero_2d(grad_refl_tmp, nx, ny);
			zero_2d(grad_refr_tmp, nx, ny);
			zero_2d(dens_refr_tmp, nx, ny);
			zero_2d(dens_refl_tmp, nx, ny);
			if (data_cal.data[ishot].nrec_refrac > 0)
			{
				_get_grad_1shot(data_cal.data[ishot].nrec_refrac, nx, ny, nx_s, nx_e, data_obs.data[ishot].err_refrac
					, data_obs.data[ishot].time_refrac, data_cal.data[ishot].time_refrac, data_obs.data[ishot].rec_refr_type, data_obs.data[ishot].rec_refr_index
					, data_obs.data[ishot].recx_refrac, data_obs.data[ishot].recy_refrac, model.xCen, model.yCen, model.intf_ynum, model.dX, model.dY, model.xTopo
					, model.Slowness, ishot, T_refrac, T_reflec, grad_refr_tmp, dens_refr_tmp, 0);
			}
			if (data_cal.data[ishot].nrec_reflec > 0)
			{
				_get_grad_1shot(data_cal.data[ishot].nrec_reflec, nx, ny, nx_s, nx_e, data_obs.data[ishot].err_reflec
					, data_obs.data[ishot].time_reflec, data_cal.data[ishot].time_reflec, data_obs.data[ishot].rec_refl_type, data_obs.data[ishot].rec_refl_index
					, data_obs.data[ishot].recx_reflec, data_obs.data[ishot].recy_reflec, model.xCen, model.yCen, model.intf_ynum, model.dX, model.dY, model.xTopo
					, model.Slowness, ishot, T_reflec, T_refrac, grad_refl_tmp, dens_refl_tmp, 1);
			}
			for (int ix = 0; ix < nx; ix++)
				for (int iy = 0; iy < ny; iy++)
				{
					if (data_cal.data[ishot].nrec_refrac > 0)grad_refr[ix][iy] -= grad_refr_tmp[ix][iy] / pow(model.Vel[ix][iy], 3);
					if (data_cal.data[ishot].nrec_reflec > 0)grad_refl[ix][iy] -= grad_refl_tmp[ix][iy] / pow(model.Vel[ix][iy], 3);
					if (precondition)
					{
						if (data_cal.data[ishot].nrec_refrac > 0)dens_refr[ix][iy] += dens_refr_tmp[ix][iy] / pow(model.Vel[ix][iy], 3);
						if (data_cal.data[ishot].nrec_reflec > 0)dens_refl[ix][iy] += dens_refl_tmp[ix][iy] / pow(model.Vel[ix][iy], 3);
					}
				}
		}
	}
	if (inv)
	{
		if (precondition)
		{
			for (int ix = 0; ix < nx; ix++)
				for (int iy = 0; iy < ny; iy++)
				{
					if (iy < model.topo_ynum[ix] + 1 || ix == nx - 1)
					{
						dens_refr[ix][iy] = 0;
						dens_refl[ix][iy] = 0;
					}
					if (iy > model.intf_ynum[ix] - 1)
					{
						dens_refl[ix][iy] = 0;
					}
				}
			normolization(dens_refr, nx, ny);
			normolization(dens_refl, nx, ny);
		}
		for (int ix = 0; ix < nx; ix++)
			for (int iy = 0; iy < ny; iy++)
			{
				if (iy < model.topo_ynum[ix] + 1 || ix == nx - 1)
				{
					grad_refr[ix][iy] = 0;
					grad_refl[ix][iy] = 0;
				}
				if (iy > model.intf_ynum[ix] - 1)
				{
					grad_refl[ix][iy] = 0;
				}
				if (precondition)
				{
					grad_refr[ix][iy] /= (dens_refr[ix][iy] + damp);
					grad_refl[ix][iy] /= (dens_refl[ix][iy] + damp);
				}
			}
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 0)
	{
		cerr << endl;
		cerr << "Forward modeling time: " << duration.count() << " s" << endl;
	}
	free_2d(grad_refl_tmp, nx);
	free_2d(grad_refr_tmp, nx);
	free_2d(dens_refr_tmp, nx);
	free_2d(dens_refl_tmp, nx);
	free_2d(dens_refl, nx);
	free_2d(dens_refr, nx);
	free_2d(T_refrac, nx);
	free_2d(T_reflec, nx);
	free_1d(tt_interface);
	///
	///Ruduce_and_Bcast
	for (int ishot = 0; ishot < data_cal.Num_source; ishot++)
	{
		for (int irec = 0; irec < data_cal.data[ishot].nrec_refrac; irec++)
		{
			void* temp_ptr = malloc(sizeof(double));
			double* new_sendbuf = (double*)temp_ptr;
			memcpy(new_sendbuf, &data_cal.data[ishot].time_refrac[irec], sizeof(double));
			MPI_Reduce(new_sendbuf, &data_cal.data[ishot].time_refrac[irec],
				1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			free(new_sendbuf);
		}
		for (int irec = 0; irec < data_cal.data[ishot].nrec_reflec; irec++)
		{
			void* temp_ptr = malloc(sizeof(double));
			double* new_sendbuf = (double*)temp_ptr;
			memcpy(new_sendbuf, &data_cal.data[ishot].time_reflec[irec], sizeof(double));
			MPI_Reduce(new_sendbuf, &data_cal.data[ishot].time_reflec[irec],
				1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			free(new_sendbuf);
		}
	}
	for (int ishot = 0; ishot < data_cal.Num_source; ishot++)
	{
		for (int irec = 0; irec < data_cal.data[ishot].nrec_refrac; irec++)
			MPI_Bcast(&data_cal.data[ishot].time_refrac[irec], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int irec = 0; irec < data_cal.data[ishot].nrec_reflec; irec++)
			MPI_Bcast(&data_cal.data[ishot].time_reflec[irec], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	if (inv)
	{
		for (int ix = 0; ix < nx; ix++)
		{
			void* temp_ptr = malloc(sizeof(double) * ny);
			void* temp_ptr2 = malloc(sizeof(double) * ny);
			double* new_sendbuf = (double*)temp_ptr;
			double* new_sendbuf2 = (double*)temp_ptr2;
			memcpy(new_sendbuf, grad_refr[ix], sizeof(double) * ny);
			memcpy(new_sendbuf2, grad_refl[ix], sizeof(double) * ny);
			MPI_Reduce(new_sendbuf, grad_refr[ix], ny, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(new_sendbuf2, grad_refl[ix], ny, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			free(new_sendbuf);
			free(new_sendbuf2);
		}
		for (int ix = 0; ix < nx; ix++)
		{
			MPI_Bcast(grad_refr[ix], ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(grad_refl[ix], ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Tomofsm2D::_Fsm_fwd_1shot_ttfield(int nx, int ny, int nrecfr, int nrecfl, int nx_s, int nx_e,
	double soux, double souy, vector< vector<int> >  sou_index, vector<int> sou_type, vector<int> intf_ynum,
	double** T1, double* tt_interface, vector<double> dX, vector<double> dY, vector< vector<double> >  Slowness,
	vector<double> xCen, vector<double> yCen, vector< vector<int> >  rec_refr_index, vector< vector<int> >  rec_refl_index,
	int raycode, int ishot)//1shot ttfield
{
	int isweep;
	int a, b, c, d;
	double Tx_min, Tz_min, T2;
	for (int ix = 0; ix < nx; ix++)
		for (int iy = 0; iy < ny; iy++)
		{
			T1[ix][iy] = 999999.;
		}
	///block///
	//////
	isweep = 0;
	while (1)
	{
		isweep++;
		if (isweep > nsweep)break;
		//first quadrant
		///add source point
		if (raycode == 0)
		{
			_add_source(nx, ny, T1, soux, souy, sou_index[0][ishot], sou_index[1][ishot]
				, sou_type[ishot], xCen, yCen, Slowness);
		}
		else
		{
			_add_interface(T1, tt_interface, nx, intf_ynum);
		}
		for (int ix = nx_s; ix < nx_e; ix++)
			for (int iy = 0; iy < ny; iy++)
			{
				_upa(&a, &b, &c, &d, ix, iy, nx, ny);
				Tx_min = min(T1[a][iy], T1[b][iy]);
				Tz_min = min(T1[ix][c], T1[ix][d]);
				T2 = _fin_min(Tx_min, Tz_min, dX[ix], dY[iy], Slowness[ix][iy]);
				T1[ix][iy] = min(T1[ix][iy], T2);
			}
		//second quadrant
		///add source point
		if (raycode == 0)
		{
			_add_source(nx, ny, T1, soux, souy, sou_index[0][ishot], sou_index[1][ishot]
				, sou_type[ishot], xCen, yCen, Slowness);
		}
		else
		{
			_add_interface(T1, tt_interface, nx, intf_ynum);
		}
		for (int ix = nx_e - 1; ix >= nx_s; ix--)
			for (int iy = 0; iy < ny; iy++)
			{
				_upa(&a, &b, &c, &d, ix, iy, nx, ny);
				Tx_min = min(T1[a][iy], T1[b][iy]);
				Tz_min = min(T1[ix][c], T1[ix][d]);
				T2 = _fin_min(Tx_min, Tz_min, dX[ix], dY[iy], Slowness[ix][iy]);
				T1[ix][iy] = min(T1[ix][iy], T2);
			}
		//third quadrant
		///add source point
		if (raycode == 0)
		{
			_add_source(nx, ny, T1, soux, souy, sou_index[0][ishot], sou_index[1][ishot]
				, sou_type[ishot], xCen, yCen, Slowness);
		}
		else
		{
			_add_interface(T1, tt_interface, nx, intf_ynum);
		}
		for (int ix = nx_e - 1; ix >= nx_s; ix--)
			for (int iy = ny - 1; iy >= 0; iy--)
			{
				_upa(&a, &b, &c, &d, ix, iy, nx, ny);
				Tx_min = min(T1[a][iy], T1[b][iy]);
				Tz_min = min(T1[ix][c], T1[ix][d]);
				T2 = _fin_min(Tx_min, Tz_min, dX[ix], dY[iy], Slowness[ix][iy]);
				T1[ix][iy] = min(T1[ix][iy], T2);
			}
		//fouth quadrant
		///add source point
		if (raycode == 0)
		{
			_add_source(nx, ny, T1, soux, souy, sou_index[0][ishot], sou_index[1][ishot]
				, sou_type[ishot], xCen, yCen, Slowness);
		}
		else
		{
			_add_interface(T1, tt_interface, nx, intf_ynum);
		}
		for (int ix = nx_s; ix < nx_e; ix++)
			for (int iy = ny - 1; iy >= 0; iy--)
			{
				_upa(&a, &b, &c, &d, ix, iy, nx, ny);
				Tx_min = min(T1[a][iy], T1[b][iy]);
				Tz_min = min(T1[ix][c], T1[ix][d]);
				T2 = _fin_min(Tx_min, Tz_min, dX[ix], dY[iy], Slowness[ix][iy]);
				T1[ix][iy] = min(T1[ix][iy], T2);
			}
	}
}

void  Tomofsm2D::_upa(int* a, int* b, int* c, int* d, int ix, int iy, int nx, int nz)
{
	if (ix == 0 && iy == 0)
	{
		*a = ix; *b = ix + 1; *c = iy; *d = iy + 1;
	}
	if (iy == 0 && ix > 0 && ix < nx - 1)
	{
		*a = ix - 1; *b = ix + 1; *c = iy; *d = iy + 1;
	}
	if (ix == 0 && iy > 0 && iy < nz - 1)
	{
		*a = ix; *b = ix + 1; *c = iy - 1; *d = iy + 1;
	}
	if (ix == 0 && iy == (nz - 1))
	{
		*a = ix; *b = ix + 1; *c = iy; *d = iy - 1;
	}
	if (ix == (nx - 1) && iy == 0)
	{
		*a = ix; *b = ix - 1; *c = iy; *d = iy + 1;
	}
	if (ix == (nx - 1) && iy > 0 && iy < nz - 1)
	{
		*a = ix; *b = ix - 1; *c = iy - 1; *d = iy + 1;
	}
	if (ix == (nx - 1) && iy == (nz - 1))
	{
		*a = ix; *b = ix - 1; *c = iy; *d = iy - 1;
	}
	if (iy == (nz - 1) && ix > 0 && ix < (nx - 1))
	{
		*a = ix + 1; *b = ix - 1; *c = iy; *d = iy - 1;
	}
	if (ix > 0 && ix < nx - 1 && iy>0 && iy < nz - 1)
	{
		*a = ix - 1; *b = ix + 1; *c = iy - 1; *d = iy + 1;
	}
}

double Tomofsm2D::_fin_min(double Tx_min, double Tz_min, double dx, double dz, double s)
{
	double tmp;
	if (abs(Tx_min - Tz_min) >= dx * s)
		tmp = min(Tx_min, Tz_min) + dx * s;
	else
		tmp = (Tx_min + Tz_min + sqrt(2 * dx * dx * s * s - pow(double(Tx_min - Tz_min), double(2.0)))) / 2.0;
	return tmp;
}

void Tomofsm2D::_add_source(int nx, int ny, double** T, double soux, double souy, int sx_index, int sy_index,
	int stype, vector<double> xCen, vector<double> yCen, vector< vector<double> >  Slowness)
{
	if (stype == 0)//incell
	{
		double t5, t6, t7, t8, t9, t10, t11, t12;
		double dis5, dis6, dis7, dis8, dis9, dis10, dis11, dis12;
		double dis1 = dist(soux, xCen[sx_index], souy, yCen[sy_index]);
		double t1 = dis1 * Slowness[sx_index][sy_index];
		double dis2 = dist(soux, xCen[sx_index], souy, yCen[sy_index + 1]);
		double t2 = dis2 * Slowness[sx_index][sy_index + 1];
		double dis3 = dist(soux, xCen[sx_index + 1], souy, yCen[sy_index + 1]);
		double t3 = dis3 * Slowness[sx_index + 1][sy_index + 1];
		double dis4 = dist(soux, xCen[sx_index + 1], souy, yCen[sy_index]);
		double t4 = dis4 * Slowness[sx_index + 1][sy_index];
		T[sx_index][sy_index] = t1;
		T[sx_index][sy_index + 1] = t2;
		T[sx_index + 1][sy_index + 1] = t3;
		T[sx_index + 1][sy_index] = t4;
		if (sx_index > 0 && sy_index - 1 > 0)
		{
			dis5 = dist(xCen[sx_index], xCen[sx_index], yCen[sy_index - 1], yCen[sy_index]);
			t5 = dis5 * Slowness[sx_index][sy_index - 1];
			T[sx_index][sy_index - 1] = t5 + t1;
		}
		if (sx_index - 1 > 0 && sy_index > 0)
		{
			dis6 = dist(xCen[sx_index - 1], xCen[sx_index], yCen[sy_index], yCen[sy_index]);
			t6 = dis6 * Slowness[sx_index - 1][sy_index];
			T[sx_index - 1][sy_index] = t6 + t1;
		}
		if (sx_index - 1 > 0 && sy_index + 1 < ny)
		{
			dis7 = dist(xCen[sx_index - 1], xCen[sx_index], yCen[sy_index + 1], yCen[sy_index + 1]);
			t7 = dis7 * Slowness[sx_index - 1][sy_index + 1];
			T[sx_index - 1][sy_index + 1] = t7 + t2;
		}
		if (sx_index > 0 && sy_index + 2 < ny)
		{
			dis8 = dist(xCen[sx_index], xCen[sx_index], yCen[sy_index + 1], yCen[sy_index + 2]);
			t8 = dis8 * Slowness[sx_index][sy_index + 2];
			T[sx_index][sy_index + 2] = t8 + t2;
		}
		if (sx_index + 1 < nx && sy_index + 2 < ny)
		{
			dis9 = dist(xCen[sx_index + 1], xCen[sx_index + 1], yCen[sy_index + 1], yCen[sy_index + 2]);
			t9 = dis9 * Slowness[sx_index + 1][sy_index + 2];
			T[sx_index + 1][sy_index + 2] = t9 + t3;
		}
		if (sx_index + 2 < nx && sy_index + 1 < ny)
		{
			dis10 = dist(xCen[sx_index + 2], xCen[sx_index + 1], yCen[sy_index + 1], yCen[sy_index + 1]);
			t10 = dis10 * Slowness[sx_index + 2][sy_index + 1];
			T[sx_index + 2][sy_index + 1] = t10 + t3;
		}
		if (sx_index + 2 < nx && sy_index < ny)
		{
			dis11 = dist(xCen[sx_index + 2], xCen[sx_index + 1], yCen[sy_index], yCen[sy_index]);
			t11 = dis11 * Slowness[sx_index + 2][sy_index];
			T[sx_index + 2][sy_index] = t11 + t4;
		}
		if (sx_index + 1 < nx && sy_index - 1 > 0)
		{
			dis12 = dist(xCen[sx_index], xCen[sx_index], yCen[sy_index], yCen[sy_index - 1]);
			t12 = dis12 * Slowness[sx_index + 1][sy_index - 1];
			T[sx_index + 1][sy_index - 1] = t12 + t4;
		}
	}
	else if (stype > 4)//inpoint
	{
		T[sx_index][sy_index] = 0;
	}
	else
	{
		if (stype == 1 || stype == 3)
		{
			double dis1 = dist(soux, xCen[sx_index], souy, yCen[sy_index]);
			double t1 = dis1 * Slowness[sx_index][sy_index];
			double dis2 = dist(soux, xCen[sx_index], souy, yCen[sy_index + 1]);
			double t2 = dis2 * Slowness[sx_index][sy_index + 1];
			T[sx_index][sy_index] = t1;
			T[sx_index][sy_index + 1] = t2;
		}
		else if (stype == 2 || stype == 4)
		{
			double dis1 = dist(soux, xCen[sx_index], souy, yCen[sy_index]);
			double t1 = dis1 * Slowness[sx_index][sy_index];
			double dis4 = dist(soux, xCen[sx_index + 1], souy, yCen[sy_index]);
			double t4 = dis4 * Slowness[sx_index + 1][sy_index];
			T[sx_index][sy_index] = t1;
			T[sx_index + 1][sy_index] = t4;
		}
		else
		{
			cerr << "Unrecognized source index type! Please check the source point" << endl;
			exit(1);
		}
	}
}

void Tomofsm2D::_add_interface(double** T, double* tt_interface, int nx, vector<int> intf_ynum)
{
	for (int ix = 0; ix < nx; ix++)
	{
		T[ix][intf_ynum[ix]] = tt_interface[ix];
	}
}

void Tomofsm2D::_ttfield_print(double** T, int ishot, int raycode, TomoMesh2D model)
{
	string filename;
	if (raycode == 0)
	{
		filename = refrac_outroot + "tt_" + to_string(ishot) + ".dat";
	}
	else
	{
		filename = reflec_outroot + "tt_" + to_string(ishot) + ".dat";
	}

	if (filetype == "binary")
	{
		write_allocdouble_2d_bin(filename, T, model.Num_xCen, model.Num_yCen);
	}
	else
	{
		ofstream fout(filename);
		if (!fout)
		{
			cerr << "Output travel-time field failed!" << endl;
			exit(1);
		}
		for (int ix = 0; ix < model.Num_xCen; ix++)
			for (int iy = 0; iy < model.Num_yCen; iy++)
			{
				fout << model.xCen[ix] / 1000. << " " << model.yCen[iy] / 1000. << " " << T[ix][iy] << endl;
			}
		fout.close();
	}
}

void Tomofsm2D::gradient_print(double** grad, int raycode, int n, TomoMesh2D model)
{
	string filename;
	if (raycode == 0)
	{
		filename = refrac_outroot + "grad_" + to_string(n) + ".dat";
	}
	else
	{
		filename = reflec_outroot + "grad_" + to_string(n) + ".dat";
	}

	if (filetype == "binary")
	{
		write_allocdouble_2d_bin(filename, grad, model.Num_xCen, model.Num_yCen);
	}
	else
	{
		ofstream fout(filename);
		if (!fout)
		{
			cerr << "Output gradient  failed!" << endl;
			exit(1);
		}
		for (int ix = 0; ix < model.Num_xCen; ix++)
			for (int iy = 0; iy < model.Num_yCen; iy++)
			{
				fout << model.xCen[ix] / 1000. << " " << model.yCen[iy] / 1000. << " " << grad[ix][iy] << endl;
			}
		fout.close();
	}
}

double Tomofsm2D::_cal_rec_tt(double** T, int ishot, int irec, int x_index, int y_index, double recx, double recy, int type,
	vector<double>xCen, vector<double>yCen)
{
	double t;
	if (type == 0)
	{
		double x_weight = (recx - xCen[x_index]) / (xCen[x_index + 1] - xCen[x_index]);
		double y_weight = (recy - yCen[y_index]) / (yCen[y_index + 1] - yCen[y_index]);
		t = (1 - x_weight) * (1 - y_weight) * T[x_index][y_index] + x_weight * (1 - y_weight) * T[x_index + 1][y_index] +
			(1 - x_weight) * y_weight * T[x_index][y_index + 1] + x_weight * y_weight * T[x_index + 1][y_index + 1];
		if (t < 0 || t > 9990)
		{
			cerr << "travel time error!" << endl;
			exit(1);
		}
		return t;
	}
	if (type == 1 || type == 3)//lb&rb
	{
		t = T[x_index][y_index] + (T[x_index][y_index + 1] - T[x_index][y_index])
			* (recy - yCen[y_index]) / (yCen[y_index + 1] - yCen[y_index]);
		if (t < 0 || t > 9990)
		{
			cerr << "travel time error!" << endl;
			exit(1);
		}
		return t;
	}
	if (type == 2 || type == 4)//db&ub
	{
		t = T[x_index][y_index] + (T[x_index + 1][y_index] - T[x_index][y_index])
			* (recx - xCen[x_index]) / (xCen[x_index + 1] - xCen[x_index]);
		if (t < 0 || t > 9990)
		{
			cerr << "travel time error!" << endl;
			exit(1);
		}
		return t;
	}
	else
	{
		t = T[x_index][y_index];
		if (t < 0 || t > 9990)
		{
			cerr << "travel time error!" << endl;
			exit(1);
		}
		return t;
	}
	return 0;
}
//inverse
void Tomofsm2D::_get_grad_1shot(int nrec, int nx, int ny, int nx_s, int nx_e, vector<double> err1, vector<double> tobs, vector<double> tcal, vector<int> type1,
	vector< vector<int> > index, vector<double> recx1, vector<double> recy1, vector<double> xCen, vector<double> yCen, vector<int> intf_ynum,
	vector<double> dX, vector<double> dY, vector<double> xTopo, vector< vector<double> > Slowness,
	int ishot, double** T, double** T2, double** grad, double** dens, int raycode)
{
	int type;
	int* x_index, * y_index;
	double** lamda_1, ** lamda_2, * lamda_tmp_1, * lamda_tmp_2;
	double** lamda_3, ** lamda_4, * lamda_tmp_3, * lamda_tmp_4;
	double** lamda_1_pre, ** lamda_2_pre, * lamda_tmp_1_pre, * lamda_tmp_2_pre;
	double** lamda_3_pre, ** lamda_4_pre, * lamda_tmp_3_pre, * lamda_tmp_4_pre;
	int** lamda_tmp_index;
	lamda_1 = alloc_double_2d(nx, ny);
	lamda_2 = alloc_double_2d(nx, ny);
	lamda_3 = alloc_double_2d(nx, ny);
	lamda_4 = alloc_double_2d(nx, ny);
	lamda_tmp_index = alloc_int_2d(2, nrec);
	lamda_tmp_1 = alloc_double_1d(nrec);
	lamda_tmp_2 = alloc_double_1d(nrec);
	lamda_tmp_3 = alloc_double_1d(nx);
	lamda_tmp_4 = alloc_double_1d(nx);
	lamda_1_pre = alloc_double_2d(nx, ny);
	lamda_2_pre = alloc_double_2d(nx, ny);
	lamda_3_pre = alloc_double_2d(nx, ny);
	lamda_4_pre = alloc_double_2d(nx, ny);
	lamda_tmp_1_pre = alloc_double_1d(nrec);
	lamda_tmp_2_pre = alloc_double_1d(nrec);
	lamda_tmp_3_pre = alloc_double_1d(nx);
	lamda_tmp_4_pre = alloc_double_1d(nx);
	zero_2d(lamda_1, nx, ny); zero_2d(lamda_2, nx, ny);
	zero_2d(lamda_3, nx, ny); zero_2d(lamda_4, nx, ny);
	zero_2d(lamda_1_pre, nx, ny); zero_2d(lamda_2_pre, nx, ny);
	zero_2d(lamda_3_pre, nx, ny); zero_2d(lamda_4_pre, nx, ny);
	zero_1d(lamda_tmp_1_pre, nrec); zero_1d(lamda_tmp_2_pre, nrec);
	zero_1d(lamda_tmp_3_pre, nx); zero_1d(lamda_tmp_4_pre, nx);
	zero_2i(lamda_tmp_index, 2, nrec);
	double delta_tx, delta_tz, lamda_tmp, lamda_tmp_pre, angle, err;
	int a1, a2, b1, b2, ix, iy, ir;
	a2 = 9999; a1 = 9999;
	///
	for (ir = 0; ir < nrec; ir++)
	{
		//type
		err = err1[ir];
		type = type1[ir];
		//
		x_index = alloc_int_1d(4);
		y_index = alloc_int_1d(4);
		int x1, y1;
		double recy, recx, t_obs, t_cal;
		x1 = index[0][ir]; y1 = index[1][ir];
		recx = recx1[ir]; recy = recy1[ir];
		t_obs = tobs[ir]; t_cal = tcal[ir];
		double dis, dismin, dt, ddt;
		dismin = 99999999.;
		x_index[0] = x1; x_index[1] = x1;
		x_index[2] = x1 + 1; x_index[3] = x1 + 1;
		y_index[0] = y1; y_index[1] = y1 + 1;
		y_index[2] = y1 + 1; y_index[3] = y1;
		for (int i = 0; i < 4; i++)
		{
			dis = dist(xCen[x_index[i]], recx, yCen[y_index[i]], recy);
			if (dis < dismin)
			{
				lamda_tmp_index[0][ir] = x_index[i];
				lamda_tmp_index[1][ir] = y_index[i];
				dismin = dis;
			}
		}
		if (t_cal - T[lamda_tmp_index[0][ir]][lamda_tmp_index[1][ir]] != 0)
			ddt = abs(t_cal - T[lamda_tmp_index[0][ir]][lamda_tmp_index[1][ir]]) / (t_cal - T[lamda_tmp_index[0][ir]][lamda_tmp_index[1][ir]]);
		else
			ddt = 0;
		dt = T[lamda_tmp_index[0][ir]][lamda_tmp_index[1][ir]]
			- (t_obs - ddt * Slowness[lamda_tmp_index[0][ir]][lamda_tmp_index[1][ir]] * dis);
		a1 = lamda_tmp_index[0][ir];
		a2 = lamda_tmp_index[1][ir];
		b1 = a2 + 1;
		b2 = a2;
		if (b1 >= ny - 1)
		{
			b1 = ny - 1;
			b2 = ny - 2;
		}
		delta_tz = (T[a1][b1] - T[a1][b2]) / (yCen[b1] - yCen[b2]);
		b1 = a1 + 1;
		b2 = a1;
		if (b1 >= nx - 1)
		{
			b1 = nx - 1;
			b2 = nx - 2;
		}
		delta_tx = (T[b1][a2] - T[b2][a2]) / (xCen[b1] - xCen[b2]);
		//calculate angle
		double ddx = xCen[b1] - xCen[b2];
		double ddy = xTopo[b1] - xTopo[b2];
		double trd = sqrt(ddx * ddx + ddy * ddy);
		double cos = ddy / trd;
		double sin = ddx / trd;
		if (cos * delta_tx + sin * delta_tz == 0)
			lamda_tmp = 0;
		else
		{
			lamda_tmp = 1. / err * dt / (cos * delta_tx + sin * delta_tz);
			lamda_tmp_pre = -1. / err / (cos * delta_tx + sin * delta_tz);
		}
		//
		if (lamda_tmp < 0)
		{
			lamda_tmp_1[ir] = lamda_tmp;
			lamda_tmp_1_pre[ir] = lamda_tmp_pre;
		}
		else
		{
			lamda_tmp_2[ir] = lamda_tmp;
			lamda_tmp_2_pre[ir] = lamda_tmp_pre;
		}
		free_1i(x_index);	free_1i(y_index);
	}
	///block///
	//////
	_calc_lamda(nrec, nx, ny, nx_s, nx_e, dX, dY, intf_ynum,
		T, lamda_1, lamda_tmp_1, lamda_tmp_index, ishot, 1, raycode, 0);
	_calc_lamda(nrec, nx, ny, nx_s, nx_e, dX, dY, intf_ynum,
		T, lamda_2, lamda_tmp_2, lamda_tmp_index, ishot, 2, raycode, 0);
	if (precondition)
	{
		_calc_lamda(nrec, nx, ny, nx_s, nx_e, dX, dY, intf_ynum,
			T, lamda_1_pre, lamda_tmp_1_pre, lamda_tmp_index, ishot, 1, raycode, 0);
		_calc_lamda(nrec, nx, ny, nx_s, nx_e, dX, dY, intf_ynum,
			T, lamda_2_pre, lamda_tmp_2_pre, lamda_tmp_index, ishot, 2, raycode, 0);
	}
	if (raycode == 1)///source point gradient
	{
		for (ix = nx_s; ix < nx_e; ix++)
		{
			lamda_tmp_3[ix] = lamda_1[ix][intf_ynum[ix] - 1];
			lamda_tmp_4[ix] = lamda_2[ix][intf_ynum[ix] - 1];
			if (precondition)
			{
				lamda_tmp_3_pre[ix] = lamda_1_pre[ix][intf_ynum[ix] - 1];
				lamda_tmp_4_pre[ix] = lamda_2_pre[ix][intf_ynum[ix] - 1];
			}
		}
		_calc_lamda(nrec, nx, ny, nx_s, nx_e, dX, dY, intf_ynum,
			T2, lamda_3, lamda_tmp_3, lamda_tmp_index, ishot, 1, raycode, 1);
		_calc_lamda(nrec, nx, ny, nx_s, nx_e, dX, dY, intf_ynum,
			T2, lamda_4, lamda_tmp_4, lamda_tmp_index, ishot, 2, raycode, 1);
		if (precondition)
		{
			_calc_lamda(nrec, nx, ny, nx_s, nx_e, dX, dY, intf_ynum,
				T2, lamda_3_pre, lamda_tmp_3_pre, lamda_tmp_index, ishot, 1, raycode, 1);
			_calc_lamda(nrec, nx, ny, nx_s, nx_e, dX, dY, intf_ynum,
				T2, lamda_4_pre, lamda_tmp_4_pre, lamda_tmp_index, ishot, 2, raycode, 1);
		}
	}
	for (ix = nx_s; ix < nx_e; ix++)
		for (iy = 0; iy < ny; iy++)
		{
			grad[ix][iy] = lamda_1[ix][iy] + lamda_2[ix][iy] + lamda_3[ix][iy] + lamda_4[ix][iy];
			if (precondition)dens[ix][iy] = abs(lamda_1_pre[ix][iy]) + abs(lamda_2_pre[ix][iy]) + abs(lamda_3_pre[ix][iy]) + abs(lamda_4_pre[ix][iy]);
		}
	free_2d(lamda_1, nx);	free_2d(lamda_2, nx);
	free_2d(lamda_3, nx);	free_2d(lamda_4, nx);
	free_1d(lamda_tmp_1);	free_1d(lamda_tmp_2);
	free_1d(lamda_tmp_3);	free_1d(lamda_tmp_4);
	free_2d(lamda_1_pre, nx);	free_2d(lamda_2_pre, nx);
	free_2d(lamda_3_pre, nx);	free_2d(lamda_4_pre, nx);
	free_1d(lamda_tmp_1_pre);	free_1d(lamda_tmp_2_pre);
	free_1d(lamda_tmp_3_pre);	free_1d(lamda_tmp_4_pre);
	free_2i(lamda_tmp_index, 2);
}

void Tomofsm2D::_calc_lamda(int nrec, int nx, int ny, int nx_s, int nx_e, vector<double>dx, vector<double>dy, vector<int> intf_ynum,
	double** T1, double** lamda, double* lamda_tmp, int** lamda_type_nidex, int ishot, int minmax, int raycode, int n2s)
{
	int isweep;
	double tmp, dx1, dx2, dy1, dy2;
	int ix, iy, ir;
	zero_2d(lamda, nx, ny);
	double a_p_1, a_p_2, a_r_1, a_r_2, b_p_1, b_p_2, b_r_1, b_r_2;
	for (isweep = 0; isweep < nsweep; isweep++)
	{
		if (n2s == 0)
		{
			for (ir = 0; ir < nrec; ir++)lamda[lamda_type_nidex[0][ir]][lamda_type_nidex[1][ir]] = lamda_tmp[ir];
		}
		else
		{
			for (ix = 0; ix < nx; ix++)lamda[ix][intf_ynum[ix]] = lamda_tmp[ix];
		}
		tmp = 0.0;
		/**********************1~I;1~J**************************************************/
		for (ix = 1 + nx_s; ix < nx_e - 1; ix++)
		{
			for (iy = 1; iy < ny - 1; iy++)
			{
				dx1 = dx[ix - 1];
				dy1 = dy[iy - 1];
				dx2 = dx[ix + 1];
				dy2 = dy[iy + 1];
				_upab(dx1, dx2, dy1, dy2, T1, a_p_1, a_p_2, a_r_1, a_r_2, b_p_1, b_p_2, b_r_1, b_r_2, ix, iy);
				_uplamda(dx1, dx2, dy1, dy2, lamda, a_p_1, a_p_2, a_r_1, a_r_2, b_p_1, b_p_2, b_r_1, b_r_2, ix, iy, minmax);
			}
		}
		if (n2s == 0)
		{
			for (ir = 0; ir < nrec; ir++)lamda[lamda_type_nidex[0][ir]][lamda_type_nidex[1][ir]] = lamda_tmp[ir];
		}
		else
		{
			for (ix = 0; ix < nx; ix++)lamda[ix][intf_ynum[ix]] = lamda_tmp[ix];
		}
		/***********************1~I;J~1**************************************************/
		for (ix = nx_s + 1; ix < nx_e - 1; ix++)
			for (iy = ny - 2; iy > 0; iy--)
			{
				dx1 = dx[ix - 1];
				dy1 = dy[iy - 1];
				dx2 = dx[ix + 1];
				dy2 = dy[iy + 1];
				_upab(dx1, dx2, dy1, dy2, T1, a_p_1, a_p_2, a_r_1, a_r_2, b_p_1, b_p_2, b_r_1, b_r_2, ix, iy);
				_uplamda(dx1, dx2, dy1, dy2, lamda, a_p_1, a_p_2, a_r_1, a_r_2, b_p_1, b_p_2, b_r_1, b_r_2, ix, iy, minmax);
			}
		if (n2s == 0)
		{
			for (ir = 0; ir < nrec; ir++)lamda[lamda_type_nidex[0][ir]][lamda_type_nidex[1][ir]] = lamda_tmp[ir];
		}
		else
		{
			for (ix = 0; ix < nx; ix++)lamda[ix][intf_ynum[ix]] = lamda_tmp[ix];
		}
		/**********************I~1;1~J***************************************************/
		for (ix = nx_e - 2; ix > nx_s; ix--)
			for (iy = 1; iy < ny - 1; iy++)
			{
				dx1 = dx[ix - 1];
				dy1 = dy[iy - 1];
				dx2 = dx[ix + 1];
				dy2 = dy[iy + 1];
				_upab(dx1, dx2, dy1, dy2, T1, a_p_1, a_p_2, a_r_1, a_r_2, b_p_1, b_p_2, b_r_1, b_r_2, ix, iy);
				_uplamda(dx1, dx2, dy1, dy2, lamda, a_p_1, a_p_2, a_r_1, a_r_2, b_p_1, b_p_2, b_r_1, b_r_2, ix, iy, minmax);
			}
		if (n2s == 0)
		{
			for (ir = 0; ir < nrec; ir++)lamda[lamda_type_nidex[0][ir]][lamda_type_nidex[1][ir]] = lamda_tmp[ir];
		}
		else
		{
			for (ix = 0; ix < nx; ix++)lamda[ix][intf_ynum[ix]] = lamda_tmp[ix];
		}
		/***********************I~1;J~1**************************************************/
		for (ix = nx_e - 2; ix > nx_s; ix--)
			for (iy = ny - 2; iy > 0; iy--)
			{
				dx1 = dx[ix - 1];
				dy1 = dy[iy - 1];
				dx2 = dx[ix + 1];
				dy2 = dy[iy + 1];
				_upab(dx1, dx2, dy1, dy2, T1, a_p_1, a_p_2, a_r_1, a_r_2, b_p_1, b_p_2, b_r_1, b_r_2, ix, iy);
				_uplamda(dx1, dx2, dy1, dy2, lamda, a_p_1, a_p_2, a_r_1, a_r_2, b_p_1, b_p_2, b_r_1, b_r_2, ix, iy, minmax);
				tmp += lamda[ix][iy];
			}
		if (n2s == 0)
		{
			for (ir = 0; ir < nrec; ir++)lamda[lamda_type_nidex[0][ir]][lamda_type_nidex[1][ir]] = lamda_tmp[ir];
		}
		else
		{
			for (ix = 0; ix < nx; ix++)lamda[ix][intf_ynum[ix]] = lamda_tmp[ix];
		}
	}
}

void Tomofsm2D::_upab(double dx1, double dx2, double dy1, double dy2, double** T1, double& a_p_1, double& a_p_2, double& a_r_1,
	double& a_r_2, double& b_p_1, double& b_p_2, double& b_r_1, double& b_r_2, int ix, int iy)
{
	double a_r, a_p, b_r, b_p;
	a_r = (T1[ix][iy] - T1[ix - 1][iy]) / dx1;
	a_p = (T1[ix + 1][iy] - T1[ix][iy]) / dx2;
	b_r = (T1[ix][iy] - T1[ix][iy - 1]) / dy1;
	b_p = (T1[ix][iy + 1] - T1[ix][iy]) / dy2;

	a_p_1 = (a_p + abs(a_p)) * 0.5;
	a_p_2 = (a_p - abs(a_p)) * 0.5;
	a_r_1 = (a_r + abs(a_r)) * 0.5;
	a_r_2 = (a_r - abs(a_r)) * 0.5;
	b_p_1 = (b_p + abs(b_p)) * 0.5;
	b_p_2 = (b_p - abs(b_p)) * 0.5;
	b_r_1 = (b_r + abs(b_r)) * 0.5;
	b_r_2 = (b_r - abs(b_r)) * 0.5;
}

void Tomofsm2D::_uplamda(double dx1, double dx2, double dy1, double dy2, double** lamda, double a_p_1, double a_p_2, double a_r_1, double a_r_2,
	double b_p_1, double b_p_2, double b_r_1, double b_r_2, int ix, int iy, int minmax)
{
	double k1;
	double lamda_e = 0;
	double dx = (dx1 + dx2) * 0.5;
	double dy = (dy1 + dy2) * 0.5;
	k1 = (a_p_2 - a_r_1) / dx + (b_p_2 - b_r_1) / dy;
	if (k1 != 0)
		lamda_e = ((a_r_2 * lamda[ix - 1][iy] - a_p_1 * lamda[ix + 1][iy]) / dx + (b_r_2 * lamda[ix][iy - 1] - b_p_1 * lamda[ix][iy + 1]) / dy) / k1;
	if (minmax == 1)//min
	{
		if ((lamda_e - lamda[ix][iy]) >= 0)// lamda[ix][iy]<0
			lamda[ix][iy] = lamda[ix][iy];
		else
			lamda[ix][iy] = lamda_e;
	}
	else
	{
		if (lamda_e - lamda[ix][iy] <= 0)// lamda[ix][iy]>0
			lamda[ix][iy] = lamda[ix][iy];
		else
			lamda[ix][iy] = lamda_e;
	}
}
