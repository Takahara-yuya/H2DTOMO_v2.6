/*
***********************************************************************

H2DTOMO.cpp	(Main program)
This file is the main program of H2DTOMO.

***********************************************************************

Dec 01, 2023
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

* ***********************************************************************
*/
#include "../include/CommonHeaders.h"
#include "../include/Tomo_mesh2D.h"
#include "../include/Tomo_srs2D.h"
#include "../include/Tomo_fsm2D.h"
#include "../include/Tomo_inverse2D.h"
#include "../include/Tomo_fileread.h"
int main(int argc, char** argv)
{
	//mpi initialize//
	int myid, np;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	//
	if (myid == 0)
	{
		// Copyright and license notice:
		std::cout << std::endl;
		std::cout << "***************************************************************************" << std::endl << std::endl;
		std::cout << "                   H2DTOMO                                 " << std::endl << std::endl;
		std::cout << "  Copyright 2023, Zuwei Huang         " << std::endl << std::endl;

		std::cout << "  H2DTOMO is free software: you can redistribute it and/or modify it under the " << std::endl <<
			"  terms of the GNU Lesser General Public License as published by the Free " << std::endl <<
			"  Software Foundation, version 3 of the License. " << std::endl << std::endl;

		std::cout << "  H2DTOMO is distributed in the hope that it will be useful, but WITHOUT ANY " << std::endl <<
			"  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS " << std::endl <<
			"  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for " << std::endl <<
			"  more details. " << std::endl << std::endl;
		std::cout << "***************************************************************************" << std::endl << std::endl;
	}
	// The user should specify the configuration file as input to H2DTOMO:
	//************************** Read Startup File *************************************
	string startup = "Startup_Tomography";
	string Fwd_Setting_File, Inverse_Setting_File, Mesh_Setting_File, Observed_data_File;
	int fwd_only, raytype;
	TomoMesh2D model; TomoSrs2D data_obs; TomoSrs2D data_cal; Tomofsm2D fsm; Tomoinv2D inv;
	TomostartupRead(startup, raytype, fwd_only, Fwd_Setting_File, Inverse_Setting_File, Mesh_Setting_File, Observed_data_File);
	///MPI_INIT_IN_FWD_AND_INV
	fsm.myid = myid; fsm.np = np;
	inv.myid = myid; inv.np = np;
	///
	if (myid == 0)cerr << "Starting constructing velocity model!" << endl;
	model.Init_model(Mesh_Setting_File);
	if (myid == 0)std::cerr << "Constructing velocity model data sucessfully!" << endl;
	if (myid == 0)std::cerr << endl;
	/////
	if (myid == 0)std::cerr << "Starting reading observed data!" << endl;
	data_obs.Read_srsData(Observed_data_File, model);
	data_cal.Read_srsData(Observed_data_File, model);
	if (myid == 0)std::cerr << "Reading observed data sucessfully!" << endl;
	if (myid == 0)std::cerr << endl;
	/////
	if (raytype == 0)
	{
		if (myid == 0)std::cerr << "Raytype=0, only refraction!" << endl;
		for (int is = 0; is < data_obs.Num_source; is++)
		{
			data_obs.data[is].nrec_reflec = 0;
			data_cal.data[is].nrec_reflec = 0;
		}
	}
	else if (raytype == 1)
	{
		if (myid == 0)std::cerr << "Raytype=1, only reflection!" << endl;
		for (int is = 0; is < data_obs.Num_source; is++)
		{
			data_obs.data[is].nrec_refrac = 0;
			data_cal.data[is].nrec_refrac = 0;
		}
	}
	//
	if (myid == 0)
	{
		std::cerr << "Observed data reading finished!" << endl;
		std::cerr << endl;
		std::cerr << "NX: " << model.Num_xCen << " " << "NY:" << model.Num_yCen << endl;
		std::cerr << "XMIN XMAX: " << model.xMin << " " << model.xMax << " " << "YMIN YMAX: " << model.yMin << " " << model.yMax << endl;
		std::cerr << "Total shot: " << data_obs.Num_source << endl;
		std::cerr << endl;
	}
	/////
	if (fwd_only == 1)
	{
		fsm.Init_fsm(Fwd_Setting_File, model);
		fsm.Fsm_fwd(model, data_obs, data_cal, false);
		if (myid == 0)fsm.output(data_cal);
	}
	else
	{
		fsm.Init_fsm(Fwd_Setting_File, model);
		inv.Init_inv(Inverse_Setting_File);
		inv.ASTinv(model, model, data_obs, data_cal, fsm);
	}
	///
	MPI_Finalize();
}
