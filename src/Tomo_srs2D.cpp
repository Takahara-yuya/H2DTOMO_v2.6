/*
***********************************************************************

Tomo_srs2D.cpp	(Constructing Tomo srs)
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

#include "../include/Tomo_srs2D.h"

void TomoSrs2D::Read_srsData(string filename, TomoMesh2D model)
{
	Num_rec_refr = 0; Num_rec_refl = 0;
	double recx, recy, time, dt;
	int raycode;
	string tmp;
	ifstream fin(filename);
	if (!fin)
	{
		cerr << "Error reding srs file!" << endl;
	}
	fin >> Num_source;//number of source 
	data = new rec_data[Num_source];
	soux.resize(Num_source); souy.resize(Num_source);
	//
	for (int is = 0; is < Num_source; is++)
	{
		int ndata;
		fin >> tmp >> soux[is] >> souy[is] >> ndata;
		data[is].nrec_reflec = 0;	data[is].nrec_refrac = 0;
		soux[is] *= 1000.; souy[is] *= 1000.;
		for (int irec = 0; irec < ndata; irec++)
		{
			fin >> tmp >> recx >> recy >> raycode >> time >> dt;
			recx *= 1000.; recy *= 1000.;
			if (recx<model.xCen[0] || recx>model.xCen[model.Num_xCen - 1])
			{
				cerr << "receiver position exceed the model bound!" << endl;
				exit(1);
			}
			if (raycode == 0)//refraction
			{
				Num_rec_refr++;
				data[is].nrec_refrac++;
				data[is].recx_refrac.push_back(recx); data[is].recy_refrac.push_back(recy);
				data[is].time_refrac.push_back(time); data[is].err_refrac.push_back(dt);
			}
			else if (raycode == 1)//reflection
			{
				Num_rec_refl++;
				data[is].nrec_reflec++;
				data[is].recx_reflec.push_back(recx); data[is].recy_reflec.push_back(recy);
				data[is].time_reflec.push_back(time); data[is].err_reflec.push_back(dt);
			}
			else
			{
				cerr << raycode << " " << is << endl;
				cerr << "Input data error! Raycode recognize failed!" << endl;
				exit(1);
			}
		}
	}
	_Source_index(model);
	_Rec_index(model);
}

void TomoSrs2D::_Source_index(TomoMesh2D model)
{
	int x_index, y_index, ix, iy;
	bool lb, ub, rb, db;
	sou_index.resize(2);//x,y
	for (int i = 0; i < 2; i++)sou_index[i].resize(Num_source);
	sou_type.resize(Num_source);
	for (int is = 0; is < Num_source; is++)
	{
		//xindex search
		ix = 0; lb = false; rb = false; ub = false; db = false;
		if (souy[is] < model.yMin)souy[is] = model.yMin;
		if (souy[is] > model.yMax)souy[is] = model.yMax;
		while (1)
		{
			if (ix > model.Num_xCen - 1)
			{
				cerr << "source index failed! check the input data!" << endl;
				break;
			}
			if (soux[is] >= model.xCen[ix] && soux[is] < model.xCen[ix + 1])
			{
				sou_index[0][is] = ix;
				if (soux[is] == model.xCen[ix])lb = true;//leftbound
				else
				{
					lb = false; rb = false;
				}
				break;
			}
			else if (soux[is] == model.xCen[model.Num_xCen - 1])
			{
				sou_index[0][is] = model.Num_xCen - 1;
				rb = true;//rightbound
				break;
			}
			else
			{
				ix++;
			}
		}
		//yindex search
		iy = 0;
		while (1)
		{
			if (iy > model.Num_yCen - 1)
			{
				//	cerr<<souy[is]<<" "<<model.yCen[model.Num_yCen-1]<<endl;
				cerr << "source index failed! check the input data!" << endl;
				break;
			}
			if (souy[is] >= model.yCen[iy] && souy[is] < model.yCen[iy + 1])
			{
				sou_index[1][is] = iy;
				if (souy[is] == model.yCen[iy])ub = true;//upbound
				else
				{
					ub = false; db = false;
				}
				break;
			}
			else if (souy[is] == model.yCen[model.Num_yCen - 1])
			{
				sou_index[1][is] = model.Num_yCen - 1;
				db = true;//downbound
				break;
			}
			else
			{
				iy++;
			}
		}
		if (model.Land && abs(model.Vel[sou_index[0][is]][sou_index[1][is]] - 340.) < 10.)
		{
	//		cerr << "Source higer than the topograpgy!" << endl;
			if (lb || rb)
			{
				sou_index[1][is] = model.topo_ynum[sou_index[0][is]];
				souy[is] = model.xTopo[sou_index[0][is]];
			}
			else
			{
				sou_index[1][is] = max(model.topo_ynum[sou_index[0][is]], model.topo_ynum[sou_index[0][is] + 1]);
				souy[is] = max(model.xTopo[sou_index[0][is]], model.xTopo[sou_index[0][is] + 1]);
			}
			ix = 0; lb = false; rb = false; ub = false; db = false;
			while (1)
			{
				if (ix > model.Num_xCen - 1)
				{
					cerr << "source index failed! check the input data!" << endl;
					break;
				}
				if (soux[is] >= model.xCen[ix] && soux[is] < model.xCen[ix + 1])
				{
					sou_index[0][is] = ix;
					if (soux[is] == model.xCen[ix])lb = true;//leftbound
					else
					{
						lb = false; rb = false;
					}
					break;
				}
				else if (soux[is] == model.xCen[model.Num_xCen - 1])
				{
					sou_index[0][is] = model.Num_xCen - 1;
					rb = true;//rightbound
					break;
				}
				else
				{
					ix++;
				}
			}
			//yindex search
			iy = 0;
			while (1)
			{
				if (iy > model.Num_yCen - 1)
				{
					cerr << "source index failed! check the input data!" << endl;
					break;
				}
				if (souy[is] >= model.yCen[iy] && souy[is] < model.yCen[iy + 1])
				{
					sou_index[1][is] = iy;
					if (souy[is] == model.yCen[iy])ub = true;//upbound
					else
					{
						ub = false; db = false;
					}
					break;
				}
				else if (souy[is] == model.yCen[model.Num_yCen - 1])
				{
					sou_index[1][is] = model.Num_yCen - 1;
					db = true;//downbound
					break;
				}
				else
				{
					iy++;
				}
			}
		}
		else if (!model.Land && abs(model.Vel[sou_index[0][is]][sou_index[1][is]] - 1500.) < 100.)
		{
		//	cerr << "Source higer than the topograpgy!" << endl;
			if (lb || rb)
			{
				sou_index[1][is] = model.topo_ynum[sou_index[0][is]];
				souy[is] = model.xTopo[sou_index[0][is]];
			}
			else
			{
				sou_index[1][is] = max(model.topo_ynum[sou_index[0][is]], model.topo_ynum[sou_index[0][is] + 1]);
				souy[is] = max(model.xTopo[sou_index[0][is]], model.xTopo[sou_index[0][is] + 1]);
			}
			ix = 0; lb = false; rb = false; ub = false; db = false;
			while (1)
			{
				if (ix > model.Num_xCen - 1)
				{
					cerr << "source xindex failed! check the input data!" << endl;
					exit(1);
					break;
				}
				if (soux[is] >= model.xCen[ix] && soux[is] < model.xCen[ix + 1])
				{
					sou_index[0][is] = ix;
					if (soux[is] == model.xCen[ix])lb = true;//leftbound
					else
					{
						lb = false; rb = false;
					}
					break;
				}
				else if (soux[is] == model.xCen[model.Num_xCen - 1])
				{
					sou_index[0][is] = model.Num_xCen - 1;
					rb = true;//rightbound
					break;
				}
				else
				{
					ix++;
				}
			}
			//yindex search
			iy = 0;
			while (1)
			{
				if (iy > model.Num_yCen - 1)
				{
					cerr << "source yindex failed! check the input data!" << endl;
					exit(1);
					break;
				}
				if (souy[is] >= model.yCen[iy] && souy[is] < model.yCen[iy + 1])
				{
					sou_index[1][is] = iy;
					if (souy[is] == model.yCen[iy])ub = true;//upbound
					else
					{
						ub = false; db = false;
					}
					break;
				}
				else if (souy[is] == model.yCen[model.Num_yCen - 1])
				{
					sou_index[1][is] = model.Num_yCen - 1;
					db = true;//downbound
					break;
				}
				else
				{
					iy++;
				}
			}
		}
		//decide type:
		//0:inbound	
		//1:leftbound		5:leftup point
		//2:downbound	6:leftdown point
		//3:rightbound		7.rightdown point
		//4:upbound			8.rightup point
		if (!lb && !ub && !rb && !db)
		{
			sou_type[is] = 0;
			continue;
		}
		if (lb && !ub && !db)
		{
			sou_type[is] = 1;
			continue;
		}
		if (db && !lb && !rb)
		{
			sou_type[is] = 2;
			continue;
		}
		if (rb && !ub && !db)
		{
			sou_type[is] = 3;
			continue;
		}
		if (ub && !lb && !rb)
		{
			sou_type[is] = 4;
			continue;
		}
		if (lb && ub)
		{
			sou_type[is] = 5;
			continue;
		}
		if (lb && db)
		{
			sou_type[is] = 6;
			continue;
		}
		if (rb && db)
		{
			sou_type[is] = 7;
			continue;
		}
		if (rb && ub)
		{
			sou_type[is] = 8;
			continue;
		}
	}
}

void TomoSrs2D::_Rec_index(TomoMesh2D model)
{
	int x_index, y_index, ix, iy;
	bool lb, ub, rb, db;
	for (int ishot = 0; ishot < Num_source; ishot++)
	{
		data[ishot].rec_refr_index.resize(2);
		for (int i = 0; i < 2; i++)data[ishot].rec_refr_index[i].resize(data[ishot].nrec_refrac);
		data[ishot].rec_refl_index.resize(2);
		for (int i = 0; i < 2; i++)data[ishot].rec_refl_index[i].resize(data[ishot].nrec_reflec);
		data[ishot].rec_refr_type.resize(data[ishot].nrec_refrac);
		data[ishot].rec_refl_type.resize(data[ishot].nrec_reflec);
		//refraction
		for (int irec = 0; irec < data[ishot].nrec_refrac; irec++)
		{
			if (data[ishot].recy_refrac[irec] < model.yMin)data[ishot].recy_refrac[irec] = model.yMin;
			if (data[ishot].recy_refrac[irec] > model.yMax)data[ishot].recy_refrac[irec] = model.yMax;
			ix = 0; lb = false; rb = false; ub = false; db = false;
			while (1)
			{
				if (ix > model.Num_xCen - 1)
				{
					cerr << "Receiver xindex failed! check the input data!" << endl;
					exit(1);
					break;
				}
				if (data[ishot].recx_refrac[irec] >= model.xCen[ix] && data[ishot].recx_refrac[irec] < model.xCen[ix + 1])
				{
					data[ishot].rec_refr_index[0][irec] = ix;
					if (data[ishot].recx_refrac[irec] == model.xCen[ix])lb = true;//leftbound
					else
					{
						lb = false; rb = false;
					}
					break;
				}
				else if (data[ishot].recx_refrac[irec] == model.xCen[model.Num_xCen - 1])
				{
					data[ishot].rec_refr_index[0][irec] = model.Num_xCen - 1;
					rb = true;//rightbound
					break;
				}
				else
				{
					ix++;
				}
			}
			//yindex search
			iy = 0;
			while (1)
			{
				if (iy > model.Num_yCen - 1)
				{
					cerr << "Receiver yindex failed! check the input data!" << endl;
					exit(1);
					break;
				}
				if (data[ishot].recy_refrac[irec] >= model.yCen[iy] && data[ishot].recy_refrac[irec] < model.yCen[iy + 1])
				{
					data[ishot].rec_refr_index[1][irec] = iy;
					if (data[ishot].recy_refrac[irec] == model.yCen[iy])ub = true;//upbound
					else
					{
						ub = false; db = false;
					}
					break;
				}
				else if (data[ishot].recy_refrac[irec] == model.yCen[model.Num_yCen - 1])
				{
					data[ishot].rec_refr_index[1][irec] = model.Num_yCen - 1;
					db = true;//downbound
					break;
				}
				else
				{
					iy++;
				}
			}
			if (model.Land && abs(model.Vel[data[ishot].rec_refr_index[0][irec]][data[ishot].rec_refr_index[1][irec]] - 340.) < 10.)
			{
		//		cerr << "Receiver higer than the topograpgy!" << endl;
				if (lb || rb)
				{
					data[ishot].rec_refr_index[1][irec] = model.topo_ynum[data[ishot].rec_refr_index[0][irec]];
					data[ishot].recy_refrac[irec] = model.xTopo[data[ishot].rec_refr_index[0][irec]];
				}
				else
				{
					data[ishot].rec_refr_index[1][irec] = max(model.topo_ynum[data[ishot].rec_refr_index[0][irec]], model.topo_ynum[data[ishot].rec_refr_index[0][irec] + 1]);
					data[ishot].recy_refrac[irec] = max(model.xTopo[data[ishot].rec_refr_index[0][irec]], model.xTopo[data[ishot].rec_refr_index[0][irec] + 1]);
				}
				ix = 0; lb = false; rb = false; ub = false; db = false;
				while (1)
				{
					if (ix > model.Num_xCen - 1)
					{
						cerr << "Receiver xindex failed! check the input data!" << endl;
						break;
					}
					if (data[ishot].recx_refrac[irec] >= model.xCen[ix] && data[ishot].recx_refrac[irec] < model.xCen[ix + 1])
					{
						data[ishot].rec_refr_index[0][irec] = ix;
						if (data[ishot].recx_refrac[irec] == model.xCen[ix])lb = true;//leftbound
						else
						{
							lb = false; rb = false;
						}
						break;
					}
					else if (data[ishot].recx_refrac[irec] == model.xCen[model.Num_xCen - 1])
					{
						data[ishot].rec_refr_index[0][irec] = model.Num_xCen - 1;
						rb = true;//rightbound
						break;
					}
					else
					{
						ix++;
					}
				}
				//yindex search
				iy = 0;
				while (1)
				{
					if (iy > model.Num_yCen - 1)
					{
						cerr << data[ishot].recy_refrac[irec] << " " << model.yCen[0] << endl;
						cerr << "Receiver yindex failed! check the input data!" << endl;
						break;
					}
					if (data[ishot].recy_refrac[irec] >= model.yCen[iy] && data[ishot].recy_refrac[irec] < model.yCen[iy + 1])
					{
						data[ishot].rec_refr_index[1][irec] = iy;
						if (data[ishot].recy_refrac[irec] == model.yCen[iy])ub = true;//upbound
						else
						{
							ub = false; db = false;
						}
						break;
					}
					else if (data[ishot].recy_refrac[irec] == model.yCen[model.Num_yCen - 1])
					{
						data[ishot].rec_refr_index[1][irec] = model.Num_yCen - 1;
						db = true;//downbound
						break;
					}
					else
					{
						iy++;
					}
				}
			}
			//decide type:
			//0:inbound	
			//1:leftbound		5:leftup point
			//2:downbound	6:leftdown point
			//3:rightbound		7.rightdown point
			//4:upbound			8.rightup point
			if (!lb && !ub && !rb && !db)
			{
				data[ishot].rec_refr_type[irec] = 0;
				continue;
			}
			if (lb && !ub && !db)
			{
				data[ishot].rec_refr_type[irec] = 1;
				continue;
			}
			if (db && !lb && !rb)
			{
				data[ishot].rec_refr_type[irec] = 2;
				continue;
			}
			if (rb && !ub && !db)
			{
				data[ishot].rec_refr_type[irec] = 3;
				continue;
			}
			if (ub && !lb && !rb)
			{
				data[ishot].rec_refr_type[irec] = 4;
				continue;
			}
			if (lb && ub)
			{
				data[ishot].rec_refr_type[irec] = 5;
				continue;
			}
			if (lb && db)
			{
				data[ishot].rec_refr_type[irec] = 6;
				continue;
			}
			if (rb && db)
			{
				data[ishot].rec_refr_type[irec] = 7;
				continue;
			}
			if (rb && ub)
			{
				data[ishot].rec_refr_type[irec] = 8;
				continue;
			}
		}
		//reflection
		for (int irec = 0; irec < data[ishot].nrec_reflec; irec++)
		{
			if (data[ishot].recy_reflec[irec] < model.yMin)data[ishot].recy_reflec[irec] = model.yMin;
			if (data[ishot].recy_reflec[irec] > model.yMax)data[ishot].recy_reflec[irec] = model.yMax;
			ix = 0; lb = false; rb = false; ub = false; db = false;
			while (1)
			{
				if (ix > model.Num_xCen - 1)
				{
					cerr << "Receiver xindex failed! check the input data!" << endl;
					exit(1);
					break;
				}
				if (data[ishot].recx_reflec[irec] >= model.xCen[ix] && data[ishot].recx_reflec[irec] < model.xCen[ix + 1])
				{
					data[ishot].rec_refl_index[0][irec] = ix;
					if (data[ishot].recx_reflec[irec] == model.xCen[ix])lb = true;//leftbound
					else
					{
						lb = false; rb = false;
					}
					break;
				}
				else if (data[ishot].recx_reflec[irec] == model.xCen[model.Num_xCen - 1])
				{
					data[ishot].rec_refl_index[0][irec] = model.Num_xCen - 1;
					rb = true;//rightbound
					break;
				}
				else
				{
					ix++;
				}
			}
			//yindex search
			iy = 0;
			while (1)
			{
				if (iy > model.Num_yCen - 1)
				{
					cerr << "Receiver yindex failed! check the input data!" << endl;
					exit(1);
					break;
				}
				if (data[ishot].recy_reflec[irec] >= model.yCen[iy] && data[ishot].recy_reflec[irec] < model.yCen[iy + 1])
				{
					data[ishot].rec_refl_index[1][irec] = iy;
					if (data[ishot].recy_reflec[irec] == model.yCen[iy])ub = true;//upbound
					else
					{
						ub = false; db = false;
					}
					break;
				}
				else if (data[ishot].recy_reflec[irec] == model.yCen[model.Num_yCen - 1])
				{
					data[ishot].rec_refl_index[1][irec] = model.Num_yCen - 1;
					db = true;//downbound
					break;
				}
				else
				{
					iy++;
				}
			}
			if (model.Land && abs(model.Vel[data[ishot].rec_refl_index[0][irec]][data[ishot].rec_refl_index[1][irec]] - 340.) < 10.)
			{
		//		cerr << "Receiver higer than the topograpgy!" << endl;
				if (lb || rb)
				{
					data[ishot].rec_refl_index[1][irec] = model.topo_ynum[data[ishot].rec_refl_index[0][irec]];
					data[ishot].recy_reflec[irec] = model.xTopo[data[ishot].rec_refl_index[0][irec]];
				}
				else
				{
					data[ishot].rec_refl_index[1][irec] = max(model.topo_ynum[data[ishot].rec_refl_index[0][irec]], model.topo_ynum[data[ishot].rec_refl_index[0][irec] + 1]);
					data[ishot].recy_reflec[irec] = max(model.xTopo[data[ishot].rec_refl_index[0][irec]], model.xTopo[data[ishot].rec_refl_index[0][irec] + 1]);
				}
				ix = 0; lb = false; rb = false; ub = false; db = false;
				while (1)
				{
					if (ix > model.Num_xCen - 1)
					{
						cerr << "Receiver xindex failed! check the input data!" << endl;
						exit(1);
						break;
					}
					if (data[ishot].recx_reflec[irec] >= model.xCen[ix] && data[ishot].recx_reflec[irec] < model.xCen[ix + 1])
					{
						data[ishot].rec_refl_index[0][irec] = ix;
						if (data[ishot].recx_reflec[irec] == model.xCen[ix])lb = true;//leftbound
						else
						{
							lb = false; rb = false;
						}
						break;
					}
					else if (data[ishot].recx_reflec[irec] == model.xCen[model.Num_xCen - 1])
					{
						data[ishot].rec_refl_index[0][irec] = model.Num_xCen - 1;
						rb = true;//rightbound
						break;
					}
					else
					{
						ix++;
					}
				}
				//yindex search
				iy = 0;
				while (1)
				{
					if (iy > model.Num_yCen - 1)
					{
						cerr << "Receiver yindex failed! check the input data!" << endl;
						exit(1);
						break;
					}
					if (data[ishot].recy_reflec[irec] >= model.yCen[iy] && data[ishot].recy_reflec[irec] < model.yCen[iy + 1])
					{
						data[ishot].rec_refl_index[1][irec] = iy;
						if (data[ishot].recy_reflec[irec] == model.yCen[iy])ub = true;//upbound
						else
						{
							ub = false; db = false;
						}
						break;
					}
					else if (data[ishot].recy_reflec[irec] == model.yCen[model.Num_yCen - 1])
					{
						data[ishot].rec_refl_index[1][irec] = model.Num_yCen - 1;
						db = true;//downbound
						break;
					}
					else
					{
						iy++;
					}
				}
			}
				//decide type:
				//0:inbound	
				//1:leftbound		5:leftup point
				//2:downbound	6:leftdown point
				//3:rightbound		7.rightdown point
				//4:upbound			8.rightup point
				if (!lb && !ub && !rb && !db)
				{
					data[ishot].rec_refl_type[irec] = 0;
					continue;
				}
				if (lb && !ub && !db)
				{
					data[ishot].rec_refl_type[irec] = 1;
					continue;
				}
				if (db && !lb && !rb)
				{
					data[ishot].rec_refl_type[irec] = 2;
					continue;
				}
				if (rb && !ub && !db)
				{
					data[ishot].rec_refl_type[irec] = 3;
					continue;
				}
				if (ub && !lb && !rb)
				{
					data[ishot].rec_refl_type[irec] = 4;
					continue;
				}
				if (lb && ub)
				{
					data[ishot].rec_refl_type[irec] = 5;
					continue;
				}
				if (lb && db)
				{
					data[ishot].rec_refl_type[irec] = 6;
					continue;
				}
				if (rb && db)
				{
					data[ishot].rec_refl_type[irec] = 7;
					continue;
				}
				if (rb && ub)
				{
					data[ishot].rec_refl_type[irec] = 8;
					continue;
				}
			}
		}
	}

