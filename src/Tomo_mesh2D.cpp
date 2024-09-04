/*
***********************************************************************

Tomo_mesh2D.cpp	(Constructing Tomo mesh)
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

#include "../include/Tomo_mesh2D.h"

#define v_air 0.34
#define v_water 1.5

void TomoMesh2D::Init_model(string filename)
{
	string tmp, uniform, gridfile, velocityfile, topofile, model_type, data_type, interfacefile;
	double dx, dy;
	ifstream fin(filename);
	if (!fin)
	{
		cerr << "Error reading meshsetting file!" << endl;
		exit(1);
	}
	fin >> tmp >> tmp >> data_type;
	fin >> tmp >> tmp >> Num_xCen >> Num_yCen;
	fin >> tmp >> tmp >> uniform;
	fin >> tmp >> tmp >> dx >> dy;
	fin >> tmp >> tmp >> gridfile;
	fin >> tmp >> tmp >> topofile;
	fin >> tmp >> tmp >> xMin >> yMin;
	fin >> tmp >> tmp >> model_type;
	fin >> tmp >> tmp >> filetype;
	fin >> tmp >> tmp >> velocityfile;
	fin >> tmp >> tmp >> vstart >> dv;
	fin >> tmp >> tmp >> interfacefile;
	fin.close();
	//
	dx = dx * 1000; dy = dy * 1000;
	xMin = xMin * 1000; yMin = yMin * 1000.;
	//read finished! construct model
	//string change:
	for (char& c : model_type) {
		c = std::tolower(c);
	}
	for (char& c : uniform) {
		c = std::tolower(c);
	}
	for (char& c : data_type) {
		c = std::tolower(c);
	}
	for (char& c : filetype) {
		c = std::tolower(c);
	}
	//data_type
	if (data_type == "land")Land = true;
	else Land = false;
	//allocated parameter
	Num_xNode = Num_xCen + 1; Num_yNode = Num_yCen + 1;
	Num_Cen = Num_xCen * Num_yCen; Num_Node = Num_xNode * Num_yNode;
	dX.resize(Num_xCen); dY.resize(Num_yCen);
	xCen.resize(Num_xCen); yCen.resize(Num_yCen);
	xTopo.resize(Num_xCen); interface.resize(Num_xCen);
	topo_ynum.resize(Num_xCen); intf_ynum.resize(Num_xCen);
	xNode.resize(Num_xNode); yNode.resize(Num_yNode);
	//
	Vel.resize(Num_xCen); Slowness.resize(Num_xCen);
	for (int i = 0; i < Vel.size(); i++)
	{
		Vel[i].resize(Num_yCen);
		Slowness[i].resize(Num_yCen);
	}
	Vel1D.resize(Num_Cen); Slowness1D.resize(Num_Cen);
	//decide dx,dy,xCen,yCen
	if (uniform == "yes")
	{
		for (int i = 0; i < Num_xCen; i++)dX[i] = dx;
		for (int i = 0; i < Num_yCen; i++)dY[i] = dx;
		for (int i = 0; i < Num_xCen; i++)xCen[i] = xMin + dx * i;
		for (int i = 0; i < Num_yCen; i++)yCen[i] = yMin + dy * i;
		for (int i = 0; i < Num_xNode; i++)xNode[i] = xMin - 0.5 * dx + dx * i;
		for (int i = 0; i < Num_yNode; i++)yNode[i] = yMin - 0.5 * dy + dy * i;
	}
	else
	{
		_Read_Gridfile(gridfile);
	}
	//model make
	if (model_type == "gradient1")
	{	//read topo
		_Read_Topofile(topofile, model_type);
		_Make1Dmodel(1);
	}
	else if (model_type == "gradient0")
	{
		_Read_Topofile(topofile, model_type);
		_Make1Dmodel(0);
	}
	else
	{
		_Read_Topofile(topofile, model_type);
		_Read_Velocityfile(velocityfile, filetype);
	}
	//make slowness
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		for (int iy = 0; iy < Num_yCen; iy++)
		{
			int index = iy + ix * Num_yCen;
			Slowness[ix][iy] = 1. / Vel[ix][iy];
			Slowness1D.coeffRef(index) = Slowness[ix][iy];
			Vel1D.coeffRef(index) = Vel[ix][iy];
		}
	}
	//interface
	_Read_Interfacefile(interfacefile);
}

void TomoMesh2D::Outputmesh(string filename, string binary)
{
	if (binary == "binary")
	{
		write_double_2d_bin(filename, Vel, Num_xCen, Num_yCen);
	}
	else
	{
		ofstream fout(filename);
		if (!fout)
		{
			cerr << "Error opening output file!" << endl;
			exit(1);
		}
		for (int ix = 0; ix < Num_xCen; ix++)
			for (int iy = 0; iy < Num_yCen; iy++)
			{
				fout << xCen[ix] / 1000. << " " << yCen[iy] / 1000. << " " << Vel[ix][iy] << endl;
			}
		fout.close();
	}
}

void TomoMesh2D::_Make1Dmodel(int type)
{
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		for (int iy = 0; iy < topo_ynum[ix]; iy++)
		{
			if (Land)Vel[ix][iy] = v_air * 1000;
			else Vel[ix][iy] = v_water * 1000;
		}
		for (int iy = topo_ynum[ix]; iy < Num_yCen; iy++)
		{
			Vel[ix][iy] = vstart + dv * (iy - topo_ynum[ix] * double(type));
			Vel[ix][iy] = Vel[ix][iy] * 1000.;
		}
	}
}

void TomoMesh2D::_Read_Gridfile(string filename)
{
	string tmp;
}

void TomoMesh2D::_Read_Velocityfile(string filename, string binary)
{
	double tmp1, tmp2;
	if (binary == "binary")
	{
		read_2d_bin(filename, Vel, Num_xCen, Num_yCen);
		for (int ix = 0; ix < Num_xCen; ix++)
			for (int iy = 0; iy < Num_yCen; iy++)
			{
				Vel[ix][iy] *= 1000.;
			}
	}
	else
	{
		ifstream fin(filename);
		if (!fin)
		{
			cerr << "Error reading Velocity file!" << endl;
			exit(1);
		}
		for (int ix = 0; ix < Num_xCen; ix++)
			for (int iy = 0; iy < Num_yCen; iy++)
			{
				fin >> tmp1 >> tmp2 >> Vel[ix][iy];
			}
		fin.close();
	}
}

void TomoMesh2D::_Read_Topofile(string filename, string model_type)
{
	double tmp1, tmp2;
	vector<double> topo_tmp, x_tmp;
	ifstream fin(filename);
	if (!fin)
	{
		cerr << "Error reading Topo file!" << endl;
		exit(1);
	}
	while (!fin.eof())
	{
		fin >> tmp1 >> tmp2;
		//if (!Land && tmp2 < 0)
		//{
		////	cerr << "Marine topo should larger than zero!" << endl;
		////	exit(1);
		//	tmp2 = 0;
		//}
		x_tmp.push_back(tmp1 * 1000.); topo_tmp.push_back(tmp2 * 1000.);
	}
	fin.close();
	//interplote
	Interpolation inter;
	inter.linearinter1D_init(x_tmp, topo_tmp);
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		xTopo[ix] = inter.linearinter1D(xCen[ix]);
	}
	//adjust yMin
//	if (Land)
//	{
	auto minElement = std::min_element(xTopo.begin(), xTopo.end());
	yMin = *minElement;
	for (int iy = 0; iy < Num_yCen; iy++)
	{
		yCen[iy] += yMin;
	}
	//define yMax xMax
	yMax = -999999.; xMax = -999999.;
	for (int iy = 0; iy < Num_yCen; iy++)
	{
		if (yMax < yCen[iy])yMax = yCen[iy];
	}
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		if (xMax < xCen[ix])xMax = xCen[ix];
	}
	//find index
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		int iy = 0;
		while (1)
		{
			if (iy > Num_yCen - 2)
			{
				cerr << "Topo index search failed!" << endl;
				cerr << xTopo[ix] << endl;
			}
			if (xTopo[ix] <= yCen[0])
			{
				topo_ynum[ix] = 0;
				break;
			}
			else
			{
				if (yCen[iy] < xTopo[ix] && xTopo[ix] <= yCen[iy + 1])
				{
					topo_ynum[ix] = iy + 1;
					break;
				}
				else
				{
					iy++;
				}
			}
		}
	}
}

void TomoMesh2D::_Read_Interfacefile(string filename)
{
	double tmp1, tmp2;
	vector<double> interface_tmp, x_tmp;
	ifstream fin(filename);
	if (!fin)
	{
		cerr << "Error reading Interface file!" << endl;
		exit(1);
	}
	while (!fin.eof())
	{
		fin >> tmp1 >> tmp2;
		if (tmp2 < 0)
		{
			cerr << "Interface should be subsurface!" << endl;
			tmp2 = fabs(tmp2);
		}
		x_tmp.push_back(tmp1 * 1000.); interface_tmp.push_back(tmp2 * 1000.);
	}
	fin.close();
	//interplote
	Interpolation inter1;
	inter1.linearinter1D_init(x_tmp, interface_tmp);
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		interface[ix] = inter1.linearinter1D(xCen[ix]);
		if (interface[ix] > yMax)interface[ix] = yMax;
	}
	//find index
	for (int ix = 0; ix < Num_xCen; ix++)
	{
		int iy = 0;
		while (1)
		{
			if (iy > Num_yCen - 2)
			{
				cerr << "Interface index search failed!" << endl;
				exit(1);
			}
			if (interface[ix] <= yCen[0])
			{
				intf_ynum[ix] = 0;
				break;
			}
			else
			{
				if (yCen[iy] < interface[ix] && interface[ix] <= yCen[iy + 1])
				{
					intf_ynum[ix] = iy + 1;
					break;
				}
				else
				{
					iy++;
				}
			}
		}
	}
}



