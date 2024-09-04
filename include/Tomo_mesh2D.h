#pragma once
#ifndef TOMO_MESH2D_H_
#define TOMO_MESH2D_H_

#include "CommonHeaders.h"
#include "Interpolation.h"
#include "tool.h"

class TomoMesh2D
{
public:
	//Basic Mesh
	TomoMesh2D() {
	}

	~TomoMesh2D() {
	}
	bool Land;

	int Num_xCen, Num_yCen, Num_Cen;
	int Num_xNode, Num_yNode, Num_Node;

	string filetype;

	double xMin, xMax, yMin, yMax;

	vector<double>  xCen, yCen, xNode, yNode, xTopo, dX, dY, interface;
	vector<int> intf_ynum, topo_ynum;

	vector< vector<double> >  Vel, Slowness;

	VectorXd Vel1D, Slowness1D;

	void Init_model(string filename);

	void Outputmesh(string filename, string binary);

protected:

	double vstart, dv;

	void _Make1Dmodel(int type);

	void _Read_Gridfile(string filename);

	void _Read_Velocityfile(string filename, string binary);

	void _Read_Topofile(string filename, string model_type);

	void _Read_Interfacefile(string filename);

};

#endif