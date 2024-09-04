#pragma once
#ifndef TOMO_SRS2D_H_
#define TOMO_SRS2D_H_
#include "CommonHeaders.h"
#include "Tomo_mesh2D.h"

class TomoSrs2D
{
public:
	TomoSrs2D() {
	}

	~TomoSrs2D() {
	}
	struct rec_data
	{
		int nrec_refrac;
		int nrec_reflec;
		vector<double> recx_refrac, recy_refrac;
		vector<double> recx_reflec, recy_reflec;
		vector<double> time_refrac, time_reflec, err_refrac, err_reflec;
		vector< vector<int> > rec_refr_index,  rec_refl_index;
		vector<int> rec_refr_type, rec_refl_type;
	};
	//Basic Mesh
	int Num_source, Num_rec_refr, Num_rec_refl;//source number
	rec_data* data;
	vector<double> soux, souy;
	vector<int> sou_type;
	vector< vector<int> >sou_index;
	void Read_srsData(string filename, TomoMesh2D model);

protected:

	void _Source_index(TomoMesh2D model);
	void _Rec_index(TomoMesh2D model);

};

#endif
