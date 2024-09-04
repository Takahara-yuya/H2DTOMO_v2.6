#pragma once
#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "CommonHeaders.h"


class Interpolation {

public:
	Interpolation() {
	}

	~Interpolation() {
	}
//	void  SplineInterpolation1D_init(vector<double>& x, vector<double>& y);

	void linearinter1D_init(vector<double>& x_values, vector<double>& y_values);

	double linearinter1D(double x_interp);

	void Bilinear_init(const std::vector<double>& x_coords,
		const std::vector<double>& y_coords,
		const std::vector< std::vector<double> >& values);

	double Bilinearinter(double x, double y);

private:
	//linear1D
	std::vector<double> x_data;
	std::vector<double> y_data;
	std::vector<double> z_data;
	std::vector<double> a, b, c, d;
	//Bilinear
	std::vector<double> x_coords_;
	std::vector<double> y_coords_;
	std::vector< std::vector<double> > values_;
};


#endif
