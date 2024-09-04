/*
***********************************************************************

Interpolation.cpp
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

#include "../include/Interpolation.h"

//void  Interpolation::SplineInterpolation1D_init(vector<double>& x, vector<double>& y)
//{
//	if (x.size() != y.size() || x.size() < 3) {
//		std::cerr << "Invalid input data for cubic spline interpolation." << std::endl;
//		return;
//	}
//
//	x_values = x;
//	y_values = y;
//
//	int n = x.size() - 1;
//
//	// Initialize vectors
//	a = y;
//	b.resize(n);
//	c.resize(n + 1, 1.0);
//	d.resize(n);
//
//	// Calculate differences
//	std::vector<double> h(n);
//	std::vector<double> alpha(n);
//	std::vector<double> l(n + 1, 1.0);
//	std::vector<double> mu(n);
//	std::vector<double> z(n + 1);
//
//	for (int i = 0; i < n; ++i) {
//		h[i] = x[i + 1] - x[i];
//		alpha[i] = 3.0 / h[i] * (a[i + 1] - a[i]) - 3.0 / h[i - 1] * (a[i] - a[i - 1]);
//	}
//
//	for (int i = 1; i < n; ++i) {
//		l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
//		mu[i] = h[i] / l[i];
//		z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
//	}
//
//	for (int j = n - 1; j >= 0; --j) {
//		c[j] = z[j] - mu[j] * c[j + 1];
//		b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
//		d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
//	}
//}
//
//double Interpolation::interpolate1D(double x_interp)
//{
//	int n = x_values.size() - 1;
//	int index = 0;
//
//	// Find the interval for interpolation
//	for (int i = 0; i < n; ++i) {
//		if (x_interp >= x_values[i] && x_interp <= x_values[i + 1]) {
//			index = i;
//			break;
//		}
//	}
//
//	// Compute interpolated value using the cubic spline formula
//	double h = x_values[index + 1] - x_values[index];
//	double t = (x_interp - x_values[index]) / h;
//	double result =
//		a[index] + b[index] * t + c[index] * t * t + d[index] * t * t * t;
//
//	return result;
//}

void Interpolation::linearinter1D_init(vector<double>& x_values, vector<double>& y_values)
{
	if (x_values.size() != y_values.size() || x_values.size() < 2) {
		throw invalid_argument("Invalid input data for interpolation");
	}

	x_data = x_values;
	y_data = y_values;
}

// 
double Interpolation::linearinter1D(double x)
{
//	 ȷ      ֵ    ֪   ݷ Χ  
	if (x < x_data.front()) {
		cerr << x << " " << x_data.front() << endl;
		cerr << "Interpolation value is outside the range of known data" << endl;
		exit(1);
	}


	int i = 0;
	while (x > x_data[i + 1])
	{
		++i;
	}

	// 1D
	if (x <= x_data.back())
	{
		double x0 = x_data[i];
		double x1 = x_data[i + 1];
		double y0 = y_data[i];
		double y1 = y_data[i + 1];

		return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
	}
	else
	{
		double x0 = x_data[i - 1];
		double x1 = x_data[i];
		double y0 = y_data[i - 1];
		double y1 = y_data[i];

		return y1 + (y1 - y0) * (x - x1) / (x1 - x0);
	}

}

void Interpolation::Bilinear_init(const std::vector<double>& x_coords,
	const std::vector<double>& y_coords,
	const std::vector< std::vector<double> >& values)
{
	if (x_coords.size() != 2 || y_coords.size() != 2 || values.size() != 2 ||
		values[0].size() != 2 || values[1].size() != 2) {
		std::cerr << "Error: Invalid input dimensions." << std::endl;
		return;
	}

	x_coords_ = x_coords;
	y_coords_ = y_coords;
	values_ = values;
}

double Interpolation::Bilinearinter(double x, double y)
{
	// 
	if (x < x_coords_[0] || x > x_coords_[1] || y < y_coords_[0] || y > y_coords_[1]) {
		std::cerr << "Error: Interpolation point is outside the provided range." << std::endl;
		cerr << x << " " << x_coords_[0] << " " << x_coords_[1] << endl;
		cerr << y << " " << y_coords_[0] << " " << y_coords_[1] << endl;
		exit(1);
	}

	double x_weight = (x - x_coords_[0]) / (x_coords_[1] - x_coords_[0]);
	double y_weight = (y - y_coords_[0]) / (y_coords_[1] - y_coords_[0]);

	double interpolated_value =
		(1 - x_weight) * (1 - y_weight) * values_[0][0] +
		x_weight * (1 - y_weight) * values_[1][0] +
		(1 - x_weight) * y_weight * values_[0][1] +
		x_weight * y_weight * values_[1][1];

	return interpolated_value;
}

