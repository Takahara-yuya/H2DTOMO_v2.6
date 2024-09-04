#ifndef TOOL_H_
#define TOOL_H_
#include "CommonHeaders.h"
inline int* alloc_int_1d(int n1)
{
	int* a = new int[n1];
	return a;
}
inline void zero_1d(int* a, int n1)
{
	for (int i = 0; i < n1; i++)
		a[i] = 0;
}
inline void free_1d(int* a)
{
	delete[] a;
}
inline double* alloc_double_1d(int n1)
{
	double* a = new double[n1];
	return a;
}
inline void zero_1d(double* a, int n1)
{
	for (int i = 0; i < n1; i++)
		a[i] = 0.0;
}
inline void free_1d(double* a)
{
	delete[] a;
}
inline void free_1i(int* a)
{
	delete[] a;
}
inline int** alloc_int_2d(int n1, int n2)
{
	int i;
	int** a;
	a = new int* [n1];
	for (i = 0; i < n1; i++)
		a[i] = new int[n2];
	return a;
}
inline void free_2i(int** a, int n1)
{
	for (int i = 0; i < n1; i++)
		delete[] a[i];
	delete[] a;
}
inline void zero_2i(int** a, int n1, int n2)
{
	int i, j;
	for (i = 0; i < n1; i++)
		for (j = 0; j < n2; j++)
			a[i][j] = 0.0;
}
inline double** alloc_double_2d(int n1, int n2)
{
	int i;
	double** a;
	a = new double* [n1];
	for (i = 0; i < n1; i++)
		a[i] = new double[n2];
	return a;
}
inline void zero_2d(double** a, int n1, int n2)
{
	int i, j;
	for (i = 0; i < n1; i++)
		for (j = 0; j < n2; j++)
			a[i][j] = 0.0;
}
inline void free_2d(double** a, int n1)
{
	for (int i = 0; i < n1; i++)
		delete[] a[i];
	delete[] a;
}
inline void write_double_1d_bin(char* file, double* a, int n1)
{
	ofstream ofile(file);
	if (!ofile)
	{
		cerr << "Can't open the output 1d bin file:" << file << endl;
		exit(0);
	}
	for (int i = 0; i < n1; i++)
		ofile.write((char*)&a[i], sizeof(double));
	ofile.close();
	ofile.clear();
}
inline void write_double_1d_asi(char* file, double* a, int n1, double dpre, double dx)
{
	ofstream ofile(file);
	if (!ofile)
	{
		cerr << "Can't open the output 1d asc file:" << file << endl;
		exit(0);
	}
	ofile << n1 << endl;
	for (int i = 0; i < n1; i++)
		ofile << (double)i * dx / dpre << "\t" << a[i] << endl;
	ofile.close();
	ofile.clear();
}
inline void write_double_1d_2d_bin(char* file, double* a, int n1, int n2)
{
	ofstream ofile(file);
	if (!ofile)
	{
		cerr << "Can't open the output 1d-2d bin file" << file << endl;
		exit(0);
	}
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
			ofile.write((char*)&a[i * n2 + j], sizeof(double));
	ofile.close();
	ofile.clear();
}
inline void write_double_1d_2d_asi(char* file, double* a, int n1, int n2, double d_pre, double dx, double dz)
{
	ofstream ofile(file);
	if (!ofile)
	{
		cerr << "Can't open the output 1d-2d asc file" << file << endl;
		exit(0);
	}
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
			ofile << (double)i * dx / d_pre << " " << (double)-fabs(j * dz / d_pre) << " " << a[i * n2 + j] << endl;
	ofile.close();
	ofile.clear();
}
inline void write_double_2d_bin(string file, vector< vector<double> >& a, int n1, int n2)
{
	ofstream ofile(file);
	if (!ofile)
	{
		cerr << "Can't open the output 2d bin file" << file << endl;
		exit(0);
	}
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
			ofile.write((char*)&a[i][j], sizeof(double));
	ofile.close();
	ofile.clear();
}
inline void write_allocdouble_2d_bin(string file, double** a, int n1, int n2)
{
	ofstream ofile(file);
	if (!ofile)
	{
		cerr << "Can't open the output 2d bin file" << file << endl;
		exit(0);
	}
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
			ofile.write((char*)&a[i][j], sizeof(double));
	ofile.close();
	ofile.clear();
}
inline void read_1d_asc(string file, double* a, int nx)
{
	ifstream ifdata(file);
	if (!ifdata)
	{
		cerr << "Can't open the input 1d asc file:" << file << endl;
		exit(0);
	}
	for (int ix = 0; ix < nx; ix++)
		ifdata >> a[ix];
	ifdata.close();
	ifdata.clear();
}
inline void read_1d_bin(string file, double* a, int nx)
{
	ifstream ifdata(file);
	if (!ifdata)
	{
		cerr << "Can't open the input 1d bin file" << file << endl;
		exit(0);
	}
	for (int ix = 0; ix < nx; ix++)
		ifdata.read((char*)&a[ix], sizeof(double));
	ifdata.close();
	ifdata.clear();
}
inline void read_2d_bin(string file, vector< vector<double> >& a, int nx, int nz)
{
	const char* charArray = file.c_str();
	double** a1;
	a1 = new double* [nx];
	for (int i = 0; i < nx; i++)
		a1[i] = new double[nz];
	FILE* f2;
	f2 = fopen(charArray, "rb");
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < nz; j++)
		{
			size_t sizeRead=fread(&a1[i][j], sizeof(double), 1, f2);
			if (sizeRead != 1) {
				printf("error!\n");
			}
		}
	fclose(f2);
	for (int ix = 0; ix < nx; ix++)
	{
		for (int iz = 0; iz < nz; iz++)
		{
			a[ix][iz] = a1[ix][iz];
		}
	}
	for (int i = 0; i < nx; i++)
		delete[] a1[i];

}
inline double dist(double x1, double x2, double y1, double y2)
{
	double dis;
	dis = pow(x1 - x2, 2) + pow(y1 - y2, 2);
	dis = sqrt(dis);
	return dis;
}
inline void MakeGauss(double** a, int h_size, double  sigma)
{
	double siz;
	int i, j;
	siz = (h_size - 1.0) / 2.0;
	double** b;
	b = alloc_double_2d(h_size, h_size);
	for (i = 0; i < h_size; i++)
		for (j = 0; j < h_size; j++)
		{
			a[i][j] = -siz + j;
			b[j][i] = a[i][j];
		}
	double h_sum = 0.0;
	for (i = 0; i < h_size; i++)
		for (j = 0; j < h_size; j++)
		{
			a[i][j] = a[i][j] * a[i][j];
			b[i][j] = b[i][j] * b[i][j];
			a[i][j] = -(a[i][j] + b[i][j]) / (2 * sigma * sigma);
			a[i][j] = exp(a[i][j]);
			h_sum += a[i][j];
		}
	free_2d(b, h_size);
	for (i = 0; i < h_size; i++)
		for (j = 0; j < h_size; j++)
			a[i][j] = a[i][j] / h_sum;
}
inline void exmodel(double** mod, double** exmod, int nx, int nz, int bond)
{
	int ix, iz;
	for (ix = 0; ix < bond; ix++)
		for (iz = 0; iz < bond; iz++)
		{
			exmod[ix][iz] = mod[0][0];
			exmod[bond + nx + ix][iz] = mod[nx - 1][0];
			exmod[ix][bond + nz + iz] = mod[0][nz - 1];
			exmod[bond + nx + ix][bond + nz + iz] = mod[nx - 1][nz - 1];
		}
	for (ix = 0; ix < nx; ix++)
		for (iz = 0; iz < bond; iz++)
		{
			exmod[ix + bond][iz] = mod[ix][0];
			exmod[ix + bond][bond + nz + iz] = mod[ix][nz - 1];
		}
	for (ix = 0; ix < bond; ix++)
		for (iz = 0; iz < nz; iz++)
		{
			exmod[ix][iz + bond] = mod[0][iz];
			exmod[bond + nx + ix][iz + bond] = mod[nx - 1][iz];
		}
	for (ix = 0; ix < nx; ix++)
		for (iz = 0; iz < nz; iz++)
			exmod[ix + bond][iz + bond] = mod[ix][iz];
}
inline void normolization(double** grad, int nx, int ny)
{
	double min, max;
	min = 99999999999999;
	max = -99999999999999;
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			if (abs(grad[i][j]) > max)max = abs(grad[i][j]);
			if (abs(grad[i][j]) < min)min = abs(grad[i][j]);
		}
	}
	if (max - min > 0)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				if (grad[i][j] != 0)
				{
					grad[i][j] = (abs(grad[i][j]) / grad[i][j]) * (abs(grad[i][j]) - min) / (max - min);
				}
				else
				{
					grad[i][j] = 0;
				}
			}
		}
	}
}

#endif
